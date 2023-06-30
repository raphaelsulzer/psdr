#include "shape_container.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


Shape_Container::Shape_Container()
{
	SD = nullptr;
}


Shape_Container::Shape_Container(Shape_Detector* _SD)
{
	SD = _SD;
}


int Shape_Container::detect()
{
    SD->_logger->info("Detect planes...");
    SD->_logger->info("epsilon = {}, min_inliers = {}, normal_threshold = {}, knn = {}",
                       SD->get_epsilon(),SD->get_min_points(),SD->get_normal_threshold(),SD->get_knn());

	SD->detect_shapes();
	SD->set_primitives();

	copy_primitives_from_detector();
	copy_support_planes_from_detector();
	SD->clear_primitives();

	discard_degenerate_primitives();

    SD->_logger->info("detected {} planes",SD->get_number_of_non_coplanar_planes());

    return 0;

}

void Shape_Container::regularize()
{
	SD->regularize_shapes();
	SD->set_primitives();

	copy_primitives_from_detector();
	copy_support_planes_from_detector();
	SD->clear_primitives();

	discard_degenerate_primitives();
}

int Shape_Container::refine(int max_iter) {

    SD->_logger->info("Refine planes...");
    SD->refine_shapes(max_iter);
	SD->set_primitives();
	copy_primitives_from_detector();
	copy_support_planes_from_detector();
	SD->clear_primitives();

	discard_degenerate_primitives();

    SD->_logger->info("refined to {} planes",SD->planes_2.size());

    return 0;
}


void Shape_Container::copy_primitives_from_detector()
{
	inexact_alpha_shapes.clear();
	inexact_best_rectangles.clear();
	inexact_convex_hulls.clear();
	convex_hulls_to_alpha_shapes.clear();

	inexact_alpha_shapes = SD->alpha_shapes_pts;
	inexact_best_rectangles = SD->best_rectangles_pts;
	inexact_convex_hulls = SD->convex_hulls_pts;
	convex_hulls_to_alpha_shapes = SD->convex_hulls_seeds;
	
}


void Shape_Container::copy_support_planes_from_detector()
{
	support_planes_of_primitives = SD->planes_2;
}


void Shape_Container::discard_degenerate_primitives()
{
	degenerate_primitives = 0;
	is_not_degenerate = std::vector<bool>(inexact_convex_hulls.size(), true);
	for (size_t i = 0 ; i < inexact_convex_hulls.size() ; ++i) {
		if (inexact_convex_hulls[i].size() < 3) {
			is_not_degenerate[i] = false;
			++degenerate_primitives;
		}
	}
}

bool Shape_Container::is_empty()
{
	return inexact_convex_hulls.empty();
}

int Shape_Container::get_number_of_primitives()
{
	int total_number_of_primitives = int(inexact_convex_hulls.size());
	return (total_number_of_primitives - degenerate_primitives);
	return inexact_convex_hulls.size();
}

int Shape_Container::get_number_of_support_planes()
{
	return SD->get_number_of_non_coplanar_planes();
}

size_t Shape_Container::size_of_inexact_convex_hulls() const
{
	return inexact_convex_hulls.size();
}

size_t Shape_Container::size_of_inexact_best_rectangles() const
{
	return inexact_best_rectangles.size();
}

const std::vector<Inexact_Point_3> & Shape_Container::get_triangle(const int triangle_index) const
{
	return inexact_alpha_shapes[triangle_index];
}

std::pair<int, int> Shape_Container::get_range_of_triangles(const size_t primitive_index) const
{
	return convex_hulls_to_alpha_shapes[primitive_index];
}

const std::vector<Inexact_Point_3> & Shape_Container::get_convex_hull_points(const size_t primitive_index) const
{
	return inexact_convex_hulls[primitive_index];
}

CGAL::Color Shape_Container::get_primitive_color(const size_t primitive_index) const
{
	return SD->planes_to_colors[primitive_index];
}

const std::vector<Inexact_Point_3> & Shape_Container::get_best_rectangle_points(const size_t primitive_index) const
{
	return inexact_best_rectangles[primitive_index];
}

bool Shape_Container::is_primitive_degenerate(const size_t primitive_index) const
{
	return !is_not_degenerate[primitive_index];
}

void Shape_Container::get_inliers(size_t primitive_index, std::vector<Point_with_normal> & inliers)
{
	SD->get_inliers(primitive_index, inliers);
}

void Shape_Container::get_inliers(size_t primitive_index, std::vector<size_t> & inliers)
{
	SD->get_inliers(primitive_index, inliers);
}

size_t Shape_Container::get_total_number_of_inliers()
{
	return SD->get_total_number_of_inliers();
}

const Inexact_Plane & Shape_Container::get_plane(const size_t p) const
{
	return support_planes_of_primitives[p];
}

void Shape_Container::inliers_to_ply(const string & filename)
{
    vector<CGAL::Color> colors;
    vector<Pwn_vector> points;
    Pwn_vector pwn;

    size_t total_convex_hulls = size_of_inexact_convex_hulls();
    pwn.reserve(total_convex_hulls);
    colors.reserve(total_convex_hulls);

    for (size_t p = 0; p < total_convex_hulls; ++p) {
        if (is_primitive_degenerate(p)) continue;

        pwn.clear();
        const CGAL::Color & C = get_primitive_color(p);
        colors.push_back(C);
        get_inliers(p,pwn);
        points.push_back(pwn);

    }

    auto pp = fs::path(filename).parent_path();
    if(!fs::is_directory(pp))
        fs::create_directories(pp);
    std::ofstream stream(filename, std::ios::out);
    if (!stream.is_open()){
        throw std::ios_base::failure("Error : cannot write into an output file");
    }

    double dims[3];
    dims[0] = SD->get_x_max() - SD->get_x_min();
    dims[1] = SD->get_y_max() - SD->get_y_min();
    dims[2] = SD->get_z_max() - SD->get_z_min();

    int lg_max = 0;
    for (int i = 0; i < 3; i++) {
        int lg = (dims[i] == 0 ? 0 : log10(dims[i]));
        if (lg > lg_max) {
            lg_max = lg;
        }
    }

    stream.precision(10 + lg_max);

    size_t nb_pts = 0;
    for (size_t i = 0; i < points.size(); i++) nb_pts += points[i].size();

    stream << "ply" << std::endl;
    stream << "format ascii 1.0" << std::endl;
    stream << "element vertex " << nb_pts << std::endl;
    stream << "property float x" << std::endl;
    stream << "property float y" << std::endl;
    stream << "property float z" << std::endl;
    stream << "property float nx" << std::endl;
    stream << "property float ny" << std::endl;
    stream << "property float nz" << std::endl;
    stream << "property uchar red" << std::endl;
    stream << "property uchar green" << std::endl;
    stream << "property uchar blue" << std::endl;
    stream << "end_header" << std::endl;

    // Vertices
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < points[i].size(); ++j) {
                stream  << get<0>(points[i][j]).x() << " " << get<0>(points[i][j]).y() << " " << get<0>(points[i][j]).z()
                << " "  << get<1>(points[i][j]).x() << " " << get<1>(points[i][j]).y() << " " << get<1>(points[i][j]).z()
                << " "  << (int)colors[i].r() << " " << (int)colors[i].g() << " " << (int)colors[i].b() << std::endl;
        }
    }

    stream << std::endl;
    stream.close();
}


void Shape_Container::to_ply(const string & filename, const string &  type)
{
    if(type == "pointcloud"){
        inliers_to_ply(filename);
        return;
    }

    vector<vector<Inexact_Point_3> > points;
    vector<CGAL::Color> colors;

    size_t total_convex_hulls = size_of_inexact_convex_hulls();
    points.reserve(total_convex_hulls);
    colors.reserve(total_convex_hulls);

    for (size_t p = 0; p < total_convex_hulls; ++p) {
        if (is_primitive_degenerate(p)) continue;

        const CGAL::Color & C = get_primitive_color(p);


        if (type == "alpha") {
            // Alpha-shapes
            std::pair<int, int> tr_range = get_range_of_triangles(p);
            for (int tr_id = tr_range.first; tr_id <= tr_range.second; ++tr_id) {
                const std::vector<Inexact_Point_3> & T = get_triangle(tr_id);
                if (T.size() == 3) {
                    points.push_back(T);
                    colors.push_back(C);
                }
            }
        }

        else if (type == "rectangle") {
            // Best rectangles
            const std::vector<Inexact_Point_3> & H = get_best_rectangle_points(p);
            points.push_back(H);
            colors.push_back(C);
        }

        else if (type == "convex") {
            // Convex hulls
            const std::vector<Inexact_Point_3> & H = get_convex_hull_points(p);
            points.push_back(H);
            colors.push_back(C);
        }
    }

    auto pp = fs::path(filename).parent_path();
    if(!fs::is_directory(pp))
        fs::create_directories(pp);
    std::ofstream stream(filename, std::ios::out);
    if (!stream.is_open()){
        throw std::ios_base::failure("Error : cannot write into an output file");
    }

    double dims[3];
    dims[0] = SD->get_x_max() - SD->get_x_min();
    dims[1] = SD->get_y_max() - SD->get_y_min();
    dims[2] = SD->get_z_max() - SD->get_z_min();

    int lg_max = 0;
    for (int i = 0; i < 3; i++) {
        int lg = (dims[i] == 0 ? 0 : log10(dims[i]));
        if (lg > lg_max) {
            lg_max = lg;
        }
    }

    stream.precision(10 + lg_max);

    size_t nb_pts = 0;
    for (size_t i = 0; i < points.size(); i++) nb_pts += points[i].size();

    stream << "ply" << std::endl;
    stream << "format ascii 1.0" << std::endl;
    stream << "element vertex " << nb_pts << std::endl;
    stream << "property float x" << std::endl;
    stream << "property float y" << std::endl;
    stream << "property float z" << std::endl;
    stream << "property uchar red" << std::endl;
    stream << "property uchar green" << std::endl;
    stream << "property uchar blue" << std::endl;
    stream << "element face " << points.size() << std::endl;
    stream << "property list uchar int vertex_index" << std::endl;
    stream << "end_header" << std::endl;

    // Vertices
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < points[i].size(); ++j) {
            stream << points[i][j].x() << " " << points[i][j].y() << " " << points[i][j].z()
                << " " << (int)colors[i].r() << " " << (int)colors[i].g() << " " << (int)colors[i].b() << std::endl;
        }
    }

    // Facets
    size_t cont = 0;
    for (size_t i = 0; i < points.size(); i++) {
        stream << points[i].size() << " ";
        for (size_t k = cont; k < cont + points[i].size(); k++) {
            stream << k << " ";
        }
        stream << std::endl;
        cont = cont + points[i].size();
    }

    stream << std::endl;
    stream.close();
}

int Shape_Container::save(const string& filename, const string& type){

    string extension = boost::filesystem::extension(filename);

    SD->_logger->info("Save planar shapes to {}", filename);

    try{
        if (extension == ".ply")
            to_ply(filename,type);
        else if (extension == ".vg" || extension == ".bvg")
            SD->to_vg(filename);
        else if (extension == ".npz")
            SD->to_npz(filename);
        else{
            SD->_logger->error("{} is not a valid extension for the plane file. Allowed are '.ply', '.vg', '.bvg' and '.npz'.",extension);
            return 1;
        }
    }
    catch (std::exception & e){
        SD->_logger->error(e.what());
        SD->_logger->error("Could not save file {}",filename);
        return 1;
    }

    return 0;

}







