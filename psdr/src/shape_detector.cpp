#include "shape_detector.h"
#include "shape_detector_index_map.h"
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Shape_regularization/regularize_planes.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Classification/property_maps.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <random>
#include <math.h>


using namespace std;

Shape_Detector::Shape_Detector()
{

    if(spdlog::get("Shape Detector")){
        _logger = spdlog::get("Shape Detector");
    }
    else{
        _logger = spdlog::stdout_color_mt("Shape Detector");
    }

    #ifdef NDEBUG
    // nondebug
    #else
    spdlog::set_level(spdlog::level::debug);
    #endif
    spdlog::set_pattern("[%H:%M:%S] [%n] [%l] %v");

    // regularization parameters
    lambda = 0.8;
    tolerance_angle = 0.5;


    path_point_cloud = "";
    path_point_cloud_extension = "";
    path_clusters = "";

    knn = -1;
    should_compute_knn = true;
    if_stop = false;

    if_constraint = false;

    lambda_r = 1.0;
    weight_mode = 0;
    lambda_c = 1.0;
    lambda_fidelity = 1.0;
    lambda_regular = 1.0;

    metric_type = 0;

    t_exlude = 0;
    t_insert = 0;
    t_merge = 0;
    t_insert = 0;
    t_split = 0;
    all_t_transfer = 0;
    all_t_exlude = 0;
    all_t_insert = 0;
    all_t_split = 0;
    all_t_merge = 0;
    t_regularization = 0;
    all_regularization = 0;

    mean_error = 0;
    mean_normal_diviation = 0;

    mean_distance_diaviation = 0;

    number_of_insert_exclude = 10;
    number_iterations = 0;
    interval_all = 0;
    if_oriented_normal = true;

}

int Shape_Detector::load_points(const std::string _filename)
{
    path_point_cloud = _filename;
    path_point_cloud_extension = boost::filesystem::extension(path_point_cloud);

    _logger->info("Load points from {}", _filename);

	points.clear();

    try{
        if (path_point_cloud_extension == ".ply")
            load_ply();
        else if (path_point_cloud_extension == ".vg" || path_point_cloud_extension == ".bvg")
            load_vg();
        else if (path_point_cloud_extension == ".npz")
            load_npz();
        else{
            _logger->error("{} is not a valid points file",path_point_cloud_extension);
            return 1;
        }
    }
    catch (std::exception & e){
        _logger->error("Could not load file {}",path_point_cloud);
        _logger->error(e.what());
        return 1;
    }
	
	inliers_to_natural_colors = std::vector<CGAL::Color>(points.size(), CGAL::black());
	spacing_is_known = false;
	set_extrema();

    return 0;
}


bool Shape_Detector::load_ply()
{
	std::ifstream streamb(path_point_cloud);
	int line = 0;
	std::string s;
	while (streamb&&line < 2) {
		if (!getline(streamb, s)) break;
		line++;
	}


	if (s == "format ascii 1.0") {

		std::ifstream stream(path_point_cloud);
        if (!stream || !CGAL::IO::read_PLY(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
            return 1;
		}
	}
	else {
		std::ifstream stream(path_point_cloud, std::ios_base::binary);


        if (!stream || !CGAL::IO::read_PLY(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
            return 1;
		}
	}
	if (points[0].second == Inexact_Vector_3(0, 0, 0)) {
        _logger->warn("input .ply file does not contain valid 'normals'. 'normals' will be estimated.");

        _logger->debug("estimate global k neighbor scale");
        // first estimate k for knn
        vector<Inexact_Point_3> tpoints;
        for(const auto pnv : points)
            tpoints.push_back(pnv.first);
        knn = CGAL::estimate_global_k_neighbor_scale(tpoints);
        should_compute_knn = false;

        _logger->debug("estimate normals using jets");

        CGAL::jet_estimate_normals<Concurrency_tag>
			(points,
				knn, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
				CGAL::parameters::point_map(Point_map())
				.normal_map(Normal_map()));
		if_oriented_normal = false;
	}
	else {
		if_oriented_normal = true;
	}

    return 0;
}


bool Shape_Detector::load_vg()
{
	bool is_binary = (path_point_cloud_extension == ".bvg");

	// Step 1.
	// Reads the file for the first time, gets positions of elements to read

	std::ifstream stream;
	if (is_binary) {
		stream.open(path_point_cloud, std::ios::binary);
	} else {
		stream.open(path_point_cloud);
	}

    if (!stream.is_open()) return 1;

	int N, G;
	int line_num_points, line_num_colors, line_num_normals, line_num_groups;
	int id_line = 1;

	while (!stream.eof()) {
		std::string line;
		std::getline(stream, line);

		std::vector<std::string> split_line;
		boost::split(split_line, line, boost::is_any_of(": "));

		if (split_line.size() > 1) {
			if (split_line[0] == "num_points") {
				line_num_points = id_line;
				N = stoi(split_line[split_line.size() - 1]);
			} else if (split_line[0] == "num_colors") {
				line_num_colors = id_line;
			} else if (split_line[0] == "num_normals") {
				line_num_normals = id_line;
			} else if (split_line[0] == "num_groups") {
				line_num_groups = id_line;
				G = stoi(split_line[split_line.size() - 1]);
			}
		}

		++id_line;
	}

	stream.close();

	// Step 2.
	// Parses elements of the file.

	if (is_binary) {
		stream.open(path_point_cloud, std::ios::binary);
	} else {
		stream.open(path_point_cloud);
	}

	std::vector<Inexact_Point_3> vg_points;
	std::vector<Inexact_Vector_3> vg_normals;
	vg_points.reserve(N);
	vg_normals.reserve(N);

	std::vector<int> vg_inliers_to_planes(N, -1);
	std::vector<std::vector<int> > vg_planes_to_inliers;
	vg_planes_to_inliers.reserve(G);

    vector<Inexact_Plane> vg_planes;
    vg_planes.reserve(G);

    vector<CGAL::Color> vg_plane_colors;
    vg_plane_colors.reserve(G);

	int points_to_read = line_num_colors - line_num_points;
	int colors_to_read = line_num_normals - line_num_colors;
	int normals_to_read = line_num_groups - line_num_normals;

	bool reading_points, reading_colors, reading_normals, reading_groups;
	reading_points = reading_colors = reading_normals = reading_groups = false;

	id_line = 1;

    std::stringstream linestream;
	while (!stream.eof()) {
		
		std::string line;
		std::vector<std::string> split_line;

		if (!reading_points && !reading_colors && !reading_normals && !reading_groups) {
			std::getline(stream, line);
			boost::split(split_line, line, boost::is_any_of(": "));

			if (split_line.size() > 1) {
				if (split_line[0] == "num_points") {
					reading_points = true; reading_colors = false; reading_normals = false; reading_groups = false;
				} else if (split_line[0] == "num_colors") {
					reading_points = false; reading_colors = true; reading_normals = false; reading_groups = false;
				} else if (split_line[0] == "num_normals") {
					reading_points = false; reading_colors = false; reading_normals = true; reading_groups = false;
				} else if (split_line[0] == "num_groups") {
					reading_points = false; reading_colors = false; reading_normals = false; reading_groups = true;
				}
			}

		} else if (reading_points) {
			// Reads points
			if (points_to_read == 2) {
				std::getline(stream, line);
				boost::split(split_line, line, boost::is_any_of(" "));
				for (int k = 0 ; k < N ; ++k) {
					double x = std::stod(split_line[3 * k + 0]);
					double y = std::stod(split_line[3 * k + 1]);
					double z = std::stod(split_line[3 * k + 2]);
					Inexact_Point_3 pt_k(x, y, z);
					vg_points.push_back(pt_k);
				}
			} else {
				for (int k = 0 ; k < N ; ++k) {
					std::getline(stream, line);
					boost::split(split_line, line, boost::is_any_of(" "));
					double x = std::stod(split_line[0]);
					double y = std::stod(split_line[1]);
					double z = std::stod(split_line[2]);
					Inexact_Point_3 pt_k(x, y, z);
					vg_points.push_back(pt_k);
				}
			}
			reading_points = false;

		} else if (reading_colors) {
			// Skips lines
			for (int k = 1 ; k < colors_to_read ; ++k) std::getline(stream, line);
			reading_colors = false;

		} else if (reading_normals) {
			// Reads normals
			if (normals_to_read == 2) {
				std::getline(stream, line);
				boost::split(split_line, line, boost::is_any_of(" "));
				for (int k = 0 ; k < N ; ++k) {
					double x = std::stod(split_line[3 * k + 0]);
					double y = std::stod(split_line[3 * k + 1]);
					double z = std::stod(split_line[3 * k + 2]);
					Inexact_Vector_3 n_k(x, y, z);
					vg_normals.push_back(n_k);
				}
			} else {
				for (int k = 0 ; k < N ; ++k) {
					std::getline(stream, line);
					boost::split(split_line, line, boost::is_any_of(" "));
					double x = std::stod(split_line[0]);
					double y = std::stod(split_line[1]);
					double z = std::stod(split_line[2]);
					Inexact_Vector_3 n_k(x, y, z);
					vg_normals.push_back(n_k);
				}
			}
			reading_normals = false;

		} else if (reading_groups) {
			for (int k = 0 ; k < G ; ++k) {
				// Lines 1 - 5 : skip
//                for (int l = 0; l < 5; ++l) std::getline(stream, line);

                // Line 1 and 2 : skip
                std::getline(stream, line);
                std::getline(stream, line);
                std::getline(stream, line);
                // Line 3 : read plane params
                boost::split(split_line, line, boost::is_any_of(": "));
                Inexact_Plane P(stod(split_line[2]), stod(split_line[3]), stod(split_line[4]), stod(split_line[5]));
                vg_planes.push_back(P);


                std::getline(stream, line);
                std::getline(stream, line);
                // Line 5 : read plane colors
                boost::split(split_line, line, boost::is_any_of(": "));
                CGAL::Color col(stoi(split_line[2]),stoi(split_line[3]),stoi(split_line[4]));
                vg_plane_colors.push_back(col);



				// Line 6 : read number of points
				std::getline(stream, line);
				boost::split(split_line, line, boost::is_any_of(": "));
				int size_group = stoi(split_line[split_line.size() - 1]);

				// Line 7 : read indices of inliers
				std::vector<int> inliers;
				inliers.reserve(size_group);
				std::getline(stream, line);
				boost::split(split_line, line, boost::is_any_of(" "));
				for (int l = 0 ; l < size_group ; ++l) {
					int id_vertex = stoi(split_line[l]);
					inliers.push_back(id_vertex);
					vg_inliers_to_planes[id_vertex] = k;
				}

				// Line 8 : skip
				std::getline(stream, line);

				// Save inliers
				vg_planes_to_inliers.push_back(inliers);
			}
			reading_groups = false;
		}
		++id_line;
	}

	stream.close();

	// Step 3.
	// Merges

	points.clear();
	points = std::vector<Point_with_normal>(N);
	for (int k = 0 ; k < N ; ++k) points[k] = std::make_pair(vg_points[k], vg_normals[k]);
	vg_points.clear();
	vg_normals.clear();

    inliers_to_planes = vg_inliers_to_planes;
	planes_to_inliers.clear();
	planes_to_inliers = vg_planes_to_inliers;

    planes_2.clear();
    planes_2 = vg_planes;
    planes_to_colors.clear();
    planes_to_colors = vg_plane_colors;

    return 0;
}

#include <xtensor-io/xnpz.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
bool Shape_Detector::load_npz()
{


    /////// with xtensor
    auto a = xt::load_npz(path_point_cloud);

    if(a.find("points") == a.end()){
        _logger->error("No points array found in .npz file {}.", path_point_cloud);
        return 1;
    }
    if(a.find("normals") == a.end()){
        _logger->error("No normals array found in .npz file {}.", path_point_cloud);
        return 1;
    }

    auto pts = a["points"].cast<double>();
    auto normals = a["normals"].cast<double>();

    points.clear();
    points = std::vector<Point_with_normal>(pts.shape()[0]);

    Inexact_Point_3 p;
    Inexact_Vector_3 n;
    for(int i = 0; i < pts.shape()[0]; i++){
        p = Inexact_Point_3(pts(i,0),pts(i,1),pts(i,2));
        n = Inexact_Vector_3(normals(i,0),normals(i,1),normals(i,2));
        points[i] = std::make_pair(p, n);
    }


    vector<int> vg_inliers_to_planes(pts.shape()[0], -1);
    vector<vector<int>> vg_planes_to_inliers;
    vector<Color> vg_planes_to_colors;
    vector<vector<double>> vg_planes;
    auto group_num_points = a["group_num_points"].cast<int32_t>();
    auto group_points = a["group_points"].cast<int32_t>();
    auto group_colors = a["group_colors"].cast<int32_t>();
    auto group_parameters = a["group_parameters"].cast<double>();
    int id = 0;
    int current_point_id = 0;
    for(int i = 0; i < group_num_points.shape()[0]; i++){

        std::vector<int> this_plane;
        for(int j = 0; j < group_num_points[i]; j++){
            current_point_id = group_points[j+id];
            vg_inliers_to_planes[current_point_id] = i;
            this_plane.push_back(current_point_id);
        }
        vg_planes_to_inliers.push_back(this_plane);
        id+=group_num_points[i];
        auto color = Color(group_colors(i,0),group_colors(i,1),group_colors(i,2));
        vg_planes_to_colors.push_back(color);
        vg_planes.push_back({group_parameters(i,0),group_parameters(i,1),group_parameters(i,2),group_parameters(i,3)});
    }

    inliers_to_planes.clear();
    inliers_to_planes = vg_inliers_to_planes;
    planes_to_inliers.clear();
    planes_to_inliers = vg_planes_to_inliers;
    planes_to_colors = vg_planes_to_colors;
    planes = vg_planes;

    return 0;

}


void Shape_Detector::compute_average_spacing_and_k_nearest_neighbors()
{
	spherical_neighborhood.clear();
	spherical_neighborhood.reserve(points.size());

	std::map<Inexact_Point_3, int> map_indice_point;
	std::list<Inexact_Point_3> list_points;

	for (int i = 0; i < points.size(); i++) {
		const Inexact_Point_3 & pt = points[i].first;
		map_indice_point[pt] = i;
		list_points.push_back(pt);
	}

	Tree tree(list_points.begin(), list_points.end());
	if (!spacing_is_known) average_spacing = 0;

    _logger->debug("Computing K nearest neighbors...");

	for (int i = 0; i < points.size(); i++) {

		Inexact_Point_3 query = points[i].first;
		Neighbor_search search(tree, query, knn + 1);

		double dist_1st_neighbor = FLT_MAX;
		std::vector<int> index_of_neighbors;
		index_of_neighbors.reserve(knn);
		
		int terms = 0;

		for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
			std::map<Inexact_Point_3, int>::iterator iter = map_indice_point.begin();
			iter = map_indice_point.find(it->first);
			if (iter != map_indice_point.end() && iter->second != i) {
				index_of_neighbors.push_back(iter->second);
				
				if (!spacing_is_known) {
					double d = sqrt((query - iter->first).squared_length());
					if (d < dist_1st_neighbor) dist_1st_neighbor = d;
				}

				++terms;
			}
		}

		spherical_neighborhood.push_back(index_of_neighbors);
		if (!spacing_is_known && terms > 0) average_spacing += dist_1st_neighbor;
	}


	if (!spacing_is_known) {
		average_spacing /= int(points.size());
        _logger->debug("Average spacing : {}",average_spacing);
		spacing_is_known = true;
	}

	should_compute_knn = false;
}



void Shape_Detector::reshuffle_points_based_on_planarity()
{
	// Step 1.
	// Compute planarity values

	typedef std::pair<size_t, double> Planarity_Value;
	std::vector<Planarity_Value> P;

	P.reserve(points.size());
	for (size_t i = 0 ; i < points.size() ; ++i) {
		const Inexact_Point_3 & p_i = points[i].first;
		const Inexact_Vector_3 & n_i = points[i].second;
		const std::vector<int> & s_i = spherical_neighborhood[i];
		
		Inexact_Plane H_i(p_i, n_i);

		double d_i = 0;
		for (size_t j = 0 ; j < s_i.size() ; ++j) {
			const Inexact_Point_3 & p_j = points[s_i[j]].first;
			const Inexact_Vector_3 & n_j = points[s_i[j]].second;
			d_i = fabs((H_i.a() * p_j.x() + H_i.b() * p_j.y() + H_i.c() * p_j.z() + H_i.d()) / (n_i * n_j));
			//d_i += sqrt(CGAL::squared_distance(H_i, p_j));
		}
		d_i /= s_i.size();

		P.push_back(std::make_pair(i, d_i));
	}

	// Step 2.
	// Sorts the planarity table

	struct _Planarity_Comparator {
		bool operator() (const Planarity_Value & L, const Planarity_Value & R) {
			return L.second < R.second;
		}
	} Planarity_Comparator;

	std::sort(P.begin(), P.end(), Planarity_Comparator);

	// Step 3.
	// Reorders points

	std::vector<size_t> new_to_old_indices(points.size(), -1);
	std::vector<size_t> old_to_new_indices(points.size(), -1);

	std::vector<Point_with_normal> ordered_points;
	ordered_points.reserve(points.size());
	
	for (size_t i = 0 ; i < P.size() ; ++i) {
		size_t j = P[i].first;
		ordered_points.push_back(points[j]);
		old_to_new_indices[j] = i;
		new_to_old_indices[i] = j;
	}

	points = ordered_points;
	ordered_points.clear();

	std::vector<std::vector<int> > ordered_spherical_neighborhood;
	ordered_spherical_neighborhood.reserve(points.size());

	for (size_t i = 0 ; i < points.size() ; ++i) {
		const std::vector<int> & N_0 = spherical_neighborhood[new_to_old_indices[i]];

		std::vector<int> N_1;
		N_1.reserve(N_0.size());

		for (size_t j = 0 ; j < N_0.size() ; ++j) {
			N_1.push_back(old_to_new_indices[N_0[j]]);
		}

		ordered_spherical_neighborhood.push_back(N_1);
	}
	
	spherical_neighborhood = ordered_spherical_neighborhood;
}



void Shape_Detector::set_extrema()
{
	for (size_t i = 0 ; i < points.size() ; ++i) {
		const Inexact_Vector_3 & n = points[i].second;
		const double &x = n.x(), &y = n.y(), &z = n.z();
		double l = sqrt(x * x + y * y + z * z);
		points[i].second = Inexact_Vector_3 (x / l, y / l, z / l);
	}

	x_min = y_min = z_min = FLT_MAX;
	x_max = y_max = z_max = -FLT_MAX;
	for (size_t i = 0 ; i < points.size() ; ++i) {
		const Inexact_Point_3 & pt = points[i].first;
		const double &x = pt.x(), &y = pt.y(), &z = pt.z();
		if (x < x_min) x_min = x;
		if (x > x_max) x_max = x;
		if (y < y_min) y_min = y;
		if (y > y_max) y_max = y;
		if (z < z_min) z_min = z;
		if (z > z_max) z_max = z;
	}
	
    _logger->debug("X : [{},{}] ", x_min ,x_max);
    _logger->debug("Y : [{},{}] ", y_min ,y_max);
    _logger->debug("Z : [{},{}] ", z_min ,z_max);

	double dx = x_max - x_min, dy = y_max - y_min, dz = z_max - z_min;
	bbox_diagonal = sqrt(dx * dx + dy * dy + dz * dz);

    _logger->debug("Bounding box diagonal : {}", bbox_diagonal);
}



int Shape_Detector::set_detection_parameters(int _rg_min_points, double _rg_epsilon, int _knn, double _rg_normal_threshold)
{

    knn = _knn;
    should_compute_knn = true; // somehow this needs to always stay on because it also computes a so called 'spherical neighborhood'

	min_points = _rg_min_points;
	epsilon = _rg_epsilon;
	normal_threshold = _rg_normal_threshold;

    // for refinement
    tolerance_coplanarity = epsilon / 0.2;

    return 0;
}


void Shape_Detector::set_knn( int _knn)
{
	if (knn != _knn) {
		knn = _knn;
		
	}

	
}



void Shape_Detector::set_regularization_parameters(
	double _lambda,
	double _tolerance_angle,
	double _tolerance_coplanarity)
{
	lambda = _lambda;
	tolerance_angle = _tolerance_angle;
	tolerance_coplanarity = _tolerance_coplanarity;

    _logger->debug( "lambda = {}", lambda);
    _logger->debug( "tolerance_angle  = {}", tolerance_angle);
    _logger->debug( "tolerance_coplanarity  = {}", tolerance_coplanarity);
}



void Shape_Detector::set_discretization_parameters(double _discretization_angle, double _discretization_distance_ratio)
{
	discretization_angle = _discretization_angle;
	discretization_distance = _discretization_distance_ratio * bbox_diagonal / 100.0;

    _logger->debug( "discretization_angle = {}", discretization_angle);
    _logger->debug( "discretization_dist  = {}", discretization_distance);
}


void Shape_Detector::detect_shapes()
{
	detect_planes(false);
	
	discretize_planes();
	//initialized the degree of planes' freedom.
	initial_regular_groups_done();

	get_coverage_and_mean_error();
	ori_all_error = all_error;
	ori_coverage = coverage;
	ori_mean_error = mean_error;
	ori_inliers_number = number_of_assigned_points;
	ori_primitives_number = primitives_number;
	ori_mean_normal_diviation = mean_normal_diviation;
	ori_freedom_of_planes = freedom_of_planes;
	get_last_information();

}


void Shape_Detector::load_shapes()
{
	detect_planes(true);
	discretize_planes();
	initial_regular_groups_done();

	get_coverage_and_mean_error();
	ori_all_error = all_error;
	ori_coverage = coverage;
	ori_mean_error = mean_error;
	ori_inliers_number = number_of_assigned_points;
	ori_primitives_number = primitives_number;
	ori_mean_normal_diviation = mean_normal_diviation;
	ori_freedom_of_planes = freedom_of_planes;
	get_last_information();

}


void Shape_Detector::regularize_shapes()
{
	regularize_planes();
	discretize_planes();
	get_coverage_and_mean_error();
}


void Shape_Detector::refine_shapes(int mi) {

    set_max_iter(mi);

	initial_regular_groups_done();
	clock_t t_start = clock();
	if (metric_type == 0) {
		planar_shape_detection_hybrid();
    }
	else if (metric_type == 1) {
		planar_shape_detection_L1();
	}
	else {
		planar_shape_detection_l2();
	}
	clock_t t_end = clock();

	double t_all = double(t_end - t_start) / CLOCKS_PER_SEC;
	
	show_result(t_all);

	if_stop = false;

	//update color
	planes_to_colors.clear();
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(100, 225);
	for (size_t i = 0; i < planes_2.size(); ++i) {
		unsigned char r = 0, g = 0, b = 0;
        r = uniform_distribution(generator);
        g = uniform_distribution(generator);
        b = uniform_distribution(generator);
		planes_to_colors.push_back(CGAL::Color(r, g, b));
	}
}

void Shape_Detector::show_result(double t) {

    _logger->debug("All operators: {} s.",t);
    _logger->debug("Transfer operator: {} times",all_t_transfer);
    _logger->debug("Merge operator: {} times",all_t_merge);

    _logger->debug("Split operator: {} times",all_t_split);

    _logger->debug("Insert operator: {} times" ,all_t_insert);
    _logger->debug("Exclude operator: {} times" ,all_t_exlude );
    _logger->debug("Regularization operator: {} times",all_regularization );

    _logger->debug("Number of iterations: {}",number_iterations);
    _logger->debug("Primitives : {}",(primitives_number));
    _logger->debug("Coverage  : {}",(coverage));
    _logger->debug("Mean error : {}" ,(mean_error) );
    _logger->debug("Mean normal deviation : {}" ,(mean_normal_diviation) );
    _logger->debug("Degree of freedom : {}" ,(freedom_of_planes) );

    _logger->debug("Primitives reducing : {}" , (last_primitives_number - primitives_number));
    _logger->debug("Coverage adding   : {}", (coverage - last_coverage));
    _logger->debug("Mean error reducing : {}" ,(last_mean_error - mean_error) );
    _logger->debug("Mean normal deviation adding : {}" , (mean_normal_diviation - last_normal_deviation) );
    _logger->debug("Degree of freedom changing: {}" , (freedom_of_planes-last_freedom_planes));
}

void Shape_Detector::calculate_center() {
	planes_centroids.clear();
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

		//Inexact_Point_3 centroid = CGAL::ORIGIN;
		double xc = 0, yc = 0, zc = 0;

		

		for (size_t j = 0; j < planes_to_inliers[i].size(); ++j) {
			const Inexact_Point_3 & pt = points[planes_to_inliers[i][j]].first;
			
			xc += pt.x(), yc += pt.y(), zc += pt.z();
		}

		xc /= planes_to_inliers[i].size();
		yc /= planes_to_inliers[i].size();
		zc /= planes_to_inliers[i].size();
		Inexact_Point_3 centroid(xc, yc, zc);
		planes_centroids.push_back(centroid);
	}

}



int Shape_Detector::change_of_degrees_of_freedom_after_add_remove(int id_plane) {
	

		if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] == -1) {
			//no regular plane
			return 0;

		}
		else if (planes_to_coplanar_done[id_plane] != -1 && planes_to_orthogonal_done[id_plane] != -1 && planes_to_parallel_done[id_plane] == -1) {
			//only parallel
			return 2;

		}
		else if (planes_to_coplanar_done[id_plane] != -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] != -1) {
			//orthogonal to do approximate, should test if orthogonal stil exist
			return 2;

		}
		else if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] != -1) {
			//orthogonal to do approximate, should test if orthogonal stil exist
			return 2+1;

		}
		else if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] != -1 && planes_to_parallel_done[id_plane] == -1) {
			//orthogonal to do approximate, should test if orthogonal stil exist
			return 2 + 1;

		}

        // added this because there was no return before
        return 0;

}

int Shape_Detector::change_of_degrees_of_freedom_after_mergy(int id_1,int id_2) {

	int re_dof = 0;
	if ((planes_to_coplanar_done[id_1] == -1 && planes_to_orthogonal_done[id_1] == -1 && planes_to_parallel_done[id_1] == -1)) {
		//no regular plane
		re_dof += 0;

	}
	else if ((planes_to_coplanar_done[id_1] != -1 && planes_to_orthogonal_done[id_1] != -1 && planes_to_parallel_done[id_1] == -1)) {
		
			//only two parallel
			re_dof += 2;
		

	}
	else if (planes_to_coplanar_done[id_1] != -1 && planes_to_orthogonal_done[id_1] == -1 && planes_to_parallel_done[id_1] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 2;

	}
	else if (planes_to_coplanar_done[id_1] == -1 && planes_to_orthogonal_done[id_1] == -1 && planes_to_parallel_done[id_1] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 3;

	}
	else if (planes_to_coplanar_done[id_1] == -1 && planes_to_orthogonal_done[id_1] != -1 && planes_to_parallel_done[id_1] == -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 3;

	}

	if (planes_to_coplanar_done[id_2] == -1 && planes_to_orthogonal_done[id_2] == -1 && planes_to_parallel_done[id_2] == -1) {
		//no regular plane
		re_dof += 0;

	}
	else if (planes_to_coplanar_done[id_2] != -1 && planes_to_orthogonal_done[id_2] != -1 && planes_to_parallel_done[id_2] == -1) {
		//only parallel
		re_dof += 2;

	}
	else if (planes_to_coplanar_done[id_2] != -1 && planes_to_orthogonal_done[id_2] == -1 && planes_to_parallel_done[id_2] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 2;

	}
	else if (planes_to_coplanar_done[id_2] == -1 && planes_to_orthogonal_done[id_2] == -1 && planes_to_parallel_done[id_2] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 3;

	}
	else if (planes_to_coplanar_done[id_2] == -1 && planes_to_orthogonal_done[id_2] != -1 && planes_to_parallel_done[id_2] == -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		re_dof += 3;

	}
	re_dof -= 3;

	return re_dof;


}


int Shape_Detector::change_of_degrees_of_freedom_after_split(int id_plane) {


	if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] == -1) {
		//no regular plane
		return 3;

	}
	else if (planes_to_coplanar_done[id_plane] != -1 && planes_to_orthogonal_done[id_plane] != -1 && planes_to_parallel_done[id_plane] == -1) {
		//only parallel
		return 5;

	}
	else if (planes_to_coplanar_done[id_plane] != -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		return 5;

	}
	else if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] == -1 && planes_to_parallel_done[id_plane] != -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		return 6;

	}
	else if (planes_to_coplanar_done[id_plane] == -1 && planes_to_orthogonal_done[id_plane] != -1 && planes_to_parallel_done[id_plane] == -1) {
		//orthogonal to do approximate, should test if orthogonal stil exist
		return 6;

	}

    // added this because there was no return before
    return 0;


}


double Shape_Detector::add_changed_error(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	for (int u : befor_ids) {
		dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
	}
	befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());
	std::vector<Inexact_Point_3> after_points;
	for (int n : befor_ids) {
		after_points.push_back(points[n].first);
	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (Inexact_Point_3 u : after_points) {
		dif_after += sqrt(CGAL::squared_distance(after_p, u));
	}

	return (dif_after - dif_befor);


}



double Shape_Detector::remove_changed_error(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	std::vector<int> after_ids;
	for (int u : befor_ids) {
		dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
		/*bool sma = false;
		for (int t : point_ids) {
			if (t == u) {
				sma = true;
				break;
			}
		}
		if (!sma)after_ids.push_back(u);*/
		if (std::find(point_ids.begin(), point_ids.end(), u) == point_ids.end()) {
			after_ids.push_back(u);

		}
	}

	//befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());

	std::vector<Inexact_Point_3> after_points;

	for (int n : after_ids) {
		after_points.push_back(points[n].first);
	}
	Inexact_Plane after_p;


    if(after_points.size() < 4)
        return 100000000.0;
    // bug here, after_points is empty, so I added the if above
	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (Inexact_Point_3 u : after_points) {
		dif_after += sqrt(CGAL::squared_distance(after_p, u));
	}

	return (dif_after - dif_befor);

}
bool Shape_Detector::merge_distance_changed_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids) {
	Inexact_Plane p_i = planes_2[i];
	Inexact_Plane p_j = planes_2[j];
	double dis_before = 0;

	std::vector<Inexact_Point_3> assigned_pts;
	std::vector<int> assigned_pts_id;

	assigned_pts.reserve(planes_to_inliers[i].size() + planes_to_inliers[j].size());
	for (int pt_index : planes_to_inliers[i]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);

		dis_before += sqrt(CGAL::squared_distance(p_i, points[pt_index].first));
	}

	for (int pt_index : planes_to_inliers[j]) {
		assigned_pts_id.push_back(pt_index);

		assigned_pts.push_back(points[pt_index].first);
		dis_before += sqrt(CGAL::squared_distance(p_j, points[pt_index].first));

	}





	size_t n_poly = assigned_pts.size();



	double dP_norm = 0;
	/*Inexact_Plane test_plane;

	linear_least_squares_fitting_3(assigned_pts.begin(), assigned_pts.end(), test_plane, CGAL::Dimension_tag<0>());
	for (int p : assigned_pts_id) {
		dP_norm += sqrt(CGAL::squared_distance(test_plane, points[p].first));

	}*/
	Inexact_Plane test_plane;

	linear_least_squares_fitting_3(assigned_pts.begin(), assigned_pts.end(), test_plane, CGAL::Dimension_tag<0>());



	Inexact_Plane new_plane;
	std::vector<Inexact_Point_3> new_assigned_pts;
	std::vector<int> new_assigned_pts_ids;
	for (int pt : assigned_pts_id) {
		if (sqrt(CGAL::squared_distance(test_plane, points[pt].first)) > epsilon) {
			move_ids.push_back(pt);
		}
		else {
			new_assigned_pts.push_back(points[pt].first);
			new_assigned_pts_ids.push_back(pt);

		}
	}

	if (new_assigned_pts_ids.size() < min_points) {
		dif = FLT_MAX;
		return false;
	}
	linear_least_squares_fitting_3(new_assigned_pts.begin(), new_assigned_pts.end(), new_plane, CGAL::Dimension_tag<0>());
	for (int p : new_assigned_pts_ids) {
		dP_norm += sqrt(CGAL::squared_distance(new_plane, points[p].first));

	}


	dif = dP_norm - dis_before;

	return true;




}



bool Shape_Detector::convert_form_merge_l2(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(i);
	o.push_back(j);
	double dif_m = 0;
	std::vector<int> merge_moved_ids;
	//we suppose that the pair of primitives can not merged,when their normal cosine smaller than normal_threshold.
	if (abs(planes_2[i].orthogonal_vector()*planes_2[j].orthogonal_vector()) < normal_threshold) {

        return false;
	}
	std::vector<double> s_d;
	bool if_good = merge_distance_changed_with_epsilon(i, j, dif_m, merge_moved_ids);
	if (!if_good) return false;
	/*if (merge_moved_ids.size() > points.size()*0.01|| merge_moved_ids.size()>planes_to_inliers[i].size()/4|| merge_moved_ids.size() > planes_to_inliers[j].size() / 4) {
		return false;
	}*/



	//	return false;
	//}
	double change_ccc = energy_changed_merge(dif_m, -1.0, merge_moved_ids.size(),i,j);
	if (change_ccc >= 0) {
		return false;
	}

	s_d.push_back(change_ccc);
	o.insert(o.end(), merge_moved_ids.begin(), merge_moved_ids.end());

	s_d.push_back(dif_m);
	one_element = std::make_pair(std::make_pair(1, o), s_d);

	return true;
}
bool Shape_Detector::separate_two_out(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif) {
	std::vector<int> this_region = planes_to_inliers[id];
	Inexact_Plane this_plane = planes_2[id];

	double max_value = FLT_MIN;
	int max_index = -1;
	double min_value = FLT_MIN;
	int min_index = -1;
	double dis_bef = 0;
	//divide into two regions according to the two farthest points.
	double this_mean_normal_diviation = 0;
	Inexact_Vector_3 plane_direction = this_plane.orthogonal_vector();

	for (int ii = 0; ii < this_region.size(); ++ii) {
		double this_dis = sqrt(CGAL::squared_distance(this_plane, points[this_region[ii]].first));
		double this_normal_d = std::abs(points[this_region[ii]].second*plane_direction);
		this_mean_normal_diviation += this_normal_d;
		dis_bef += this_dis;
		if (this_plane.has_on_positive_side(points[this_region[ii]].first)) {
			if (this_dis > max_value) {
				max_value = this_dis;
				max_index = this_region[ii];

			}
		}
		else {
            if (this_dis > min_value) {
//            if (this_dis >= min_value) {
                min_value = this_dis;
				min_index = this_region[ii];
			}

		}
	}
	this_mean_normal_diviation /= double(this_region.size());
	//we do not split the primitive that has smaller normal diviation than ori mean.




	if (this_mean_normal_diviation > ori_mean_normal_diviation) {
		min_list.clear();
		max_list.clear();
		dif = 0;
		return false;
	}

    if (min_index < 0 || max_index < 0) {
        min_list.clear();
        max_list.clear();
        dif = 0;
        return false;
    }

	std::map<int, int> label_points;
	Inexact_Point_3 max_point = points[max_index].first;
    // bug here, min_index = -1, so I added the if above this block
	Inexact_Point_3 min_point = points[min_index].first;
	std::vector<Inexact_Point_3> max_point_list;
	std::vector<Inexact_Point_3> min_point_list;
	Inexact_Point_2 max_point_2d = this_plane.to_2d(max_point);
	Inexact_Point_2 min_point_2d = this_plane.to_2d(min_point);


	for (int j = 0; j < this_region.size(); ++j) {
		if ((max_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length() < (min_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length()) {
			max_list.push_back(this_region[j]);
			label_points[this_region[j]] = 1;
			max_point_list.push_back(points[this_region[j]].first);
		}
		else {
			label_points[this_region[j]] = -1;
			min_point_list.push_back(points[this_region[j]].first);
			min_list.push_back(this_region[j]);
		}
	}

	//transfer between max and min lists.
	Inexact_Plane plane_max;
	Inexact_Plane plane_min;

	if (max_point_list.size() < 3 || min_point_list.size() < 3) {
		return false;
	}


	linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
	linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
	bool propo = false;
	int Nb_neigh = 10;
	int refine_time = 0;

	do {
		refine_time++;
		int moving_n = 0;
		propo = true;
		std::vector<int> if_moved = std::vector<int>(this_region.size(), -1);
		for (int f = 0; f < this_region.size(); ++f) {

			if (if_moved[f] > 2) continue;
			int this_label = label_points[this_region[f]];
			if ((int)spherical_neighborhood[this_region[f]].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[this_region[f]].size();
			}
			bool bb = false;
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[this_region[f]][it];
				if (inliers_to_planes[neighbor_index] != id) continue; //the neighbor point should be in the splitted region
				if (label_points[neighbor_index] != this_label) {
					bb = true;
					break;
				}


			}
			if (bb == false) continue;
			Inexact_Point_3 this_point = points[this_region[f]].first;



			if (this_label == -1) {

				if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) > normal_threshold) {
					if (CGAL::squared_distance(plane_max, this_point) < CGAL::squared_distance(plane_min, this_point)) {																																														 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f] ++;
					}
				}
			}
			else {

				if (abs(points[this_region[f]].second * plane_min.orthogonal_vector()) > normal_threshold) {
					if (CGAL::squared_distance(plane_min, this_point) < CGAL::squared_distance(plane_max, this_point)) {																																													 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f] ++;
					}
				}
			}

		}

		//update the two regions
		std::map<int, int>::iterator iter;
		iter = label_points.begin();
		max_list.clear();
		min_list.clear();
		max_point_list.clear();
		min_point_list.clear();
		while (iter != label_points.end()) {
			if (iter->second == 1) {
				max_list.push_back(iter->first);
				max_point_list.push_back(points[iter->first].first);

			}
			else {
				min_list.push_back(iter->first);
				min_point_list.push_back(points[iter->first].first);
			}
			iter++;
		}
		if (max_point_list.size() < min_points || min_point_list.size() < min_points) {

			return false;
		}
		linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
		linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
		if (moving_n < 5 || refine_time>25) {
			propo = false;

		}

	} while (propo);

	double dis_after = 0;

	for (int i = 0; i < max_list.size(); ++i) {

		dis_after += sqrt(CGAL::squared_distance(plane_max, points[max_list[i]].first));

	}

	for (int i = 0; i < min_list.size(); ++i) {

		dis_after += sqrt(CGAL::squared_distance(plane_min, points[min_list[i]].first));

	}

	dif = dis_after - dis_bef;

	return true;


}


bool Shape_Detector::convert_form_split_l2(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element) {
	std::vector<int> max_list;
	std::vector<int> min_list;
	double dif_s = 0;
	bool if_splite = separate_two_out(i, max_list, min_list, dif_s);
	double en_s = energy_changed(dif_s, 1.0, i);
	if (en_s >= 0 || max_list.size() < min_points || min_list.size() < min_points || !if_splite) {
		return false;
	}
	max_list_vector.push_back(max_list);
	min_list_vector.push_back(min_list);

	std::vector<int> o;
	o.push_back(i);
	o.push_back(max_list_vector.size() - 1);
	std::vector<double> s_d;
	s_d.push_back(en_s);
	s_d.push_back(dif_s);
	one_element = std::make_pair(std::make_pair(2, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_exclude_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(bad_points_shape[i].first.first);

	o.insert(o.end(), bad_points_shape[i].first.second.begin(), bad_points_shape[i].first.second.end());

	one_element = std::make_pair(std::make_pair(3, o), bad_points_shape[i].second);

	return true;
}

bool Shape_Detector::convert_form_insert_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;
	o.push_back(good_points_shape[i].first.first);
	o.insert(o.end(), good_points_shape[i].first.second.begin(), good_points_shape[i].first.second.end());


	one_element = std::make_pair(std::make_pair(4, o), good_points_shape[i].second);
	return true;
}

bool Shape_Detector::convert_form_regularization_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;//id of the parallel cluster + all the planes will be affected if this operation implied.
	o.push_back(i);
	
	if (parallel_id_changed_normal_after_orthogonal[i].size() > 0) {
		for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[i]) {
			o.insert(o.end(), parallel_clusters_to_planes[orth.first].begin(), parallel_clusters_to_planes[orth.first].end());
		}
	}
	else {
		o.insert(o.end(), parallel_clusters_to_planes[i].begin(), parallel_clusters_to_planes[i].end());
	}

	


	one_element = std::make_pair(std::make_pair(5, o), energy_changing_initial_regularization[i]);
	return true;
}
bool Shape_Detector::convert_form_regularization_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;//all the planes will be affected if this operation implied.
	o.push_back(i);
	o.insert(o.end(), parallel_clusters_to_planes[i].begin(), parallel_clusters_to_planes[i].end());
	/*for (std::vector<int> this_coplanar_ids : parallel_cluster_to_coplanr_cases[i]) {
		o.insert(o.end(), this_coplanar_ids.begin(), this_coplanar_ids.end());
	}*/


	for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[i]) {
		o.insert(o.end(), parallel_clusters_to_planes[orth.first].begin(), parallel_clusters_to_planes[orth.first].end());


	}




	one_element = std::make_pair(std::make_pair(5, o), energy_changing_initial_regularization[i]);
	return true;
}

bool Shape_Detector::update_bad_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {



	std::vector<std::pair<int, double>> ids_bad_inliers;
	double bigest_normal = mean_normal_diviation;
	double bigest_dis = 1.5 * mean_distance_diaviation;


	const Inexact_Plane & H = planes_2[id_shape];


	for (int j : planes_to_inliers[id_shape]) {
		const Point_with_normal & pt = points[j];
		if (points_if_added[j] > 2) continue;

		if (sqrt(CGAL::squared_distance(H, pt.first)) > bigest_dis) {
			ids_bad_inliers.push_back(std::make_pair(j, sqrt(CGAL::squared_distance(H, pt.first))));
		}

	}

	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator> p_q;
	for (std::pair<int, double> ppp : ids_bad_inliers) {
		p_q.push(ppp);

	}
	std::vector<int> removed_points_ids;
	double en = 0;
	double dif_e = 0;

	double n_b = 0;
	while (!p_q.empty() && n_b < number_of_insert_exclude) {
		removed_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_b++;
	}
	dif_e = remove_changed_error(id_shape, removed_points_ids);
	en = energy_changed_second(dif_e, n_b,id_shape);

	if (removed_points_ids.size() != 0 && (planes_to_inliers[id_shape].size() - removed_points_ids.size()) >= min_points && en < 0) {
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);
		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), removed_points_ids.begin(), removed_points_ids.end());


		one_element = std::make_pair(std::make_pair(3, o), dd);
		return true;
	}
	else { return false; }



}
bool Shape_Detector::update_good_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	int Nb_neigh = 10;
	std::vector<std::pair<int, double>> ids_good_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (if_insert_candidate[i] == 0) continue;
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;//avoide infinite loop
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;

		bool if_correspond = false;
		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);
				if (inliers_to_planes[neighbor_index] == id_shape) {
					if_correspond = true;
				}
			}


		}

		if (one_neight_id.empty() || !if_correspond) { continue; }

		double mm_dis = epsilon;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > normal_threshold) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < mm_dis) {

					mm_dis = sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first));
					changed_plane_id = neight_id;

				}
			}
		}


		if (changed_plane_id == id_shape) {
			ids_good_outliers.push_back(std::make_pair(i, mm_dis));//point_id, mean distance.

		}

	}


	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator> p_q;

	for (std::pair<int, double> ppp : ids_good_outliers) {
		p_q.push(ppp);
	}

	std::vector<int> added_points_ids;
	double en = 0;
	double dif_e = 0;
	double n_a = 0;

	while (!p_q.empty() && n_a < number_of_insert_exclude) {
		added_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_a++;
	}
	dif_e = add_changed_error(id_shape, added_points_ids);
	en = energy_changed_second(dif_e, -n_a,id_shape);

	if (added_points_ids.size() != 0 && en < 0) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);

		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), added_points_ids.begin(), added_points_ids.end());


		one_element = std::make_pair(std::make_pair(4, o), dd);


		return true;

	}
	else {
		return false;
	}


}


void Shape_Detector::update_regular_done_group_by_planes() {
	//update group
	int number_parallel_group = parallel_done_to_planes.size();
	int number_orthogonal_group = orthogonal_done_to_planes.size();
	int number_coplanar_group = coplanar_done_to_planes.size();

	parallel_done_to_planes.clear();
	orthogonal_done_to_planes.clear();
	coplanar_done_to_planes.clear();
	parallel_done_to_planes.resize(number_parallel_group);
	orthogonal_done_to_planes.resize(number_orthogonal_group);
	coplanar_done_to_planes.resize(number_coplanar_group);

	for (int i = 0; i < planes_to_parallel_done.size(); ++i) {

		if (planes_to_parallel_done[i] != -1 && planes_to_orthogonal_done[i] != -1) {
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wrong";
			getchar();
		}
		if (planes_to_parallel_done[i] != -1) {
			parallel_done_to_planes[planes_to_parallel_done[i]].push_back(i);
		}
		if (planes_to_orthogonal_done[i] != -1) {
			orthogonal_done_to_planes[planes_to_orthogonal_done[i]].push_back(i);
		}
		if (planes_to_coplanar_done[i] != -1) {

			coplanar_done_to_planes[planes_to_coplanar_done[i]].push_back(i);
		}
	}
	std::vector<std::vector<int>> parallel_done_to_planes_new;
	for (std::vector<int> this_parallel_done : parallel_done_to_planes) {
		if (this_parallel_done.size() > 1) {

			parallel_done_to_planes_new.push_back(this_parallel_done);
		}

	}

	if (parallel_done_to_planes_new.size() != parallel_done_to_planes.size()) {
		parallel_done_to_planes = parallel_done_to_planes_new;
		planes_to_parallel_done.clear();
		planes_to_parallel_done = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < parallel_done_to_planes.size(); ++i) {
			for (int id_p : parallel_done_to_planes[i]) {
				planes_to_parallel_done[id_p] = i;
			}

		}

	}


	std::vector<std::vector<int>> orthogonal_done_to_planes_new;
	for (std::vector<int> this_orthogonal_done : orthogonal_done_to_planes) {
		if (this_orthogonal_done.size() > 1) {

			orthogonal_done_to_planes_new.push_back(this_orthogonal_done);
		}

	}

	if (orthogonal_done_to_planes_new.size() != orthogonal_done_to_planes.size()) {
		orthogonal_done_to_planes = orthogonal_done_to_planes_new;
		planes_to_orthogonal_done.clear();
		planes_to_orthogonal_done = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < orthogonal_done_to_planes.size(); ++i) {
			for (int id_p : orthogonal_done_to_planes[i]) {
				planes_to_orthogonal_done[id_p] = i;
			}

		}

	}

	std::vector<std::vector<int>> coplanar_done_to_planes_new;
	for (std::vector<int> this_coplanar_done : coplanar_done_to_planes) {
		if (this_coplanar_done.size() > 1) {

			coplanar_done_to_planes_new.push_back(this_coplanar_done);
		}

	}

	if (coplanar_done_to_planes_new.size() != coplanar_done_to_planes.size()) {
		coplanar_done_to_planes = coplanar_done_to_planes_new;
		planes_to_coplanar_done.clear();
		planes_to_coplanar_done = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < coplanar_done_to_planes.size(); ++i) {
			for (int id_p : coplanar_done_to_planes[i]) {
				planes_to_coplanar_done[id_p] = i;
			}

		}

	}


}

//update detected regular relationships
void Shape_Detector::update_regular_relations_after_merge(int id_1,int id_2, std::vector<int> respective_planes) {



	
	// parallel
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		std::vector<int> this_parallel_cluser;
		
		for (int para_id : parallel_clusters_to_planes[i]) {
			if (para_id != id_1 && para_id != id_2) {
				this_parallel_cluser.push_back(respective_planes[para_id]);
			}
			else {
				if_regularization_can_conducted[i] = false;
			}


		}
		//update parallel_clusters_to_planes
		parallel_clusters_to_planes[i] = this_parallel_cluser;
	

	}

	//orthogonal
	for (int i = 0; i < parallel_id_changed_normal_after_orthogonal.size(); ++i) {
		bool if_this_influence = false;
	
		for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[i]) {
			if (!if_regularization_can_conducted[orth.first])
			{
				if_regularization_can_conducted[i] = false;

			}


		}
	
	}

	//coplanar
	for (int i = 0; i < parallel_cluster_to_coplanr_cases.size(); ++i) {
		std::vector<std::vector<int>> this_coplanar_clusters;
		for (std::vector<int> this_coplanar_ids : parallel_cluster_to_coplanr_cases[i]) {
			std::vector<int> this_coplanar_cluser;
		
			for (int this_coplanar_id : this_coplanar_ids) {
				if (this_coplanar_id != id_1 && this_coplanar_id != id_2) {
					this_coplanar_cluser.push_back(respective_planes[this_coplanar_id]);
				}

			}
			if (this_coplanar_cluser.size() > 1) {
				this_coplanar_clusters.push_back(this_coplanar_cluser);
			}

		}
		parallel_cluster_to_coplanr_cases[i] = this_coplanar_clusters;
	

	}


}

void Shape_Detector::update_regular_relations_after_splite(int one_splite, std::vector<int> respective_planes) {

	// parallel
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		std::vector<int> this_parallel_cluser;
	
		for (int para_id : parallel_clusters_to_planes[i]) {
			if (para_id != one_splite) {
				this_parallel_cluser.push_back(respective_planes[para_id]);
			}
			else {
				if_regularization_can_conducted[i] = false;
			}


		}
		//update parallel_clusters_to_planes
		parallel_clusters_to_planes[i] = this_parallel_cluser;

		

	}
	//orthogonal
	for (int i = 0; i < parallel_id_changed_normal_after_orthogonal.size(); ++i) {


		for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[i]) {
			if (if_regularization_can_conducted[orth.first] == false) {
				if_regularization_can_conducted[i] = false;
			}



		
	}
	
		
	}
	//coplanar
	for (int i = 0; i < parallel_cluster_to_coplanr_cases.size(); ++i) {
		std::vector<std::vector<int>> this_coplanar_clusters;
		for (std::vector<int> this_coplanar_ids : parallel_cluster_to_coplanr_cases[i]) {
			std::vector<int> this_coplanar_cluser;
			
			for (int this_coplanar_id : this_coplanar_ids) {
				if (this_coplanar_id != one_splite) {
					this_coplanar_cluser.push_back(respective_planes[this_coplanar_id]);
				}
				else {
					if_regularization_can_conducted[i] = false;
				}
			}
			if (this_coplanar_cluser.size() > 1) {
				this_coplanar_clusters.push_back(this_coplanar_cluser);
			}

			
		}
		parallel_cluster_to_coplanr_cases[i] = this_coplanar_clusters;
		

	}


}

void Shape_Detector::update_regular_relations_after_add_remove(int one_change_plane) {
	//only parallel, should test there is no orthogonal
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		

		for (int para_id : parallel_clusters_to_planes[i]) {
			if (para_id == one_change_plane) {
				if_regularization_can_conducted[i] = false;
			}
			


		}
		



	}
	//orthogonal
	for (int i = 0; i < parallel_id_changed_normal_after_orthogonal.size(); ++i) {


		for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[i]) {
			if (if_regularization_can_conducted[orth.first] == false) {
				if_regularization_can_conducted[i] = false;
			}




		}


	}
	//coplanar
	for (int i = 0; i < parallel_cluster_to_coplanr_cases.size(); ++i) {
		
		for (std::vector<int> this_coplanar_ids : parallel_cluster_to_coplanr_cases[i]) {
		
			for (int this_coplanar_id : this_coplanar_ids) {
				if (this_coplanar_id == one_change_plane) {
					if_regularization_can_conducted[i] = false;
				}
			}
			


		}
		

	}

}


int Shape_Detector::get_the_recent_degrees_of_freedom() {
	int dof = 0;
	for (int i = 0; i < planes_to_coplanar_done.size(); ++i) {
		//if this plane is not regularized, its dof is 3.
		if (planes_to_coplanar_done[i] == -1 && planes_to_orthogonal_done[i] == -1 && planes_to_parallel_done[i] == -1) {
			dof += 3;

		}
		//if the plane is only paralleled(no coplanar), its sign distance dof is 1.
		else if (planes_to_coplanar_done[i] == -1 && planes_to_orthogonal_done[i] == -1 && planes_to_parallel_done[i] != -1) {
			dof += 1;

		}
		//if the plane is orthogonaled(no coplanar), its sign distance dof is 1.
		else if (planes_to_coplanar_done[i] == -1 && planes_to_orthogonal_done[i] != -1 && planes_to_parallel_done[i] == -1) {
			dof += 1;

		}

	}
	//add the orientation dof of orthogonal groups(including coplanar cases).
	dof += 3 * orthogonal_done_to_planes.size();
	//add the orientation dof of parallel groups(including coplanar cases).
	dof += 2 * parallel_done_to_planes.size();
	//add the sign distance dof of coplanar groups.
	dof +=  coplanar_done_to_planes.size();

	return dof;
}

//initialize the vectors that record the conducted regularity operations (regularized primitives).
void Shape_Detector::initial_regular_groups_done() {
	planes_if_regularized = std::vector<bool>(planes_2.size(), false);
	planes_to_parallel_done = std::vector<int>(planes_2.size(), -1);
	planes_to_orthogonal_done = std::vector<int>(planes_2.size(), -1);
	planes_to_coplanar_done = std::vector<int>(planes_2.size(), -1);
	parallel_done_to_planes.clear();
	orthogonal_done_to_planes.clear();
	coplanar_done_to_planes.clear();


}


//**************

//**************

//**************

//************Functions for calculating the energy changing for 'any' operation.

//energy changing for regularity operation (L2).
double Shape_Detector::energy_changed_regularization_euclidean_distance(double dis, double changed_freedom) {

	double maxixum_freedom = 3.0* double(points.size()) / double(min_points);

	if (weight_mode == 2) {


		double term_distance = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (mean_distance_current);
		double term_freedom = lambda_regular * changed_freedom / (double(freedom_of_planes_current));


		return(term_distance + term_freedom);
	}
	else if (weight_mode == 1) {
		double term_distance = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (ori_mean_error);
		double term_freedom = lambda_regular * changed_freedom / (double(ori_freedom_of_planes));

		return(term_distance + term_freedom);

	}
	else {
		double term_distance = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (epsilon);

		double term_freedom = lambda_regular * changed_freedom / (double(maxixum_freedom));

		return(term_distance + term_freedom);


	}


}

double Shape_Detector::energy_changed_regularization_normal_deviation(double dis, double changed_freedom) {

	double maxixum_freedom = 3.0* double(points.size()) / double(min_points);

	if (weight_mode == 2) {


		double term_distance = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (mean_normal_current);


		double term_freedom = lambda_regular * changed_freedom / (double(freedom_of_planes_current));


		return(term_distance + term_freedom);
	}
	else if (weight_mode == 1) {
		double term_distance = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (ori_mean_error);


		double term_freedom = lambda_regular * changed_freedom / (double(ori_freedom_of_planes));

		return(term_distance + term_freedom);

	}
	else {
		double term_distance = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (epsilon);

		double term_freedom = lambda_regular * changed_freedom / (double(maxixum_freedom));

		return(term_distance + term_freedom);


	}


}

double Shape_Detector::energy_changed_second_normal(double dis, double numb, int id_plane) {
	double change_dof = double(change_of_degrees_of_freedom_after_add_remove(id_plane));
	if (weight_mode == 2) {
		double change_fedilite = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (mean_normal_current);


		double change_completness = double(lambda_c)*numb / (double(number_inlier_before_opers));
		
		double change_freedom = double(lambda_regular)*change_dof / (double(freedom_of_planes_current));


		return (change_fedilite + change_completness + change_freedom);
	}
	else if (weight_mode == 1) {
		double change_fedilite = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (ori_mean_normal_diviation);

	

		double change_completness = double(lambda_c)*numb / (double(ori_inliers_number));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));


		return (change_fedilite + change_completness + change_freedom);
	}
	else {
		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double change_fedilite = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (1.0);

		
		double change_completness = double(lambda_c)*numb / (double(points.size()));
		
		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));


		return (change_fedilite + change_completness + change_freedom);

	}


}

double Shape_Detector::energy_changed_normal(double dis, double numb, int id_plane) {

	double change_dof = double(change_of_degrees_of_freedom_after_split(id_plane));

	if (weight_mode == 2) {
		double term1 = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (mean_normal_current);

		

		double term2 = lambda_r * numb / (double(size_current_primitives));
		
		double change_freedom = double(lambda_regular)*(change_dof) / (double(freedom_of_planes_current));
		return(term1 + term2 + change_freedom);
	}
	else if (weight_mode == 1) {
		double term1 = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (ori_mean_normal_diviation);


		double term2 = lambda_r * numb / (double(ori_primitives_number));
		
		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));
		return(term1 + term2 + change_freedom);
	}
	else {

		double fix_number_prim = double(points.size()) / double(min_points);
		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double term1 = -double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (1.0);




		double term2 = lambda_r * numb / double(fix_number_prim);
		
		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));

		return(term1 + term2 + change_freedom);
	}

}

double Shape_Detector::energy_changed_normal_merge(double dis, double numb, int nmoves, int id_plane_1, int id_plane_2) {
	double change_dof = double(change_of_degrees_of_freedom_after_mergy(id_plane_1, id_plane_2));
	double fix_number_prim = double(points.size()) / double(min_points);

	if (weight_mode == 2) {
		double term1 = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / (mean_normal_current);
		double term2 = lambda_r * numb / (double(size_current_primitives));
		double term3 = double(lambda_c)*nmoves / (double(number_inlier_before_opers));
		double change_freedom = double(lambda_regular)*(change_dof) / (double(freedom_of_planes_current));
		return(term1 + term2 + term3 + change_freedom);
	}
	else if (weight_mode == 1) {
		double term1 = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / (ori_mean_normal_diviation);
		double term2 = lambda_r * numb / (double(ori_primitives_number));
		double term3 = double(lambda_c)*nmoves / (double(ori_inliers_number));
		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));
		return(term1 + term2 + term3 + change_freedom);
	}
	else {
		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double term1 = -double(lambda_fidelity)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / 1.0;
		double term2 = lambda_r * numb / double(fix_number_prim);
		double term3 = double(lambda_c)*nmoves / (double(points.size()));
		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));
		return(term1 + term2 + term3 + change_freedom);
	}

}

double Shape_Detector::energy_changed_merge(double dis, double numb, int nmoves, int id_plane_1, int id_plane_2) {
	double change_dof = double(change_of_degrees_of_freedom_after_mergy(id_plane_1, id_plane_2));

	double dem_number_shape = double(points.size()) / double(min_points);

	if (weight_mode == 2) {
		double term1 = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (mean_distance_current);


		double term2 = lambda_r * numb / (double(size_current_primitives));
		double term3 = double(lambda_c)*nmoves / (double(number_inlier_before_opers));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(freedom_of_planes_current));


		return(term1 + term2 + term3 + change_freedom);
	}
	else if (weight_mode == 1) {
		double term1 = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (ori_mean_error);


		double term2 = lambda_r * numb / (double(ori_primitives_number));
		double term3 = double(lambda_c)*nmoves / (double(ori_inliers_number));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));


		return(term1 + term2 + term3 + change_freedom);

	}
	else {

		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double term1 = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (epsilon);

		double term2 = lambda_r * numb / (double(dem_number_shape));
		double term3 = double(lambda_c)*nmoves / (double(points.size()));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));

		return(term1 + term2 + term3 + change_freedom);

	}
}
double Shape_Detector::energy_changed(double dis, double numb, int id_plane) {

	double change_dof = double(change_of_degrees_of_freedom_after_split(id_plane));

	if (weight_mode == 2) {
		double term1 = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (mean_distance_current);


		double term2 = lambda_r * numb / (double(size_current_primitives));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(freedom_of_planes_current));

		return(term1 + term2 + change_freedom);
	}
	else if (weight_mode == 1) {
		double term1 = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (ori_mean_error);

		double term2 = lambda_r * numb / (double(ori_primitives_number));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));

		return(term1 + term2 + change_freedom);

	}
	else {
		double dem_number_shape = double(points.size()) / double(min_points);
		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double term1 = double(lambda_fidelity)*(dis / double(number_of_assigned_points)) / (epsilon);

		double term2 = lambda_r * numb / (double(dem_number_shape));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));

		return(term1 + term2 + change_freedom);
	}

}

double Shape_Detector::energy_changed_second(double dis, double numb, int id_plane) {
	double change_dof = double(change_of_degrees_of_freedom_after_add_remove(id_plane));

	if (weight_mode == 2) {
		double change_fedilite = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (mean_distance_current);



		double change_completness = double(lambda_c)*(numb / (double(number_inlier_before_opers)));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(freedom_of_planes_current));




		return (change_fedilite + change_completness + change_freedom);
	}
	else if (weight_mode == 1) {
		double change_fedilite = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (ori_mean_error);



		double change_completness = double(lambda_c)*(numb / (double(ori_inliers_number)));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(ori_freedom_of_planes));

		return (change_fedilite + change_completness + change_freedom);
	}
	else {
		double maxixum_freedom = 3.0* double(points.size()) / double(min_points);
		double change_fedilite = double(lambda_fidelity)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (epsilon);



		double change_completness = double(lambda_c)*(numb / (double(ori_inliers_number)));

		double change_freedom = double(lambda_regular)*(change_dof) / (double(maxixum_freedom));

		return (change_fedilite + change_completness + change_freedom);
	}


}

//***********Functions for refining planar shape configurations with dfferent metrics: L2, L1.1, Hybrid.
//In each function, local operations(merging, splitting, insertion and exclusion) 
//and global transfert operator are alternatively conducted until 
//the energy can not be reduced anymore.

//L2.
void Shape_Detector::planar_shape_detection_l2() {

	int timett = 0;//time of loop
	bool propagation2;

	//get the current informations
	get_last_information();
	old_size_current_primitives = last_primitives_number;
	old_mean_normal_diaviation = last_normal_deviation;
	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	old_freedom_of_planes = last_freedom_planes;

	//a vector to show how many times each plane has been splitted.
	region_type = std::vector<int>(planes_to_inliers.size(), 0);

	//a vector to show how many times each point has been inserted to a plane.
	points_if_added = std::vector<int>(points.size(), -1);

	//number of operations that have been conducted.
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	all_regularization = 0;
	clock_t t_start_all = clock();

	//In each loop, first: global transder, then local operations.
    number_iterations=0;
	do {
		timett++;
		propagation2 = true;

		//number of maximum points that could be considered in each insertion and exclusion operation. 
		number_of_insert_exclude *= ((timett) / 7 + 1);

		clock_t t_start = clock();
		//get the mean distance before transfer operator.
		mean_distance_diaviation_before_transfer = mean_distance_diaviation;

		//do the transfer operator and take L2 as the criterion.
		transfer_operator_l2();

		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		
		//show the information of the current configuration.
		get_distance_diviation_show_normal_info(trsfer_time);

		//Before 7 times, transfer operator does not skip regulared primitives, so we think that

		//the regularized primitives are affected by trnasfer operator, we initial the primitive relartionships.
		if (all_t_transfer <= 7) {
			//initialize the vectors that record the conducted regularity operations (regularized primitives).
			initial_regular_groups_done();
            calculate_center();//get the center of each plane
		}
		all_t_transfer++;


		clock_t t_start1 = clock();
		//vector: save all the insertion and exclusion operations and their energy changing.
		good_points_shape.clear();
		bad_points_shape.clear();
		//get the information before local operations.
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;
		if (lambda_c > 0) {
			//find all the insertion operations and thier energy changing.
			get_good_points();
		}
		if (timett > 1) {
			//find all the exclusion operations and thier energy changing.
			get_bad_points();
		}
		//find all the regularity operations and their energy changing.
		calculate_energy_changing_for_regularization_initialization();
		//find all the adjancent pairs of primitives.
		test_connected_primitives();
		//conduct the local operations under the guidance of a dynamic energy.
		local_operators();

		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;

		//get the information after local operations.
		get_distance_diviation_show_merge_info(local_time);

		number_iterations++;

        if ((t_m == 0 && timett > 1)||(t_m==t_regularization&&timett>10)
                || timett > 800 || (timett > 500 && t_split == t_merge) || if_stop
                || (max_iter>0 && number_iterations > max_iter))
            propagation2 = false;
	} while (propagation2);

	orthogonal_numbers = 0;
	parallel_numbers = 0;
	coplanar_numbers = 0;

	get_distance_diviation();
	clock_t t_end_all = clock();
	interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

	get_coverage_and_mean_error_pure();
	planes_1 = planes_2;
	planes_0 = planes_2;


}

//L1.1.
void Shape_Detector::planar_shape_detection_L1() {

	int timett = 0;
	bool propagation2;

	get_last_information();

	region_type = std::vector<int>(planes_to_inliers.size(), 0);

	points_if_added = std::vector<int>(points.size(), -1);

	old_size_current_primitives = last_primitives_number;
	old_mean_normal_diaviation = last_normal_deviation;
	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	old_freedom_of_planes = last_freedom_planes;
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	all_regularization = 0;
	clock_t t_start_all = clock();
    number_iterations=0;
	do {
		timett++;

		number_of_insert_exclude *= ((timett) / 7 + 1);
		propagation2 = true;

		clock_t t_start = clock();


		//get_distance_diviation();
		mean_normal_diviation_befor_transfer = mean_normal_current;
		transfer_operator_normal();



		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		
		get_distance_diviation_show_normal_info(trsfer_time);
		if (all_t_transfer <= 7) {
			
			initial_regular_groups_done();
            calculate_center();//get the center of each plane
		}
		all_t_transfer++;



		clock_t t_start1 = clock();

		good_points_shape.clear();
		bad_points_shape.clear();
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;
		if (lambda_c > 0) {
			get_good_points_normal();
		}
		//the first cycle, we do not exclude the point, since the initialized configuration is bad, there are plenty of points would be excluded befor transfer.
		//And the excluded points will be inserted later. 
		if (timett > 1) {
			get_bad_points_normal();
		}
		

		calculate_energy_changing_for_regularization_initialization_normal();
		

		test_connected_primitives();
		local_operators_normal();
	


		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;

		get_distance_diviation_show_merge_info(local_time);

		number_iterations++;

        if ((t_m == 0 && timett > 1)||(t_m==t_regularization&&timett>10)
                || timett > 800 || (timett > 500 && t_split == t_merge) || if_stop
                || (max_iter>0 && number_iterations > max_iter))
            propagation2 = false;
	} while (propagation2);

	clock_t t_end_all = clock();
	interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

	get_distance_diviation();


	get_coverage_and_mean_error_pure();
	planes_1 = planes_2;
	planes_0 = planes_2;


}

//Hybrid: L2 in energy function, L1.1 in global transfert operator.
void Shape_Detector::planar_shape_detection_hybrid() {

	//test_curvature();
	int timett = 0;
	bool propagation2;
	//SD->set_primitives(true);
	get_last_information();


	region_type = std::vector<int>(planes_to_inliers.size(), 0);

	points_if_added = std::vector<int>(points.size(), -1);

	old_size_current_primitives = last_primitives_number;
	old_mean_normal_diaviation = last_normal_deviation;
	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	old_freedom_of_planes = last_freedom_planes;

	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	all_regularization = 0;
	clock_t t_start_all = clock();
    number_iterations=0;
    do{
		timett++;
		propagation2 = true;

		number_of_insert_exclude *= ((timett) / 7 + 1);

		clock_t t_start = clock();

		mean_distance_diaviation_before_transfer = mean_distance_diaviation;

		transfer_operator_normal_for_hybrid();

		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		
		get_distance_diviation_show_normal_info(trsfer_time);
		if (all_t_transfer <= 7) {
			
			initial_regular_groups_done();
            calculate_center();//get the center of each plane

		}
		all_t_transfer++;
		clock_t t_start1 = clock();
		good_points_shape.clear();
		bad_points_shape.clear();
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;

		if (lambda_c > 0) {
			get_good_points();
		}
		if (timett > 1) {
			get_bad_points();

		}
		
		calculate_energy_changing_for_regularization_initialization();
		test_connected_primitives();
		local_operators();


		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;

		get_distance_diviation_show_merge_info(local_time);
		
		number_iterations++;

        if ((t_m == 0 && timett > 1)||(t_m==t_regularization&&timett>10)
                || timett > 800 || (timett > 500 && t_split == t_merge) || if_stop
                || (max_iter>0 && number_iterations > max_iter))
            propagation2 = false;
    }while (propagation2);


    orthogonal_numbers=0;
    parallel_numbers=0;
    coplanar_numbers=0;

    get_distance_diviation();
    clock_t t_end_all = clock();
    interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

    get_coverage_and_mean_error_pure();
    planes_1 = planes_2;
    planes_0 = planes_2;

}

//************Functions for conducting transfer operators with different metrics: L2, L1.1
//Lloyd's algorithm is used to solve the optimization problem.
//In each epoch, the boundary points of adjacent primitives are moved to a better primitive.
//Better means: smaller distance error (L2) or smaller normal diviation(L1.1).

//transfert: L2, energy: L2
void Shape_Detector::transfer_operator_l2() {

	bool propagation1;
	int times = 0;
	int number_moving;

	t_l = 0;
	//A vector to count the transfered times of each point.
	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);
	//A vector to record the point number of each primitive.
	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back(planes_to_inliers[mm].size());
	}
	//K used for finding neigbors.
	int Nb_neigh = knn;
	do {
		//backup.
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;
		std::vector<int> region_type_last = region_type;
		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;

		times++;
		number_moving = 0;
		propagation1 = true;

		//For each point, first check if it's boundary point, then transfer it to the better primitive.

		for (size_t i = 0; i < points.size(); ++i) {
			//if this point is outlier, skip it.
			if (inliers_to_planes[i] == -1) continue;

			//if this point has been transfered more than 5 times, skip it.
			if (if_moved[i] >= 5) continue;

			//if this point belongs to a primitive that has too small inlier number, skip it.
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;

			//skip regularized primitives after 7 iterations.
			if (planes_if_regularized[inliers_to_planes[i]] && all_t_transfer > 7) continue;


			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}
			std::set<int> one_neight_id;

			//check if this point is on the boundary between primitives
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];
				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {
					//record the neigbor primitives' id.
					one_neight_id.insert(inliers_to_planes[neighbor_index]);
				}
			}
			//if it does not have neighbor primitive, then skip it.
			if (one_neight_id.empty()) { continue; }

			double min_distance = sqrt(CGAL::squared_distance(planes_2[inliers_to_planes[i]], points[i].first));

			int changed_plane_id = -1;
			//transfer the point to its best primitive, best means: smallest distance error and normal diviation is also under threshold.
			for (int neight_id : one_neight_id) {

				double neight_distance = sqrt(CGAL::squared_distance(planes_2[neight_id], points[i].first));

				if (neight_distance < min_distance) {

					if (abs(points[i].second * planes_2[neight_id].orthogonal_vector()) >= normal_threshold) {
						min_distance = neight_distance;
						changed_plane_id = neight_id;
					}
				}
			}
			if (changed_plane_id != -1) {
				//When the trnasfer operator is conducted more than 7 times, the inliers of regularized primitives are not trnsferred.
				if (!planes_if_regularized[changed_plane_id] || all_t_transfer <= 7) {
					//conduct the transfer and update the inlier number of impacted primitives.
					update_primitive_size[inliers_to_planes[i]]--;
					inliers_to_planes[i] = changed_plane_id;
					update_primitive_size[changed_plane_id]++;
					//update the trnasfer number of this point.
					if_moved[i] ++;
					number_moving++;
				}
			}
		}
		//update planes_to_inliers.
		int plane_size = planes_to_inliers.size();
		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}

		//update planes_2.
		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {

			if (!planes_if_regularized[ii] || all_t_transfer <= 7) {
				std::vector<Inexact_Point_3> inliers_i;
				inliers_i.reserve(planes_to_inliers[ii].size());

				for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
					const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
					inliers_i.push_back(pt);

				}



				if (inliers_i.size() < 3) {
					getchar();
				}
				Inexact_Plane plane;
				linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());

				planes_2[ii]=plane;
			}

		}
		//get the informations of the current configuration.
		get_distance_diviation();
		//after 10 times of transfer, make sure that the transfer operation reducing the energy, which can avoid the infinity loop between local and global operations.
		if (mean_distance_diaviation > mean_distance_diaviation_before_transfer&&all_t_transfer > 10) {
			//if mean distance improved, go back to the configuration before transfer.
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;
			region_type = region_type_last;
			inliers_to_planes = inliers_to_planes_last;
		}

		//if constraint is set, the mean distance should be smaller than the original one (the configuration of region growing).
		if (if_constraint) {

			if (mean_distance_diaviation > ori_mean_error) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;
				region_type = region_type_last;
				inliers_to_planes = inliers_to_planes_last;

			}
		}
		t_l++;

		if (number_moving == 0 || times > 100) { propagation1 = false; }
	} while (propagation1);

}
//transfert: L1.1, energy: L2
//******* comment please check transfer_operator_l2()
void Shape_Detector::transfer_operator_normal_for_hybrid() {
	
	bool propagation1;
	int times = 0;
	int number_moving;

	t_l = 0;


	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);


	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back((int)planes_to_inliers[mm].size());
	}


	int Nb_neigh = knn;
	do {
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;
		std::vector<int> region_type_last = region_type;


		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;

		times++;
		number_moving = 0;
		propagation1 = true;

		for (size_t i = 0; i < points.size(); ++i) {

			if (inliers_to_planes[i] == -1) continue;
			if (if_moved[i] >= 5) continue;
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;
			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}
			if (planes_if_regularized[inliers_to_planes[i]] && all_t_transfer > 7) continue;
			std::set<int> one_neight_id;//connected primitives for the point

			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];


				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {

					one_neight_id.insert(inliers_to_planes[neighbor_index]);
				}
			}
			if (one_neight_id.empty()) { continue; }

			Inexact_Point_3 this_point = points[i].first;


			//detect the best primitive
			double max_cos = abs(points[i].second * planes_2[inliers_to_planes[i]].orthogonal_vector());

			int changed_plane_id = -1;
			for (int neight_id : one_neight_id) {


				double that_cos = abs(points[i].second * planes_2[neight_id].orthogonal_vector());
				if (that_cos > max_cos) {
					if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_point)) <= epsilon) {

						max_cos = that_cos;
						changed_plane_id = neight_id;
					}
				}

			}
			if (changed_plane_id != -1) {
				if (!planes_if_regularized[changed_plane_id] || all_t_transfer <= 7) {

					update_primitive_size[inliers_to_planes[i]]--;
					inliers_to_planes[i] = changed_plane_id;

					update_primitive_size[changed_plane_id]++;
					if_moved[i]++;
					number_moving++;
				}
			}

		}
		//update the configuration
		int plane_size = planes_to_inliers.size();

		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}

		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {
			if (!planes_if_regularized[ii] || all_t_transfer <= 7) {

				std::vector<Inexact_Point_3> inliers_i;
				inliers_i.reserve(planes_to_inliers[ii].size());

				for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
					const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
					inliers_i.push_back(pt);

				}


				if (inliers_i.size() < 3) {


					getchar();
				}
				Inexact_Plane plane;
				linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());


				planes_2[ii] = plane;
			}

		}

		get_distance_diviation();
		//make sure that the transfer operation reducing the energy.

		if (mean_distance_diaviation > mean_distance_diaviation_before_transfer&&all_t_transfer > 10) {
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;
			region_type = region_type_last;
			inliers_to_planes = inliers_to_planes_last;


		}

		if (if_constraint) {

			if (mean_distance_diaviation > ori_mean_error) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;
				region_type = region_type_last;

				inliers_to_planes = inliers_to_planes_last;
			
			}
		}
		t_l++;

		if (number_moving == 0 || times > 100) { propagation1 = false; }
	} while (propagation1);

}

//transfert: L1.1, energy: L1.1
//******* comment please check transfer_operator_l2()
void Shape_Detector::transfer_operator_normal() {

	bool propagation1;
	int times = 0;
	int number_moving;

	t_l = 0;

	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);

	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back((int)planes_to_inliers[mm].size());
	}

	int Nb_neigh = knn;
	do {
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;
		std::vector<int> region_type_last = region_type;

		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;

		times++;
		number_moving = 0;
		propagation1 = true;

		for (size_t i = 0; i < points.size(); ++i) {
			if (inliers_to_planes[i] == -1) continue;
			if (if_moved[i] >= 50) continue;
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;
			if (planes_if_regularized[inliers_to_planes[i]] && all_t_transfer > 7) continue;
			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}

			std::set<int> one_neight_id;//connected primitives for the point
			//check if this point is on the boundary of primitives.
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];


				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {

					one_neight_id.insert(inliers_to_planes[neighbor_index]);

				}


			}

			if (one_neight_id.empty()) { continue; }

			Inexact_Point_3 this_point = points[i].first;


			//detect the best primitive
			double max_cos = abs(points[i].second * planes_2[inliers_to_planes[i]].orthogonal_vector());

			int changed_plane_id = -1;
			for (int neight_id : one_neight_id) {


				double that_cos = abs(points[i].second * planes_2[neight_id].orthogonal_vector());
				if (that_cos > max_cos) {
					if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_point)) <= epsilon) {

						max_cos = that_cos;
						changed_plane_id = neight_id;
					}
				}

			}

			if (changed_plane_id != -1) {
				if (!planes_if_regularized[changed_plane_id] || all_t_transfer <= 7) {
					update_primitive_size[inliers_to_planes[i]]--;
					inliers_to_planes[i] = changed_plane_id;
					update_primitive_size[changed_plane_id]++;
					if_moved[i]++;
					number_moving++;
				}

			}

		}
		//update the configuration
		int plane_size = planes_to_inliers.size();

		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}

		



		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {

			if (!planes_if_regularized[ii] || all_t_transfer <= 7) {
				std::vector<Inexact_Point_3> inliers_i;
				inliers_i.reserve(planes_to_inliers[ii].size());

				for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
					const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
					inliers_i.push_back(pt);

				}


				if (inliers_i.size() < 3) {

					getchar();
				}
				Inexact_Plane plane;
				linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());


				planes_2[ii]=plane;
			}
		}

		get_distance_diviation();
		//make sure that the transfer operation reducing the energy.
		if (mean_normal_diviation < mean_normal_diviation_befor_transfer && all_t_transfer>10) {
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;
			region_type = region_type_last;
			inliers_to_planes = inliers_to_planes_last;


		}

		if (if_constraint) {

			if (mean_normal_diviation < ori_mean_normal_diviation) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;
				region_type = region_type_last;
				inliers_to_planes = inliers_to_planes_last;


			}
		}


		t_l++;

		if (number_moving == 0 || times > 1000) { propagation1 = false; }
	} while (propagation1);







}


//**************

//**************

//**************

//**************

//**************

//**************


//**********Function for local operations: splitting, merging, insertion, exclusion and regularization.

//Finding all the insertion operation when fidelity metric is L2.

void Shape_Detector::get_good_points() {
	if_insert_candidate = std::vector<int>(points.size(), -1);
	good_points_shape.clear();
	int Nb_neigh = 10;
	std::map<int, std::vector<std::pair<int, double>>> plane_outliers;
	//For each outlier, find its best neighbor primitives. Best: closest and normal deviation under a threshold.
	for (int i = 0; i < points.size(); ++i) {
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;

		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;



		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);

			}


		}

		if (one_neight_id.empty()) { continue; }

		double mm_dis = epsilon;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {

			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > normal_threshold) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < mm_dis) {

					mm_dis = sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first));
					changed_plane_id = neight_id;

				}
			}

		}

		if (changed_plane_id != -1) {
			plane_outliers[changed_plane_id].push_back(std::make_pair(i, mm_dis));//plane_id, point_id, distance.

		}

	}
	std::map<int, std::vector<std::pair<int, double>>>::iterator g_i;
	g_i = plane_outliers.begin();

	//For each primitive, get its 'number_of_insert_exclude' insertable outliers. And calculate the disatnce error changing and energy changing for each insertion operation.

	//Note that an insertion operation is conducted on one primitive and several outliers.
	while (g_i != plane_outliers.end()) {
		std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator> p_q;

		for (std::pair<int, double> ppp : g_i->second) {
			p_q.push(ppp);
		}

		std::vector<int> added_points_ids;
		double en = 0;
		double dif_e = 0;
		double n_a = 0;
		while (!p_q.empty() && n_a < number_of_insert_exclude) {
			added_points_ids.push_back(p_q.top().first);
			p_q.pop();
			n_a++;
		}
		dif_e = add_changed_error(g_i->first, added_points_ids);
		en = energy_changed_second(dif_e, -n_a, g_i->first);
		if (added_points_ids.size() != 0 && en < 0) {
			for (int idd : added_points_ids) {
				if_insert_candidate[idd] = 0;
			}
			std::vector<double> dd;
			dd.push_back(en);
			dd.push_back(dif_e);
			good_points_shape.push_back(std::make_pair(std::make_pair(g_i->first, added_points_ids), dd));
		}
		g_i++;
	}

}

//Finding all the exclusion operation when fidelity metric is L2.
void Shape_Detector::get_bad_points() {

	bad_points_shape.clear();
	if (old_coverage > ori_coverage) {
		//For each primitive, find its father inliers. farther: distance more than 1.5 * mean_distance_diaviation.

		std::map<int, std::vector<std::pair<int, double>>> plane_inliers;
		double bigest_dis = 1.5 * mean_distance_diaviation;

		for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

			const Inexact_Plane & H = planes_2[i];


			for (int j : planes_to_inliers[i]) {
				const Point_with_normal & pt = points[j];
				if (points_if_added[j] > 2) continue;
				if (sqrt(CGAL::squared_distance(H, pt.first)) > bigest_dis) {
					plane_inliers[i].push_back(std::make_pair(j, sqrt(CGAL::squared_distance(H, pt.first))));
				}

			}
		}

		std::map<int, std::vector<std::pair<int, double>>>::iterator b_o;
		b_o = plane_inliers.begin();
		//For each primitive, get its 'number_of_insert_exclude' fathest inliers. And calculate the disatnce error changing and energy changing for each exclusion operation.

		//Note that an exclusion operation is conducted on one primitive and several inliers.
		while (b_o != plane_inliers.end()) {


			std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator> p_q;
			for (std::pair<int, double> ppp : b_o->second) {
				p_q.push(ppp);

			}
			std::vector<int> removed_points_ids;
			double en = 0;
			double dif_e = 0;

			double n_b = 0;
			while (!p_q.empty() && n_b < number_of_insert_exclude) {
				removed_points_ids.push_back(p_q.top().first);


				p_q.pop();
				n_b++;
			}
			dif_e = remove_changed_error(b_o->first, removed_points_ids);
			en = energy_changed_second(dif_e, n_b, b_o->first);

			if (removed_points_ids.size() != 0 && (planes_to_inliers[b_o->first].size() - removed_points_ids.size()) >= min_points && en < 0) {
				std::vector<double> dd;
				dd.push_back(en);
				dd.push_back(dif_e);
				bad_points_shape.push_back(std::make_pair(std::make_pair(b_o->first, removed_points_ids), dd));
			}
			b_o++;
		}

	}
}

//Finding all the regularity operation when fidelity metric is L2.
void Shape_Detector::calculate_energy_changing_for_regularization_initialization() {
	//find all the nearly parallel primitives using mean shift. We denote each nearly paralle cluster by 'parallel cluster'.
	initialize_parallel_clusters();

	//find all the nearly orthogonal 'parallel cluster'.
	initialize_orthogonal_clusters_after_parallel();

	//find all the nearly coplanar primitives in each 'parallel cluster'.
	initialize_coplanar_clusters_after_parallel();

	

	//Note that, the regular relationships are dependant and complicated, thus we do not conduct them individually.

	//We make the following constraint:

	//If we make the primitives of a 'parallel cluster' to be exactly parallel, we should also make the primitivs of its nearly orthogonal 'paralle cluster's to be exactly paralle

	//and make the 'parallel cluster's be exact orthogonal. The nearly coplanar primitives in these 'parallel cluster's should also be made exactly coplanar.

	energy_changing_initial_regularization.clear();//save energy changing and euclidean changing
	energy_changing_initial_regularization.resize(parallel_clusters_to_planes.size());
	if_regularization_can_conducted.clear();
	if_regularization_can_conducted = std::vector<bool>(parallel_clusters_to_planes.size(), true);

	//We enumerate all the situations of each 'parallel cluster' and find all the possible regularity operations and the energy changing.

	//Note that the regularity operations are conducted on several primitives.
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		

		if (parallel_clusters_to_planes[i].size() == 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 && parallel_cluster_to_coplanr_cases[i].size() == 0) {
			//nothing happend.
			energy_changing_initial_regularization[i].push_back(10);
			energy_changing_initial_regularization[i].push_back(0);

		}
		else if (parallel_clusters_to_planes[i].size() > 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 ) {
			// parallel or paralle + coplanar
			double euclidean_distance_change = euclidean_distance_change_after_a_coplanar_or_parallel(parallel_clusters_to_planes[i],i);
			double degree_freedom_change = get_changed_degrees_of_freedom_after_a_regular_operation(i);
			double energy_change = energy_changed_regularization_euclidean_distance(euclidean_distance_change, degree_freedom_change);
			energy_changing_initial_regularization[i].push_back(energy_change);
			energy_changing_initial_regularization[i].push_back(euclidean_distance_change);


		}
		else if ( parallel_id_changed_normal_after_orthogonal[i].size() > 0 ) {
		//orthogonal or orthogonal + coplanar. Note that parallel is a part of orthogonal.

			double euclidean_distance_change = euclidean_distance_change_after_orthogonal(parallel_id_changed_normal_after_orthogonal[i]);
			double degree_freedom_change = get_changed_degrees_of_freedom_after_a_regular_operation(i);

			double energy_change = energy_changed_regularization_euclidean_distance(euclidean_distance_change, degree_freedom_change);
			energy_changing_initial_regularization[i].push_back(energy_change);
			energy_changing_initial_regularization[i].push_back(euclidean_distance_change);
		}
		else {
			std::cout << "wrong situation!!!!!!!!!!!" << std::endl;
			getchar();

		}

	}





}


//Finding all the regularity operation when fidelity metric is L1.1.
void Shape_Detector::calculate_energy_changing_for_regularization_initialization_normal() {

	initialize_parallel_clusters();



	initialize_orthogonal_clusters_after_parallel();




	initialize_coplanar_clusters_after_parallel();




	//todo:energy
	energy_changing_initial_regularization.clear();
	energy_changing_initial_regularization.resize(parallel_clusters_to_planes.size());
	if_regularization_can_conducted.clear();
	if_regularization_can_conducted = std::vector<bool>(parallel_clusters_to_planes.size(), true);
	
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		if (parallel_clusters_to_planes[i].size() == 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 && parallel_cluster_to_coplanr_cases[i].size() == 0) {
			//nothing can changed
			//std::cout << "nothing can changed" << std::endl;
			energy_changing_initial_regularization[i].push_back(10);
			energy_changing_initial_regularization[i].push_back(0);

		}
		else if (parallel_clusters_to_planes[i].size() > 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0) {
			// parallel or parallel + coplanar.
			//std::cout << "only parallel operator." << std::endl;
			double nomral_deviation_change = normal_deviation_change_after_a_coplanar_or_parallel(parallel_clusters_to_planes[i],i);
			double degree_freedom_change = get_changed_degrees_of_freedom_after_a_regular_operation(i);
			double energy_change = energy_changed_regularization_normal_deviation(nomral_deviation_change, degree_freedom_change);
			energy_changing_initial_regularization[i].push_back(energy_change);
			energy_changing_initial_regularization[i].push_back(nomral_deviation_change);


		}
		else if (parallel_id_changed_normal_after_orthogonal[i].size() > 0) {
			// orthogonal or orthogonal + coplanar.
			//std::cout << "only orthogonal operator." << std::endl;

			double nomral_deviation_change = normal_deviation_change_after_orthogonal(parallel_id_changed_normal_after_orthogonal[i]);
			double degree_freedom_change = get_changed_degrees_of_freedom_after_a_regular_operation(i);
			double energy_change = energy_changed_regularization_normal_deviation(nomral_deviation_change, degree_freedom_change);
			energy_changing_initial_regularization[i].push_back(energy_change);
			energy_changing_initial_regularization[i].push_back(nomral_deviation_change);
		}
		
		
		else {
			std::cout << "wrong situation!!!!!!!!!!!" << std::endl;
			getchar();

		}

	}


}


//**************

//**************

//**************

//**************

//**************

//**************Functions for conducting the local functions with differnt fidelity metrics: L2, L1.1.


//Conduct the local operations under the guidance of a dynamic priority queue with fidelity metric: L2.

//All the local operations are first put into a priority queue and sorted in the descending order of their energy changing.

//Conduct the first operation in the priority queue if it reduces the energy function and under some constraint.

//Update the impacted local operations.

//Loop until there is no more satisfied operation that can reduce the energy.



//****Fidelity metirc: L2.
void Shape_Detector::local_operators() {
	//points_changed = std::vector<int>(points.size(), 0);
	bool bb = false;
	all_t_merge += t_merge;
	all_t_split += t_split;
	all_t_insert += t_insert;
	all_t_exlude += t_exlude;
	all_regularization += t_regularization;
	t_m = 0;
	t_merge = 0;
	t_split = 0;
	t_exlude = 0;
	t_insert = 0;
	t_regularization = 0;




	
	std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy> p_q;

	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> v_p_q;

	std::vector<std::vector<int>> max_list_vector;
	std::vector<std::vector<int>> min_list_vector;

	//Put all the merging operations to the priority queue.
	for (int i = 0; i < primitive_connection.size(); ++i) {

		std::vector<int> one = primitive_connection[i];
		for (int j = i + 1; j < one.size(); ++j) {
			if (one[j] >= 1) {


				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
				//check if i and j primitives can be merged.
				bool if_satisfied = convert_form_merge_l2(i, j, one_element);
				if (!if_satisfied) continue;

				p_q.push(one_element);
				v_p_q.push_back(one_element);

			}
		}
	}

	//Put all the splitting operations to the priority queue.

	for (int i = 0; i < planes_to_inliers.size(); ++i) {

		if (region_type[i] > 5) continue;
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		//check if i primitives can be merged.

		bool if_satisfied = convert_form_split_l2(i, max_list_vector, min_list_vector, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);




	}

	//Put all the exclusion operations to the priority queue.

	for (int i = 0; i < bad_points_shape.size(); ++i) {
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_exclude_l2(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}

	//Put all the insersion operations to the priority queue.

	for (int i = 0; i < good_points_shape.size(); ++i) {




		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_insert_l2(i, one_element);
		if (!if_satisfied) continue;

		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}

	//Put all the regularity operations to the priority queue.

	for (int i = 0; i < energy_changing_initial_regularization.size(); ++i) {
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_regularization_l2(i, one_element);
		if (!if_satisfied) continue;

		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}





	do {

		bb = true;

		std::vector<int> if_merged;
		if_merged = std::vector<int>(planes_to_inliers.size(), -1);



		std::vector<int> one_merged;

		std::vector<int> max_list;
		std::vector<int> min_list;
		std::vector<int> id_of_moved_merge;


		int one_splite = -1;
		std::vector<int> one_remove;
		int one_remove_plane = -1;
		std::vector<int> one_add;
		int one_add_plane = -1;
		std::vector<int> regularization_planes;
		int id_parallel_cluster;
		//test the top operation in the priority queue, if it's satisfied, conduct this operation and update the priority queue.
		while (!p_q.empty()) {


			if (p_q.top().second[0] >= 0) break;
			if (p_q.top().first.first == 1) {
				//merging operation
				std::vector<int> this_pair = p_q.top().first.second;

				int f = this_pair[0];
				int s = this_pair[1];





				double after_operation_distance_visual = all_distance_diaviation + p_q.top().second[1];
				//if constraint, the changed configuration should have higher fidelity, completenness and simplicity than original one.
				if (if_constraint) {
					
					if (number_of_assigned_points + 2 - this_pair.size() <= ori_inliers_number) {
						p_q.pop();
						continue;
					}

					if (after_operation_distance_visual / double(number_of_assigned_points) > ori_mean_error) {

						p_q.pop();
						continue;
					}
					else {
						p_q.pop();
						id_of_moved_merge.clear();
						if (this_pair.size() > 2) {
							for (int k = 2; k < this_pair.size(); ++k) {
								id_of_moved_merge.push_back(this_pair[k]);
							}
						}
						number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();
						if_merged[f] = 0;
						if_merged[s] = 0;
						one_merged.push_back(f);
						one_merged.push_back(s);
						//after_operation_distance = after_operation_distance_visual;
						all_distance_diaviation = after_operation_distance_visual;
						mean_distance_current = all_distance_diaviation / number_of_assigned_points;

						t_m++;
						t_merge++;
						//std::cout << "merging: "<<f<<","<<s<<" "<< planes_if_regularized[f]<<" "<< planes_if_regularized[s];

						break;
					}
				}
				else {

					

					p_q.pop();
					id_of_moved_merge.clear();
					if (this_pair.size() > 2) {
						for (int k = 2; k < this_pair.size(); ++k) {
							id_of_moved_merge.push_back(this_pair[k]);
						}
					}
					number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();

					if_merged[f] = 0;
					if_merged[s] = 0;
					one_merged.push_back(f);
					one_merged.push_back(s);

					all_distance_diaviation = after_operation_distance_visual;
					mean_distance_current = all_distance_diaviation / number_of_assigned_points;
					t_m++;
					t_merge++;


					break;
					//}
				}

			}
			else if (p_q.top().first.first == 2) {
				//splitting operation

				int splited_index = p_q.top().first.second[0];
				int info_index = p_q.top().first.second[1];


				double after_splited_distance = all_distance_diaviation + p_q.top().second[1];



				max_list = max_list_vector[info_index];
				min_list = min_list_vector[info_index];

				if (max_list.size() < min_points || min_list.size() < min_points) {
					p_q.pop();
					continue;
				}
				if (if_constraint) {

					if (after_splited_distance / double(number_of_assigned_points) > ori_mean_error || (planes_2.size() + 1) > ori_primitives_number) {

						p_q.pop();
						continue;

					}
					else {

						p_q.pop();
						if_merged[splited_index] = 0;

						one_splite = splited_index;

						all_distance_diaviation = after_splited_distance;
						mean_distance_current = all_distance_diaviation / number_of_assigned_points;

						t_m++;
						t_split++;

						break;
					}
				}
				else {

					p_q.pop();
					if_merged[splited_index] = 0;

					one_splite = splited_index;

					all_distance_diaviation = after_splited_distance;
					mean_distance_current = all_distance_diaviation / number_of_assigned_points;

					t_m++;
					t_split++;


					break;
				}

			}
			else if (p_q.top().first.first == 3) {

				//exclusion operation

				int plane_id = p_q.top().first.second[0];

				
				double after_delete_distance = all_distance_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (number_of_assigned_points + 1 - p_q.top().first.second.size() <= ori_inliers_number || after_delete_distance / double(number_of_assigned_points + 1 - p_q.top().first.second.size()) >= ori_mean_error) {

						
						p_q.pop();
						continue;
					}
				}


				if (planes_to_inliers[plane_id].size() + 1 - p_q.top().first.second.size() < min_points) {

					p_q.pop();
					continue;
				}
				else {

					if_merged[plane_id] = 0;
					one_remove_plane = plane_id;
					std::vector<int> points_ids = p_q.top().first.second;
					points_ids.erase(points_ids.begin());
					one_remove = points_ids;
					
					all_distance_diaviation = after_delete_distance;

					number_of_assigned_points = number_of_assigned_points - one_remove.size();
					
					p_q.pop();
					t_m++;
					t_exlude++;
					mean_distance_current = all_distance_diaviation / number_of_assigned_points;

					break;
				}

			}
			else if (p_q.top().first.first == 4) {
				//insertion operation.



				int plane_id = p_q.top().first.second[0];

			
				double after_add_distance = all_distance_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (after_add_distance / double(number_of_assigned_points - 1 + p_q.top().first.second.size()) >= ori_mean_error) {
						
						p_q.pop();
						continue;
					}
				}
				if (planes_to_inliers[plane_id].size() - 1 + p_q.top().first.second.size() < min_points) {

					p_q.pop();
				}
				else {

					if_merged[plane_id] = 0;
					one_add_plane = plane_id;
					std::vector<int> points_ids = p_q.top().first.second;
					points_ids.erase(points_ids.begin());
					one_add = points_ids;
				
					all_distance_diaviation = after_add_distance;

					number_of_assigned_points = number_of_assigned_points + one_add.size();
					
					mean_distance_current = all_distance_diaviation / number_of_assigned_points;

					p_q.pop();
					t_m++;
					t_insert++;

					break;
				}

			}
			else if (p_q.top().first.first == 5) {
				//regularity operation


				if (lambda_regular == 0) {
					p_q.pop();

					continue;

				}


				id_parallel_cluster = p_q.top().first.second[0];
				
				
				// the operation is skipped if it containes regularized primitives.
			
				if (!if_regularization_can_conducted[id_parallel_cluster]) {
					p_q.pop();
					continue;
				}
				

				



				double after_add_distance = all_distance_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (after_add_distance / double(number_of_assigned_points) >= ori_mean_error) {
						
						p_q.pop();
						continue;
					}
				}

				regularization_planes = p_q.top().first.second;
				regularization_planes.erase(regularization_planes.begin());
			
				all_distance_diaviation = after_add_distance;
				mean_distance_current = all_distance_diaviation / number_of_assigned_points;


				p_q.pop();
				t_m++;
				t_regularization++;

				break;


			}
		}
		//conduct the merging operation
		if (one_merged.size() != 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() == 0) {

			int last_size = planes_2.size();
			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			int id_1 = one_merged[0];
			int id_2 = one_merged[1];


			

			std::vector<int> one_merge_points_real;
			std::vector<int> one_merge_points = planes_to_inliers[id_1];
			one_merge_points.insert(one_merge_points.end(), planes_to_inliers[id_2].begin(), planes_to_inliers[id_2].end());
			for (int id : one_merge_points) {
				if (std::find(id_of_moved_merge.begin(), id_of_moved_merge.end(), id) == id_of_moved_merge.end()) {
					one_merge_points_real.push_back(id);
				}
			}

			one_merge_points = one_merge_points_real;
			for (int m_id : id_of_moved_merge) {
				inliers_to_planes[m_id] = -1;
			}


			if (one_merge_points.size() >= min_points) {

				if (id_2 < id_1) {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_2) {
							respective_planes[t] = t;
						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_1) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_1);
					region_type.erase(region_type.begin() + id_2);



					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);

					planes_2.erase(planes_2.begin() + id_1);
					planes_2.erase(planes_2.begin() + id_2);

					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_1);
					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_2);

				
					

					planes_centroids.erase(planes_centroids.begin() + id_1);
					planes_centroids.erase(planes_centroids.begin() + id_2);


					planes_if_regularized.erase(planes_if_regularized.begin() + id_1);
					planes_if_regularized.erase(planes_if_regularized.begin() + id_2);

					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_1);
					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_2);
					planes_to_coplanar_done.push_back(-1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_2);
					planes_to_parallel_done.push_back(-1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_2);
					planes_to_orthogonal_done.push_back(-1);





				}
				else {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_1) {
							respective_planes[t] = t;
						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_2) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_2);
					region_type.erase(region_type.begin() + id_1);

					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);

					planes_2.erase(planes_2.begin() + id_2);
					planes_2.erase(planes_2.begin() + id_1);

					planes_centroids.erase(planes_centroids.begin() + id_2);
					planes_centroids.erase(planes_centroids.begin() + id_1);

					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_2);
					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_1);

					
			

					planes_if_regularized.erase(planes_if_regularized.begin() + id_2);
					planes_if_regularized.erase(planes_if_regularized.begin() + id_1);

					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_2);
					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_1);
					planes_to_coplanar_done.push_back(-1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_2);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_1);
					planes_to_parallel_done.push_back(-1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_2);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_1);
					planes_to_orthogonal_done.push_back(-1);


				}

				planes_if_regularized.push_back(false);


				
				//update detected regular relationships, because id of plane has been changed.
				update_regular_relations_after_merge(id_1, id_2, respective_planes);

			



				planes_to_inliers.push_back(one_merge_points);
				planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
				planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 1));

				std::vector<Inexact_Point_3> inliers;
				inliers.reserve(one_merge_points.size());

				for (int jj = 0; jj < one_merge_points.size(); ++jj) {

					inliers.push_back(points[one_merge_points[jj]].first);

				}
				Inexact_Plane plane;
				if (inliers.size() < 3) {

					continue;
				}
				linear_least_squares_fitting_3(inliers.begin(), inliers.end(), plane, CGAL::Dimension_tag<0>());

				planes_2.push_back(plane);
				
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}
				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;

				//update the record of already regularized primitives
				update_regular_done_group_by_planes();
				//update the priority queue.
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				
				p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();

				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {//when the operator is merge.
						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//when that is the merged pair.
						{

							continue;
						}
						else if ((if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) || (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0)) {//when only one primitive is merged.



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(respective_planes[this_pair[1]], respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

						else {



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o = this_pair;
							o[0] = respective_planes[this_pair[0]];
							o[1] = respective_planes[this_pair[1]];
							one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;
						if (if_merged[this_reg[0]] == 0)//when the primitive has been merged.
						{

							continue;
						}

						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o.push_back(respective_planes[this_reg[0]]);
							o.push_back(this_reg[1]);
							one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if ( !if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_l2(this_p_p[0], one_element);

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}
				}

				//add new operator, split the new merged primitive
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
				bool if_satisfied = convert_form_split_l2(last_size - 2, max_list_vector, min_list_vector, one_element_1);
				if (if_satisfied) {
					p_q.push(one_element_1);
					local_v_p_q.push_back(one_element_1);
				}
				//add new operat, exlude and inerte for the new merged primitive
				if (all_t_transfer > 0) {
					std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
					bool if_satisfied_b = update_bad_points(last_size - 2, one_element_b);
					if (if_satisfied_b) {
						p_q.push(one_element_b);
						local_v_p_q.push_back(one_element_b);
					}
				}

				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
				bool if_satisfied_g = update_good_points(last_size - 2, one_element_g);
				if (if_satisfied_g) {
					p_q.push(one_element_g);
					local_v_p_q.push_back(one_element_g);
				}

				v_p_q.clear();
				v_p_q = local_v_p_q;


			}

		}
		//conduct splitting operation
		else if (one_merged.size() == 0 && one_splite != -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() == 0) {
			int last_size = planes_2.size();

			std::vector<int> respective_planes = std::vector<int>(last_size, -1);
			


			for (int t = 0; t < last_size; ++t) {
				if (t < one_splite) {
					respective_planes[t] = t;
				}
				else if (t == one_splite) {
					respective_planes[t] = last_size;
				}

				else {
					respective_planes[t] = t - 1;

				}


			}

		


			planes_to_inliers.erase(planes_to_inliers.begin() + one_splite);
			region_type.push_back(region_type[one_splite]);
			region_type.push_back(region_type[one_splite]);
			region_type.erase(region_type.begin() + one_splite);




			planes_2.erase(planes_2.begin() + one_splite);


			planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + one_splite);
			planes_to_coplanar_done.push_back(-1);
			planes_to_coplanar_done.push_back(-1);
			planes_to_parallel_done.erase(planes_to_parallel_done.begin() + one_splite);
			planes_to_parallel_done.push_back(-1);
			planes_to_parallel_done.push_back(-1);
			planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + one_splite);
			planes_to_orthogonal_done.push_back(-1);
			planes_to_orthogonal_done.push_back(-1);




			planes_to_inliers.push_back(max_list);
			planes_to_inliers.push_back(min_list);




			planes_if_regularized.erase(planes_if_regularized.begin() + one_splite);
			planes_if_regularized.push_back(false);
			planes_if_regularized.push_back(false);

		
			planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + one_splite);
			planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
			planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 2));

			planes_centroids.erase(planes_centroids.begin() + one_splite);
			planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
			planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 2));


			update_regular_relations_after_splite(one_splite, respective_planes);
		


			std::vector<Inexact_Point_3> max_inliers;
			max_inliers.reserve(max_list.size());

			for (int jj = 0; jj < max_list.size(); ++jj) {
				
				max_inliers.push_back(points[max_list[jj]].first);

			}
			Inexact_Plane max_plane;
			if (max_inliers.size() < 3) {

				getchar();
			}
			linear_least_squares_fitting_3(max_inliers.begin(), max_inliers.end(), max_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(max_plane);
		


			std::vector<Inexact_Point_3> min_inliers;
			min_inliers.reserve(min_list.size());

			for (int jj = 0; jj < min_list.size(); ++jj) {
				
				min_inliers.push_back(points[min_list[jj]].first);

			}
			Inexact_Plane min_plane;
			if (min_inliers.size() < 3) {

				getchar();
			}
			linear_least_squares_fitting_3(min_inliers.begin(), min_inliers.end(), min_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(min_plane);
		

			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			update_regular_done_group_by_planes();
				
			
			//update energy function
			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();

			
			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {
					std::vector<int> this_pair = one_p_q.first.second;
					//double this_weight = p_q_c.top().second;
					if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//not possible case
					{
						continue;
					}
					else if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 1, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 2, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;


							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}



					}
					else if (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0) {
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 1, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 2, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}


					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o = this_pair;
						o[0] = respective_planes[this_pair[0]];
						o[1] = respective_planes[this_pair[1]];
						one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);

					}
				}
				else if (one_p_q.first.first == 2) {
					std::vector<int> this_reg = one_p_q.first.second;
					//double this_weight = p_q_c.top().second;
					if (if_merged[this_reg[0]] == 0)
					{

						continue;
					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o.push_back(respective_planes[this_reg[0]]);
						o.push_back(this_reg[1]);
						one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 3) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 4) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 5) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if ( !if_regularization_can_conducted[this_p_p[0]]) {
						continue;
					}
					else {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						convert_form_regularization_l2(this_p_p[0], one_element);

						p_q.push(one_element);
						local_v_p_q.push_back(one_element);



					}
				}
			}

			//add new operators, splite the new two primitives.
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
			bool if_satisfied = convert_form_split_l2(planes_to_inliers.size() - 1, max_list_vector, min_list_vector, one_element_1);
			if (if_satisfied) {
				p_q.push(one_element_1);
				local_v_p_q.push_back(one_element_1);
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_2;
			bool if_satisfied_2 = convert_form_split_l2(planes_to_inliers.size() - 2, max_list_vector, min_list_vector, one_element_2);
			if (if_satisfied_2) {
				p_q.push(one_element_2);
				local_v_p_q.push_back(one_element_2);
			}
			//add new operators, exlude and insert for two new primitives.
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
				bool if_satisfied_b = update_bad_points(planes_to_inliers.size() - 2, one_element_b);
				if (if_satisfied_b) {
					p_q.push(one_element_b);
					local_v_p_q.push_back(one_element_b);
				}
			}
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
			bool if_satisfied_g = update_good_points(planes_to_inliers.size() - 2, one_element_g);
			if (if_satisfied_g) {
				p_q.push(one_element_g);
				local_v_p_q.push_back(one_element_g);
			}
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b_2;
				bool if_satisfied_b_2 = update_bad_points(planes_to_inliers.size() - 1, one_element_b_2);
				if (if_satisfied_b_2) {
					p_q.push(one_element_b_2);
					local_v_p_q.push_back(one_element_b_2);
				}
			}
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g_2;
			bool if_satisfied_g_2 = update_good_points(planes_to_inliers.size() - 1, one_element_g_2);
			if (if_satisfied_g_2) {
				p_q.push(one_element_g_2);
				local_v_p_q.push_back(one_element_g_2);
			}

			v_p_q.clear();
			v_p_q = local_v_p_q;



		}
		//conduct exclusion operation
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() != 0 && one_add.size() == 0 && regularization_planes.size() == 0) {
			for (int iii : one_remove) {
				points_if_added[iii]++;
			}
			std::vector<int> chaged_in_p;
			for (int id : planes_to_inliers[one_remove_plane]) {

				if (std::find(one_remove.begin(), one_remove.end(), id) == one_remove.end()) {
					chaged_in_p.push_back(id);
				}
			}
			

			planes_to_inliers[one_remove_plane] = chaged_in_p;




			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(chaged_in_p.size());

			for (int jj = 0; jj < chaged_in_p.size(); ++jj) {
				
				chaged_in_p_inliers.push_back(points[chaged_in_p[jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {
				
				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {

					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_remove_plane] = changed_plane;
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;


				
				planes_centroids_coplanar[one_remove_plane] = get_specific_centroids(one_remove_plane);
				planes_centroids[one_remove_plane] = get_specific_centroids(one_remove_plane);

				planes_if_regularized[one_remove_plane] = false;
				planes_to_coplanar_done[one_remove_plane] = -1;

				planes_to_parallel_done[one_remove_plane] = -1;

				planes_to_orthogonal_done[one_remove_plane] = -1;



				update_regular_relations_after_add_remove(one_remove_plane);
				
				update_regular_done_group_by_planes();
				

				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				
				p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();

				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {
						std::vector<int> this_pair = one_p_q.first.second;
					
						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;


							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}

						else {


							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);





						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;
						
						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_l2(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;
						
						if (if_merged[this_p_p[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}




						}


						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;


						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;

						if ( !if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_l2(this_p_p[0], one_element);

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}
				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		//conduct insertion operation		
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() != 0 && regularization_planes.size() == 0) {

			planes_to_inliers[one_add_plane].insert(planes_to_inliers[one_add_plane].end(), one_add.begin(), one_add.end());
			for (int iii : one_add) {
				points_if_added[iii]++;
			}

			

			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(planes_to_inliers[one_add_plane].size());

			for (int jj = 0; jj < planes_to_inliers[one_add_plane].size(); ++jj) {
				
				chaged_in_p_inliers.push_back(points[planes_to_inliers[one_add_plane][jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {

				
				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {


					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_add_plane] = changed_plane;


				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;


				planes_if_regularized[one_add_plane] = false;


				planes_centroids_coplanar[one_add_plane] = get_specific_centroids(one_add_plane);
				planes_centroids[one_add_plane] = get_specific_centroids(one_add_plane);


				planes_to_coplanar_done[one_add_plane] = -1;

				planes_to_parallel_done[one_add_plane] = -1;

				planes_to_orthogonal_done[one_add_plane] = -1;


				update_regular_relations_after_add_remove(one_add_plane);
				
				update_regular_done_group_by_planes();
			
				//update priority queue

				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				
				p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {


						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}

						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;



							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}
					}
					else if (one_p_q.first.first == 2) {


						std::vector<int> this_reg = one_p_q.first.second;
						
						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_l2(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {

						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}
						}
						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {

						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}


						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}

					}
					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (!if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_l2(this_p_p[0], one_element);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		//conduct regularity operation
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() != 0) {
		

			if_regularization_can_conducted[id_parallel_cluster] = false;



			//change planes function
		
			std::vector<std::pair<int, Inexact_Vector_3>> this_orthogonal_cluster = parallel_id_changed_normal_after_orthogonal[id_parallel_cluster];
			
			if (this_orthogonal_cluster.size() > 0) {
				//orthogonal
				int number_of_orth = orthogonal_done_to_planes.size();
				std::vector<int> this_or_p_ids;
				for (int i = 0; i < this_orthogonal_cluster.size(); ++i) {
					Inexact_Vector_3 this_normal = this_orthogonal_cluster[i].second;
					int parallel_cluster_id = this_orthogonal_cluster[i].first;
					std::vector<int> planes_ids = parallel_clusters_to_planes[parallel_cluster_id];
					if_regularization_can_conducted[parallel_cluster_id] = false;
					
					for (int id : planes_ids) {
					
						if_merged[id] = 0;

						planes_to_orthogonal_done[id] = number_of_orth;
						planes_to_parallel_done[id] = -1;
						this_or_p_ids.push_back(id);
						//orthogonal + coplanar
						planes_2[id] = Inexact_Plane(planes_centroids_coplanar[id], this_normal);
						planes_centroids[id] = planes_centroids_coplanar[id];
						planes_if_regularized[id] = true;
						
						
					}
				}
				orthogonal_done_to_planes.push_back(this_or_p_ids);

				//coplanar
				for (std::pair<int, Inexact_Vector_3> a_or : parallel_id_changed_normal_after_orthogonal[id_parallel_cluster]) {
					if (parallel_cluster_to_coplanr_cases[a_or.first].size() > 0) {
						for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[a_or.first]) {
							int id_copl_group = coplanar_done_to_planes.size();
							coplanar_done_to_planes.push_back(a_co);
							for (int id_c : a_co) {
								if_merged[id_c] = 0;

								planes_to_coplanar_done[id_c] = id_copl_group;


							}
						}
					}
				}

			}
			else {
			//parallel
				int number_of_parallel = parallel_done_to_planes.size();
				
				for (int id : parallel_clusters_to_planes[id_parallel_cluster]) {
					planes_to_parallel_done[id] = number_of_parallel;
					planes_to_orthogonal_done[id] = -1;

					if_merged[id] = 0;

					//parallel + coplanar.
					planes_2[id] = Inexact_Plane(planes_centroids_coplanar[id], normal_parallel_clusters[id_parallel_cluster]);
					planes_centroids[id] = planes_centroids_coplanar[id];
					planes_if_regularized[id] = true;
					


				}
				parallel_done_to_planes.push_back(parallel_clusters_to_planes[id_parallel_cluster]);

			//coplanar
				if (parallel_cluster_to_coplanr_cases[id_parallel_cluster].size() > 0) {
					for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[id_parallel_cluster]) {
						int id_copl_group = coplanar_done_to_planes.size();
						coplanar_done_to_planes.push_back(a_co);
						for (int id_c : a_co) {
							if_merged[id_c] = 0;

							planes_to_coplanar_done[id_c] = id_copl_group;


						}
					}
				}
				

			}

			for (int iii = 0; iii < parallel_id_changed_normal_after_orthogonal.size(); ++iii) {

				for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[iii]) {
					if (if_regularization_can_conducted[orth.first] == false) {
						if_regularization_can_conducted[iii] = false;
					}
				}


			}
			
			

		
			update_regular_done_group_by_planes();
			
			


			//update priority queue
			
			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			
			p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();

			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {


					std::vector<int> this_pair = one_p_q.first.second;


					if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						bool if_satisfied = convert_form_merge_l2(this_pair[0], this_pair[1], one_element);
						if (!if_satisfied) continue;
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);




					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;



						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);

					}
				}
				else if (one_p_q.first.first == 2) {


					std::vector<int> this_reg = one_p_q.first.second;
					
					if (if_merged[this_reg[0]] == 0)
					{
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
						bool if_satisfied = convert_form_split_l2(this_reg[0], max_list_vector, min_list_vector, one_element_1);
						if (!if_satisfied) continue;

						p_q.push(one_element_1);
						local_v_p_q.push_back(one_element_1);


					}

					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);


					}

				}
				else if (one_p_q.first.first == 3) {

					std::vector<int> this_p_p = one_p_q.first.second;

					if (if_merged[this_p_p[0]] == 0) {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
						bool if_satisfied_b = update_bad_points(this_p_p[0], one_element_b);
						if (if_satisfied_b) {
							p_q.push(one_element_b);
							local_v_p_q.push_back(one_element_b);
						}
					}
					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);


					}

				}
				else if (one_p_q.first.first == 4) {

					std::vector<int> this_p_p = one_p_q.first.second;

					if (if_merged[this_p_p[0]] == 0) {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
						bool if_satisfied_g = update_good_points(this_p_p[0], one_element_g);
						if (if_satisfied_g) {
							p_q.push(one_element_g);
							local_v_p_q.push_back(one_element_g);
						}

					}


					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);

					}

				}
				else if (one_p_q.first.first == 5) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if ( !if_regularization_can_conducted[this_p_p[0]]) {
						continue;
					}
					else {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						convert_form_regularization_l2(this_p_p[0], one_element);

						p_q.push(one_element);
						local_v_p_q.push_back(one_element);



					}

				}
			}


			v_p_q.clear();
			v_p_q = local_v_p_q;



		}
		else {
			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			bb = false;
		}
		if (t_m >= 1000) bb = false;

	} while (bb == true);

}


//****Fidelity metirc: L1.1.
void Shape_Detector::local_operators_normal() {

	bool bb = false;
	all_t_merge += t_merge;
	all_t_split += t_split;
	all_t_insert += t_insert;
	all_t_exlude += t_exlude;
	all_regularization += t_regularization;
	t_m = 0;
	t_merge = 0;
	t_split = 0;
	t_exlude = 0;
	t_insert = 0;
	t_regularization = 0;



	std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy> p_q;
	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> v_p_q;
	std::vector<std::vector<int>> max_list_vector;
	std::vector<std::vector<int>> min_list_vector;

	//add the merging operation to priority queue
	for (int i = 0; i < primitive_connection.size(); ++i) {


		std::vector<int> one = primitive_connection[i];
		for (int j = i + 1; j < one.size(); ++j) {
			if (one[j] >= 1) {

				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
				bool if_satisfied = convert_form_merge_normal(i, j, one_element);
				if (!if_satisfied)continue;

				p_q.push(one_element);
				v_p_q.push_back(one_element);

			}
		}
	}

	//add the spliting operation to priority queue
	for (int i = 0; i < planes_to_inliers.size(); ++i) {
		if (region_type[i] > 5) continue;
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_split_normal(i, max_list_vector, min_list_vector, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);


	}
	//add the exclusion operation to priority queue
	for (int i = 0; i < bad_points_shape.size(); ++i) {

		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

		bool if_satisfied = convert_form_exclude_normal(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}
	//add the insertion operation to priority queue
	for (int i = 0; i < good_points_shape.size(); ++i) {




		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_insert_normal(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}
	//add the regularity operation to priority queue

	for (int i = 0; i < energy_changing_initial_regularization.size(); ++i) {
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_regularization_normal(i, one_element);
		if (!if_satisfied) continue;

		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}





	do {


		bb = true;

		std::vector<int> if_merged;
		if_merged = std::vector<int>(planes_to_inliers.size(), -1);



		std::vector<int> one_merged;//ids of merged primitives.
		std::vector<int> id_of_moved_merge;//ids of the ecluded points because of merging.
		std::vector<int> max_list;
		std::vector<int> min_list;


		int one_splite = -1;//id of split primitive.
		std::vector<int> one_remove;//ids of excluded inliers.
		int one_remove_plane = -1;//id of the primitive, from which the inliers are excluded.
		std::vector<int> one_add;//ids of inserted outliers.
		int one_add_plane = -1;//id of the primitive, to which the outliers are inserted.
		std::vector<int> regularization_planes;
		int id_parallel_cluster;

		//test if the top operation is available, if it's valable, the operation is used to update the configuration and the priority queue is

		//also updated, otherwise, the top operation is poped and the next operation is tested until the top opeartion can not decrease the energy.

		while (!p_q.empty()) {


			if (p_q.top().second[0] > 0) break;//when the top operation can not decrease the energy, stop.
			if (p_q.top().first.first == 1) {//merge operation in the top of priority queue.

				std::vector<int> this_pair = p_q.top().first.second;//ids of the pair of merged primitives + ids of excluded points.

				int f = this_pair[0];
				int s = this_pair[1];
				if ((planes_to_inliers[f].size() + planes_to_inliers[s].size() + 2 - this_pair.size()) < min_points) {//when the new primitives have less point than threshold, the operation will be ignoed.
					p_q.pop();
					continue;
				}



				double after_operation_distance_visual = all_normal_diaviation + p_q.top().second[1];
				if (if_constraint) {

					if (number_of_assigned_points + 2 - this_pair.size() < ori_inliers_number) {//when the current inlier number is smaller than the original one.

						p_q.pop();
						continue;
					}

					if (after_operation_distance_visual / double(number_of_assigned_points) < ori_mean_normal_diviation) {//when the current normal deviation is bigger than orig

						p_q.pop();
						continue;


					}

				}

				p_q.pop();
				id_of_moved_merge.clear();
				//saving the points that far away from the new plane.
				if (this_pair.size() > 2) {
					for (int k = 2; k < this_pair.size(); ++k) {
						id_of_moved_merge.push_back(this_pair[k]);
					}
				}
				number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();

				if_merged[f] = 0;
				if_merged[s] = 0;
				one_merged.push_back(f);
				one_merged.push_back(s);

				all_normal_diaviation = after_operation_distance_visual;
				mean_normal_current = all_normal_diaviation / number_of_assigned_points;

				t_m++;
				t_merge++;


				break;


			}
			else if (p_q.top().first.first == 2) {//split operator in the top of priority queue.


				int splited_index = p_q.top().first.second[0];
				int info_index = p_q.top().first.second[1];

				double after_splited_distance = all_normal_diaviation + p_q.top().second[1];



				max_list = max_list_vector[info_index];
				min_list = min_list_vector[info_index];

				if (if_constraint) {
					//when current normal deviation is bigger than original one or current number of primitives is bigger than original one.
					if (after_splited_distance / double(number_of_assigned_points) < ori_mean_normal_diviation || (planes_2.size() + 1) > ori_primitives_number) {

						p_q.pop();
						continue;
					}


				}

				p_q.pop();
				if_merged[splited_index] = 0;

				one_splite = splited_index;

				all_normal_diaviation = after_splited_distance;
				mean_normal_current = all_normal_diaviation / number_of_assigned_points;
				t_m++;
				t_split++;


				break;

			}
			else if (p_q.top().first.first == 3) {//exclusion operator in the top of priority queue.

				int plane_id = p_q.top().first.second[0];

				double after_delete_distance = all_normal_diaviation + p_q.top().second[1];

				if (if_constraint) {
					//when the current number of inliers is smaller than original one or the current normal devidation is bigger than original one.
					if (number_of_assigned_points + 1 - p_q.top().first.second.size() < ori_inliers_number || after_delete_distance / double(number_of_assigned_points + 1 - p_q.top().first.second.size()) < ori_mean_normal_diviation) {

						p_q.pop();
						continue;
					}
				}
				//when the inliers number of the primitive is smaller than threhsold.
				if (planes_to_inliers[plane_id].size() + 1 - p_q.top().first.second.size() < min_points) {

					p_q.pop();
					continue;
				}


				if_merged[plane_id] = 0;
				one_remove_plane = plane_id;
				std::vector<int> points_ids = p_q.top().first.second;
				points_ids.erase(points_ids.begin());
				one_remove = points_ids;

				all_normal_diaviation = after_delete_distance;

				number_of_assigned_points = number_of_assigned_points - one_remove.size();
				mean_normal_current = all_normal_diaviation / number_of_assigned_points;

				p_q.pop();
				t_m++;
				t_exlude++;


				break;

			}
			else if (p_q.top().first.first == 4) {//insert operator in the top of priority queue.

				int plane_id = p_q.top().first.second[0];


				double after_add_distance = all_normal_diaviation + p_q.top().second[1];
				if (if_constraint) {
					//when the current normal deviation is bigger than original one.
					if (after_add_distance / double(number_of_assigned_points - 1 + p_q.top().first.second.size()) < ori_mean_normal_diviation) {
						p_q.pop();
						continue;
					}
				}
				if_merged[plane_id] = 0;
				one_add_plane = plane_id;
				std::vector<int> points_ids = p_q.top().first.second;
				points_ids.erase(points_ids.begin());
				one_add = points_ids;

				all_normal_diaviation = after_add_distance;

				number_of_assigned_points = number_of_assigned_points + one_add.size();
				mean_normal_current = all_normal_diaviation / number_of_assigned_points;

				p_q.pop();
				t_m++;
				t_insert++;



				break;


			}

			else if (p_q.top().first.first == 5) {

				//regularity operation
				if (lambda_regular == 0) {
					p_q.pop();

					continue;

				}
				id_parallel_cluster = p_q.top().first.second[0];
				std::vector<std::pair<int, Inexact_Vector_3>> this_orthogonal_cluster = parallel_id_changed_normal_after_orthogonal[id_parallel_cluster];
				bool if_regulared = false;
				for (int i = 0; i < this_orthogonal_cluster.size(); ++i) {

					int parallel_cluster_id = this_orthogonal_cluster[i].first;
					if (!if_regularization_can_conducted[parallel_cluster_id]) {
						if_regulared = true;
						break;


					}
				}

				if (if_regulared) {
					p_q.pop();
					continue;
				}


				double after_add_distance = all_normal_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (after_add_distance / double(number_of_assigned_points) <= ori_mean_normal_diviation) {
					
						p_q.pop();
						continue;
					}
				}




				regularization_planes = p_q.top().first.second;
				regularization_planes.erase(regularization_planes.begin());
				
				all_normal_diaviation = after_add_distance;
				mean_normal_current = all_normal_diaviation / number_of_assigned_points;


				p_q.pop();
				t_m++;
				t_regularization++;

				break;


			}
		}



		//implement the operation and update the priority queue.

		//conduct merging operation
		if (one_merged.size() != 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() == 0) {//merge operator.
			int last_size = planes_2.size();
			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			int id_1 = one_merged[0];
			int id_2 = one_merged[1];
			std::vector<int> one_merge_points_real;
			std::vector<int> one_merge_points = planes_to_inliers[id_1];
			one_merge_points.insert(one_merge_points.end(), planes_to_inliers[id_2].begin(), planes_to_inliers[id_2].end());
			for (int id : one_merge_points) {
				if (std::find(id_of_moved_merge.begin(), id_of_moved_merge.end(), id) == id_of_moved_merge.end()) {
					one_merge_points_real.push_back(id);
				}
			}

			one_merge_points = one_merge_points_real;
			for (int m_id : id_of_moved_merge) {//exclude the far away inliers.
				inliers_to_planes[m_id] = -1;

			}
			if (one_merge_points.size() >= min_points) {

				bool if_regular_done = false;
				/*if (planes_to_parallel_done[id_2] != -1 || planes_to_parallel_done[id_1] != -1 || planes_to_orthogonal_done[id_1] != -1 || planes_to_orthogonal_done[id_2] != -1 || planes_to_coplanar_done[id_1] != -1 || planes_to_coplanar_done[id_2] != -1) {
					if_regular_done = true;

				}*/
				if (id_2 < id_1) {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_2) {
							respective_planes[t] = t;
						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_1) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_1);
					region_type.erase(region_type.begin() + id_2);



					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);

					planes_2.erase(planes_2.begin() + id_1);
					planes_2.erase(planes_2.begin() + id_2);

					



					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_1);
					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_2);

					planes_centroids.erase(planes_centroids.begin() + id_1);
					planes_centroids.erase(planes_centroids.begin() + id_2);

					planes_if_regularized.erase(planes_if_regularized.begin() + id_1);
					planes_if_regularized.erase(planes_if_regularized.begin() + id_2);

					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_1);
					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_2);
					planes_to_coplanar_done.push_back(-1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_2);
					planes_to_parallel_done.push_back(-1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_2);
					planes_to_orthogonal_done.push_back(-1);
				}
				else {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_1) {
							respective_planes[t] = t;
						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_2) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_2);
					region_type.erase(region_type.begin() + id_1);

					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);

					planes_2.erase(planes_2.begin() + id_2);
					planes_2.erase(planes_2.begin() + id_1);


					

					

					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_2);

					planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + id_1);

					planes_centroids.erase(planes_centroids.begin() + id_2);

					planes_centroids.erase(planes_centroids.begin() + id_1);


					planes_if_regularized.erase(planes_if_regularized.begin() + id_2);


					planes_if_regularized.erase(planes_if_regularized.begin() + id_1);

					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_2);
					planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + id_1);
					planes_to_coplanar_done.push_back(-1);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_2);
					planes_to_parallel_done.erase(planes_to_parallel_done.begin() + id_1);
					planes_to_parallel_done.push_back(-1);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_2);
					planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + id_1);
					planes_to_orthogonal_done.push_back(-1);


				}


				planes_if_regularized.push_back(false);

				//update regularization cluster id
				update_regular_relations_after_merge(id_1, id_2, respective_planes);




				planes_to_inliers.push_back(one_merge_points);

				planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
				planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 1));

				
				


				std::vector<Inexact_Point_3> inliers;
				inliers.reserve(one_merge_points.size());

				for (int jj = 0; jj < one_merge_points.size(); ++jj) {

					inliers.push_back(points[one_merge_points[jj]].first);

				}
				Inexact_Plane plane;
				if (inliers.size() < 3) {

					continue;
				}
				linear_least_squares_fitting_3(inliers.begin(), inliers.end(), plane, CGAL::Dimension_tag<0>());

				planes_2.push_back(plane);
				

				//update inliers_to_planes.
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;
				//if (if_regular_done) {
					update_regular_done_group_by_planes();
				//}
				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();//clean p_q.
				//update the existing operations
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {//when the operator is merge.
						std::vector<int> this_pair = one_p_q.first.second;
						if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//when that is the merged pair.
						{


							continue;
						}
						else if ((if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) || (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0)) {//when only one primitive is merged.


							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(respective_planes[this_pair[1]], respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
						else {//when the operation is not affected.



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o = this_pair;
							o[0] = respective_planes[this_pair[0]];
							o[1] = respective_planes[this_pair[1]];

							one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
					}
					else if (one_p_q.first.first == 2) {//when the operator is split.
						std::vector<int> this_reg = one_p_q.first.second;
						if (if_merged[this_reg[0]] == 0)//when the primitive has been merged.
						{
							continue;
						}

						else {//when the operation is not affected.
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o.push_back(respective_planes[this_reg[0]]);
							o.push_back(this_reg[1]);
							one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}

					}
					else if (one_p_q.first.first == 3) {//when the operator is exclude.
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 4) {//when the operator is insert.
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];
							one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}

					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (!if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_normal(this_p_p[0], one_element);

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}
				}


				//add new operator, split the new merged primitive
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
				bool if_satisfied = convert_form_split_normal(last_size - 2, max_list_vector, min_list_vector, one_element_1);
				if (if_satisfied) {
					p_q.push(one_element_1);
					local_v_p_q.push_back(one_element_1);
				}
				//add new operat, exlude and inerte for the new merged primitive
				if (all_t_transfer > 0) {
					std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
					bool if_satisfied_b = update_bad_points_normal(last_size - 2, one_element_b);
					if (if_satisfied_b) {
						p_q.push(one_element_b);
						local_v_p_q.push_back(one_element_b);
					}
				}
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
				bool if_satisfied_g = update_good_points_normal(last_size - 2, one_element_g);
				if (if_satisfied_g) {
					p_q.push(one_element_g);
					local_v_p_q.push_back(one_element_g);
				}

				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		//conduct splitting operation
		else if (one_merged.size() == 0 && one_splite != -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() == 0) {//split operator.
			int last_size = planes_2.size();


			std::vector<int> respective_planes = std::vector<int>(last_size, -1);
			bool if_regular_done = false;
			/*if (planes_to_parallel_done[one_splite] != -1 || planes_to_orthogonal_done[one_splite] != -1 || planes_to_coplanar_done[one_splite] != -1) {
				if_regular_done = true;

			}*/

			for (int t = 0; t < last_size; ++t) {
				if (t < one_splite) {
					respective_planes[t] = t;
				}
				else if (t == one_splite) {
					respective_planes[t] = last_size;
				}

				else {
					respective_planes[t] = t - 1;

				}


			}
			planes_to_inliers.erase(planes_to_inliers.begin() + one_splite);
			region_type.push_back(region_type[one_splite]);
			region_type.push_back(region_type[one_splite]);
			region_type.erase(region_type.begin() + one_splite);



			planes_2.erase(planes_2.begin() + one_splite);

			planes_to_coplanar_done.erase(planes_to_coplanar_done.begin() + one_splite);
			planes_to_coplanar_done.push_back(-1);
			planes_to_coplanar_done.push_back(-1);
			planes_to_parallel_done.erase(planes_to_parallel_done.begin() + one_splite);
			planes_to_parallel_done.push_back(-1);
			planes_to_parallel_done.push_back(-1);
			planes_to_orthogonal_done.erase(planes_to_orthogonal_done.begin() + one_splite);
			planes_to_orthogonal_done.push_back(-1);
			planes_to_orthogonal_done.push_back(-1);

			planes_to_inliers.push_back(max_list);
			planes_to_inliers.push_back(min_list);


			planes_if_regularized.erase(planes_if_regularized.begin() + one_splite);
			planes_if_regularized.push_back(false);
			planes_if_regularized.push_back(false);

			

			planes_centroids_coplanar.erase(planes_centroids_coplanar.begin() + one_splite);
			planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
			planes_centroids_coplanar.push_back(get_specific_centroids(planes_to_inliers.size() - 2));

			planes_centroids.erase(planes_centroids.begin() + one_splite);
			planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 1));
			planes_centroids.push_back(get_specific_centroids(planes_to_inliers.size() - 2));

			update_regular_relations_after_splite(one_splite, respective_planes);


			std::vector<Inexact_Point_3> max_inliers;
			max_inliers.reserve(max_list.size());

			for (int jj = 0; jj < max_list.size(); ++jj) {

				max_inliers.push_back(points[max_list[jj]].first);

			}
			Inexact_Plane max_plane;
			if (max_inliers.size() < 3) {
				getchar();
			}
			linear_least_squares_fitting_3(max_inliers.begin(), max_inliers.end(), max_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(max_plane);
		

			std::vector<Inexact_Point_3> min_inliers;
			min_inliers.reserve(min_list.size());

			for (int jj = 0; jj < min_list.size(); ++jj) {

				min_inliers.push_back(points[min_list[jj]].first);

			}
			Inexact_Plane min_plane;
			if (min_inliers.size() < 5) {
				getchar();
			}
			linear_least_squares_fitting_3(min_inliers.begin(), min_inliers.end(), min_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(min_plane);


			//update inliers to plane
			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			//if (if_regular_done) {
				update_regular_done_group_by_planes();
			//}
			//update priority queue
			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			while (!p_q.empty()) p_q.pop();
			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {
					std::vector<int> this_pair = one_p_q.first.second;

					if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//not possible case
					{
						continue;
					}
					else if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 1, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 2, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}



					}
					else if (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 1, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 2, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}


					}
					else {// is not affected
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o = this_pair;
						o[0] = respective_planes[this_pair[0]];
						o[1] = respective_planes[this_pair[1]];


						one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);




					}
				}
				else if (one_p_q.first.first == 2) {
					std::vector<int> this_reg = one_p_q.first.second;

					if (if_merged[this_reg[0]] == 0)
					{
						continue;

					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o.push_back(respective_planes[this_reg[0]]);
						o.push_back(this_reg[1]);
						one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 3) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 4) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}

				else if (one_p_q.first.first == 5) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if ( !if_regularization_can_conducted[this_p_p[0]]) {
						continue;
					}
					else {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						convert_form_regularization_normal(this_p_p[0], one_element);

						p_q.push(one_element);
						local_v_p_q.push_back(one_element);



					}
				}
			}

			//add new operators, splite the new two primitives.
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
			bool if_satisfied = convert_form_split_normal(planes_to_inliers.size() - 1, max_list_vector, min_list_vector, one_element_1);
			if (if_satisfied) {
				p_q.push(one_element_1);
				local_v_p_q.push_back(one_element_1);
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_2;
			bool if_satisfied_2 = convert_form_split_normal(planes_to_inliers.size() - 2, max_list_vector, min_list_vector, one_element_2);
			if (if_satisfied_2) {
				p_q.push(one_element_2);
				local_v_p_q.push_back(one_element_2);
			}
			//add new operators, exlude and insert for two new primitives.
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
				bool if_satisfied_b = update_bad_points_normal(planes_to_inliers.size() - 2, one_element_b);
				if (if_satisfied_b) {
					p_q.push(one_element_b);
					local_v_p_q.push_back(one_element_b);
				}
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
			bool if_satisfied_g = update_good_points_normal(planes_to_inliers.size() - 2, one_element_g);
			if (if_satisfied_g) {
				p_q.push(one_element_g);
				local_v_p_q.push_back(one_element_g);
			}
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b_2;
				bool if_satisfied_b_2 = update_bad_points_normal(planes_to_inliers.size() - 1, one_element_b_2);
				if (if_satisfied_b_2) {
					p_q.push(one_element_b_2);
					local_v_p_q.push_back(one_element_b_2);
				}
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g_2;
			bool if_satisfied_g_2 = update_good_points_normal(planes_to_inliers.size() - 1, one_element_g_2);
			if (if_satisfied_g_2) {
				p_q.push(one_element_g_2);
				local_v_p_q.push_back(one_element_g_2);
			}

			v_p_q.clear();
			v_p_q = local_v_p_q;



		}
		//conduct exclusion operation
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() != 0 && one_add.size() == 0 && regularization_planes.size() == 0) {
			for (int iii : one_remove) {
				points_if_added[iii]++;
			}
			std::vector<int> chaged_in_p;
			for (int id : planes_to_inliers[one_remove_plane]) {

				if (std::find(one_remove.begin(), one_remove.end(), id) == one_remove.end()) {
					chaged_in_p.push_back(id);
				}
			}


			planes_to_inliers[one_remove_plane] = chaged_in_p;




			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(chaged_in_p.size());

			for (int jj = 0; jj < chaged_in_p.size(); ++jj) {

				chaged_in_p_inliers.push_back(points[chaged_in_p[jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {

				bool if_regular_done = false;
				/*if (planes_to_parallel_done[one_remove_plane] != -1 || planes_to_orthogonal_done[one_remove_plane] != -1 || planes_to_coplanar_done[one_remove_plane] != -1) {
					if_regular_done = true;

				}*/
				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {

					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_remove_plane] = changed_plane;
				//update inliers to plane
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;


			
				planes_centroids_coplanar[one_remove_plane] = get_specific_centroids(one_remove_plane);
				planes_centroids[one_remove_plane] = get_specific_centroids(one_remove_plane);

				planes_if_regularized[one_remove_plane] = false;
				planes_to_coplanar_done[one_remove_plane] = -1;

				planes_to_parallel_done[one_remove_plane] = -1;

				planes_to_orthogonal_done[one_remove_plane] = -1;

				
				update_regular_relations_after_add_remove(one_remove_plane);
				//if (if_regular_done) {
					update_regular_done_group_by_planes();
				//}


				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {
						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;

						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_normal(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points_normal(this_p_p[0], one_element_b);
							if (if_satisfied_b) {

								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;


						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points_normal(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if ( !if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_normal(this_p_p[0], one_element);

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;


			}
		}
		//conduct insertion operation
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() != 0 && regularization_planes.size() == 0) {

			planes_to_inliers[one_add_plane].insert(planes_to_inliers[one_add_plane].end(), one_add.begin(), one_add.end());
			for (int iii : one_add) {
				points_if_added[iii]++;
			}



			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(planes_to_inliers[one_add_plane].size());

			for (int jj = 0; jj < planes_to_inliers[one_add_plane].size(); ++jj) {

				chaged_in_p_inliers.push_back(points[planes_to_inliers[one_add_plane][jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {

				bool if_regular_done = false;
				/*if (planes_to_parallel_done[one_add_plane] != -1 || planes_to_orthogonal_done[one_add_plane] != -1 || planes_to_coplanar_done[one_add_plane] != -1) {
					if_regular_done = true;

				}*/
				Inexact_Plane changed_plane;

				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_add_plane] = changed_plane;




				//update inlier to plane
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;


				planes_if_regularized[one_add_plane] = false;
		
				planes_centroids_coplanar[one_add_plane] = get_specific_centroids(one_add_plane);
				planes_centroids[one_add_plane] = get_specific_centroids(one_add_plane);


				planes_to_coplanar_done[one_add_plane] = -1;

				planes_to_parallel_done[one_add_plane] = -1;

				planes_to_orthogonal_done[one_add_plane] = -1;

				//update regularization cluster id



				update_regular_relations_after_add_remove(one_add_plane);
				//if (if_regular_done) {
					update_regular_done_group_by_planes();
				//}
				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {


						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}
					}
					else if (one_p_q.first.first == 2) {


						std::vector<int> this_reg = one_p_q.first.second;

						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_normal(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {


						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points_normal(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}
						}
						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {


						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points_normal(this_p_p[0], one_element_g);
							if (if_satisfied_g) {

								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}

					}
					else if (one_p_q.first.first == 5) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (!if_regularization_can_conducted[this_p_p[0]]) {
							continue;
						}
						else {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							convert_form_regularization_normal(this_p_p[0], one_element);

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		//update regularity operation
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0 && regularization_planes.size() != 0) {



			if_regularization_can_conducted[id_parallel_cluster] = false;
			//change planes function
			std::vector<std::pair<int, Inexact_Vector_3>> this_orthogonal_cluster = parallel_id_changed_normal_after_orthogonal[id_parallel_cluster];


			if (this_orthogonal_cluster.size() > 0) {
				//orthogonal
				int number_of_orth = orthogonal_done_to_planes.size();
				std::vector<int> this_or_p_ids;
				for (int i = 0; i < this_orthogonal_cluster.size(); ++i) {
					Inexact_Vector_3 this_normal = this_orthogonal_cluster[i].second;
					int parallel_cluster_id = this_orthogonal_cluster[i].first;
					std::vector<int> planes_ids = parallel_clusters_to_planes[parallel_cluster_id];
					if_regularization_can_conducted[parallel_cluster_id] = false;
					for (int id : planes_ids) {
						if_merged[id] = 0;

						planes_to_orthogonal_done[id] = number_of_orth;
						planes_to_parallel_done[id] = -1;
						this_or_p_ids.push_back(id);
						
						planes_2[id] = Inexact_Plane(planes_centroids_coplanar[id], this_normal);
						planes_centroids[id] = planes_centroids_coplanar[id];
						planes_if_regularized[id] = true;

					}
				}
				orthogonal_done_to_planes.push_back(this_or_p_ids);
				//coplanar
				for (std::pair<int, Inexact_Vector_3> a_or : this_orthogonal_cluster) {
					if (parallel_cluster_to_coplanr_cases[a_or.first].size() > 0) {
						for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[a_or.first]) {
							int id_copl_group = coplanar_done_to_planes.size();
							coplanar_done_to_planes.push_back(a_co);
							for (int id_c : a_co) {
								if_merged[id_c] = 0;

								planes_to_coplanar_done[id_c] = id_copl_group;


							}
						}
					}
				}

			}
			else {
				//parallel
				int number_of_parallel = parallel_done_to_planes.size();
				for (int id : parallel_clusters_to_planes[id_parallel_cluster]) {
					planes_to_parallel_done[id] = number_of_parallel;
					planes_to_orthogonal_done[id] = -1;
					if_merged[id] = 0;

					
					planes_2[id] = Inexact_Plane(planes_centroids_coplanar[id], normal_parallel_clusters[id_parallel_cluster]);
					planes_centroids[id] = planes_centroids_coplanar[id];
					planes_if_regularized[id] = true;

				}
				parallel_done_to_planes.push_back(parallel_clusters_to_planes[id_parallel_cluster]);

				//coplanar
				if (parallel_cluster_to_coplanr_cases[id_parallel_cluster].size() > 0) {
					for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[id_parallel_cluster]) {
						int id_copl_group = coplanar_done_to_planes.size();
						coplanar_done_to_planes.push_back(a_co);
						for (int id_c : a_co) {
							if_merged[id_c] = 0;

							planes_to_coplanar_done[id_c] = id_copl_group;
							



						}
					}
				}

			}


			
			for (int iii = 0; iii < parallel_id_changed_normal_after_orthogonal.size(); ++iii) {

				for (std::pair<int, Inexact_Vector_3> orth : parallel_id_changed_normal_after_orthogonal[iii]) {
					if (if_regularization_can_conducted[orth.first] == false) {
						if_regularization_can_conducted[iii] = false;
					}
				}


			}

			

			update_regular_done_group_by_planes();


			//update priority queue
			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			
			p_q = std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy>();

			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {


					std::vector<int> this_pair = one_p_q.first.second;


					if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						bool if_satisfied = convert_form_merge_normal(this_pair[0], this_pair[1], one_element);
						if (!if_satisfied) continue;
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);




					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;



						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);

					}
				}
				else if (one_p_q.first.first == 2) {


					std::vector<int> this_reg = one_p_q.first.second;
					
					if (if_merged[this_reg[0]] == 0)
					{
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
						bool if_satisfied = convert_form_split_normal(this_reg[0], max_list_vector, min_list_vector, one_element_1);
						if (!if_satisfied) continue;

						p_q.push(one_element_1);
						local_v_p_q.push_back(one_element_1);


					}

					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);


					}

				}
				else if (one_p_q.first.first == 3) {

					std::vector<int> this_p_p = one_p_q.first.second;

					if (if_merged[this_p_p[0]] == 0) {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
						bool if_satisfied_b = update_bad_points_normal(this_p_p[0], one_element_b);
						if (if_satisfied_b) {
							p_q.push(one_element_b);
							local_v_p_q.push_back(one_element_b);
						}
					}
					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);


					}

				}
				else if (one_p_q.first.first == 4) {

					std::vector<int> this_p_p = one_p_q.first.second;

					if (if_merged[this_p_p[0]] == 0) {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
						bool if_satisfied_g = update_good_points_normal(this_p_p[0], one_element_g);
						if (if_satisfied_g) {
							p_q.push(one_element_g);
							local_v_p_q.push_back(one_element_g);
						}

					}


					else {

						p_q.push(one_p_q);
						local_v_p_q.push_back(one_p_q);

					}

				}
				else if (one_p_q.first.first == 5) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if ( !if_regularization_can_conducted[this_p_p[0]]) {
						continue;
					}
					else {

						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
						convert_form_regularization_normal(this_p_p[0], one_element);

						p_q.push(one_element);
						local_v_p_q.push_back(one_element);



					}

				}
			}


			v_p_q.clear();
			v_p_q = local_v_p_q;



		}


		else {


			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			bb = false;
		}
		if (t_m >= 1000) bb = false;

	} while (bb == true);

}

//**************

//**************

//**************


void Shape_Detector::get_last_information() {
	last_coverage = coverage;
	last_mean_error = mean_error;
	last_primitives_number = primitives_number;
	last_normal_deviation = mean_normal_diviation;
	last_freedom_planes = freedom_of_planes;


}
//find all the adjancent pairs of primitives.
void Shape_Detector::test_connected_primitives() {
	std::vector<int> one_connect;
	one_connect = std::vector<int>(planes_to_inliers.size(), 0);
	primitive_connection = std::vector<std::vector<int>>(planes_to_inliers.size(), one_connect);
	int Nb_neigh = 10;
	for (size_t i = 0; i < points.size(); i++) {
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;




		for (int it = 0; it < Nb_neigh; it++) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {
				primitive_connection[inliers_to_planes[neighbor_index]][inliers_to_planes[i]] ++;
				primitive_connection[inliers_to_planes[i]][inliers_to_planes[neighbor_index]] ++;

			}


		}

	}

	for (int i = 0; i < primitive_connection.size(); ++i) {
		std::vector<int> this_raw = primitive_connection[i];
		for (int j = i; j < this_raw.size(); ++j) {
			if (this_raw[j] >= 10) {
				primitive_connection[i][j] = 1;
				primitive_connection[j][i] = 1;

			}
			else {
				primitive_connection[i][j] = 0;
				primitive_connection[j][i] = 0;


			}

		}

	}
}
void Shape_Detector::get_distance_diviation() {
	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;
	freedom_of_planes = 0;
	freedom_of_planes_current = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		
		const Inexact_Plane & H = planes_2[i];


		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());
		}
	}
	freedom_of_planes = get_the_recent_degrees_of_freedom();

	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;
	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_distance_current = mean_distance_diaviation;
	mean_normal_current = mean_normal_diviation;
	freedom_of_planes_current = freedom_of_planes;

}
void Shape_Detector::get_good_points_normal() {
	if_insert_candidate = std::vector<int>(points.size(), -1);
	good_points_shape.clear();
	int Nb_neigh = 10;
	std::map<int, std::vector<std::pair<int, double>>> plane_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;


		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);

			}


		}

		if (one_neight_id.empty()) { continue; }

		double mm_normal = normal_threshold;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > mm_normal) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < epsilon) {

					mm_normal = abs(this_p.second * planes_2[neight_id].orthogonal_vector());
					changed_plane_id = neight_id;

				}

			}
		}


		if (changed_plane_id != -1) {
			plane_outliers[changed_plane_id].push_back(std::make_pair(i, mm_normal));//plane_id, point_id, normal_deviation.

		}

	}
	std::map<int, std::vector<std::pair<int, double>>>::iterator g_i;
	g_i = plane_outliers.begin();
	while (g_i != plane_outliers.end()) {
		std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator_normal> p_q;

		for (std::pair<int, double> ppp : g_i->second) {
			p_q.push(ppp);
		}

		std::vector<int> added_points_ids;
		double en = 0;
		double dif_e = 0;
		double n_a = 0;

		while (!p_q.empty() && n_a < number_of_insert_exclude) {
			added_points_ids.push_back(p_q.top().first);
			p_q.pop();
			n_a++;
		}
		dif_e = add_changed_error_normal(g_i->first, added_points_ids);
		en = energy_changed_second_normal(dif_e, -n_a, g_i->first);

		if (added_points_ids.size() != 0 && en < 0) {
			for (int idd : added_points_ids) {
				if_insert_candidate[idd] = 0;
			}
			std::vector<double> dd;
			dd.push_back(en);
			dd.push_back(dif_e);
			good_points_shape.push_back(std::make_pair(std::make_pair(g_i->first, added_points_ids), dd));//plane_id,point_ids,energy_change,distance
		}
		g_i++;
	}

}





double Shape_Detector::add_changed_error_normal(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	for (int u : befor_ids) {

		dif_befor += abs(befor_p.orthogonal_vector()*points[u].second);

	}
	befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());
	std::vector<Inexact_Point_3> after_points;

	for (int n : befor_ids) {
		after_points.push_back(points[n].first);

	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (int u : befor_ids) {
		dif_after += abs(after_p.orthogonal_vector()*points[u].second);
	}

	return (dif_after - dif_befor);


}

void Shape_Detector::get_bad_points_normal() {

	bad_points_shape.clear();
	if (old_coverage > ori_coverage) {
		std::map<int, std::vector<std::pair<int, double>>> plane_inliers;
		double bigest_normal = mean_normal_diviation;
		//double bigest_normal = 0.01;

		for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

			const Inexact_Plane & H = planes_2[i];


			for (int j : planes_to_inliers[i]) {
				const Point_with_normal & pt = points[j];
				if (points_if_added[j] > 2) continue;

				if (abs(H.orthogonal_vector()*pt.second) < bigest_normal) {



					plane_inliers[i].push_back(std::make_pair(j, abs(H.orthogonal_vector()*pt.second)));
				}

			}
		}

		std::map<int, std::vector<std::pair<int, double>>>::iterator b_o;
		b_o = plane_inliers.begin();

		while (b_o != plane_inliers.end()) {


			std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator_normal> p_q;
			for (std::pair<int, double> ppp : b_o->second) {
				p_q.push(ppp);

			}
			std::vector<int> removed_points_ids;
			double en = 0;
			double dif_e = 0;

			double n_b = 0;
			while (!p_q.empty() && n_b < number_of_insert_exclude) {
				removed_points_ids.push_back(p_q.top().first);
				p_q.pop();
				n_b++;
			}
			dif_e = remove_changed_error_normal(b_o->first, removed_points_ids);
			en = energy_changed_second_normal(dif_e, n_b, b_o->first);

			if (removed_points_ids.size() != 0 && (planes_to_inliers[b_o->first].size() - removed_points_ids.size()) >= min_points && en < 0) {
				std::vector<double> dd;
				dd.push_back(en);
				dd.push_back(dif_e);
				bad_points_shape.push_back(std::make_pair(std::make_pair(b_o->first, removed_points_ids), dd));
			}
			b_o++;
		}

	}
}
double Shape_Detector::remove_changed_error_normal(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	std::vector<int> after_ids;
	for (int u : befor_ids) {
		//dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
		dif_befor += abs(points[u].second * befor_p.orthogonal_vector());

		if (std::find(point_ids.begin(), point_ids.end(), u) == point_ids.end()) {
			after_ids.push_back(u);

		}
	}




	std::vector<Inexact_Point_3> after_points;
	std::vector<int> after_points_ids;
	for (int n : after_ids) {
		after_points.push_back(points[n].first);
		after_points_ids.push_back(n);
	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (int u : after_points_ids) {
		//dif_after += sqrt(CGAL::squared_distance(after_p, u));
		dif_after += abs(points[u].second * after_p.orthogonal_vector());
	}

	return (dif_after - dif_befor);


}


bool Shape_Detector::separate_two_out_normal(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif) {
	std::vector<int> this_region = planes_to_inliers[id];
	Inexact_Plane this_plane = planes_2[id];

	double max_value = FLT_MIN;
	int max_index = -1;
	double min_value = FLT_MIN;
	int min_index = -1;
	double dis_bef = 0;
	double mean_distance = 0;
	//divide into two regions

	Inexact_Vector_3 plane_direction = this_plane.orthogonal_vector();

	for (int ii = 0; ii < this_region.size(); ++ii) {

		double this_dis = sqrt(CGAL::squared_distance(this_plane, points[this_region[ii]].first));
		double this_normal_d = std::abs(points[this_region[ii]].second*plane_direction);
		mean_distance += this_dis;
		dis_bef += this_normal_d;
		if (this_plane.has_on_positive_side(points[this_region[ii]].first)) {
			if (this_dis > max_value) {
				max_value = this_dis;
				max_index = this_region[ii];

			}
		}
		else {
			if (this_dis > min_value) {
				min_value = this_dis;
				min_index = this_region[ii];
			}

		}
	}
	mean_distance /= double(this_region.size());
	//we do not split the primitive that has smaller distance error than the ori mean.
	if (mean_distance < ori_mean_error) {
		min_list.clear();
		max_list.clear();
		dif = 0;
		return false;
	}
	std::map<int, int> label_points;
	Inexact_Point_3 max_point = points[max_index].first;
	Inexact_Point_3 min_point = points[min_index].first;
	std::vector<Inexact_Point_3> max_point_list;
	std::vector<Inexact_Point_3> min_point_list;
	Inexact_Point_2 max_point_2d = this_plane.to_2d(max_point);
	Inexact_Point_2 min_point_2d = this_plane.to_2d(min_point);


	for (int j = 0; j < this_region.size(); ++j) {
		if ((max_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length() < (min_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length()) {
			max_list.push_back(this_region[j]);
			label_points[this_region[j]] = 1;
			max_point_list.push_back(points[this_region[j]].first);
		}
		else {
			label_points[this_region[j]] = -1;
			min_point_list.push_back(points[this_region[j]].first);
			min_list.push_back(this_region[j]);
		}
	}

	//transformation
	Inexact_Plane plane_max;
	Inexact_Plane plane_min;
	if (max_point_list.size() < 3 || min_point_list.size() < 3) {

		return false;
	}
	linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
	linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
	bool propo = false;
	int Nb_neigh = 10;
	int refine_time = 0;
	do {
		refine_time++;
		int moving_n = 0;
		propo = true;
		std::vector<int> if_moved = std::vector<int>(this_region.size(), -1);
		for (int f = 0; f < this_region.size(); ++f) {
			if (if_moved[f] > 2) continue;
			int this_label = label_points[this_region[f]];
			if ((int)spherical_neighborhood[this_region[f]].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[this_region[f]].size();
			}
			bool bb = false;
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[this_region[f]][it];
				if (inliers_to_planes[neighbor_index] != id) continue; //the neighbor point should be in the splitted region
				if (label_points[neighbor_index] != this_label) {
					bb = true;
					break;
				}


			}
			if (bb == false) continue;
			Inexact_Point_3 this_point = points[this_region[f]].first;



			if (this_label == -1) {
				if (CGAL::squared_distance(plane_max, this_point) < epsilon) {
					if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) > abs(points[this_region[f]].second * plane_min.orthogonal_vector())) {
						//if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f]++;
					}
				}
			}
			else {
				if (CGAL::squared_distance(plane_min, this_point) < epsilon) {


					if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) < abs(points[this_region[f]].second * plane_min.orthogonal_vector())) {//need to changed for smooth																																														 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f]++;
					}
				}
			}







		}
		std::map<int, int>::iterator iter;
		iter = label_points.begin();
		max_list.clear();
		min_list.clear();
		max_point_list.clear();
		min_point_list.clear();
		while (iter != label_points.end()) {
			if (iter->second == 1) {
				max_list.push_back(iter->first);
				max_point_list.push_back(points[iter->first].first);

			}
			else {
				min_list.push_back(iter->first);
				min_point_list.push_back(points[iter->first].first);
			}
			iter++;
		}
		if (max_point_list.size() < min_points || min_point_list.size() < min_points) {

			return false;
		}
		linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
		linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
		if (moving_n < 5 || refine_time>25) {
			propo = false;
		}

	} while (propo);

	double dis_after = 0;

	for (int i = 0; i < max_list.size(); ++i) {
		dis_after += abs(points[max_list[i]].second * plane_max.orthogonal_vector());

	}

	for (int i = 0; i < min_list.size(); ++i) {
		dis_after += abs(points[min_list[i]].second * plane_min.orthogonal_vector());
	}

	dif = dis_after - dis_bef;

	return true;


}


bool Shape_Detector::merge_distance_changed_normal_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids) {
	Inexact_Plane p_i = planes_2[i];
	Inexact_Plane p_j = planes_2[j];
	double dis_before = 0;

	std::vector<Inexact_Point_3> assigned_pts;
	std::vector<int> assigned_pts_id;

	assigned_pts.reserve(planes_to_inliers[i].size() + planes_to_inliers[j].size());
	for (int pt_index : planes_to_inliers[i]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);
		dis_before += abs(points[pt_index].second * p_i.orthogonal_vector());
	}

	for (int pt_index : planes_to_inliers[j]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);

		dis_before += abs(points[pt_index].second * p_j.orthogonal_vector());


	}





	size_t n_poly = assigned_pts.size();
	//dis_before /= double(n_poly);


	double dP_norm = 0;

	Inexact_Plane test_plane;
	linear_least_squares_fitting_3(assigned_pts.begin(), assigned_pts.end(), test_plane, CGAL::Dimension_tag<0>());
	/*for (int id_in : assigned_pts_id) {

		dP_norm += abs(points[id_in].second * test_plane.orthogonal_vector());

		move_ids.clear();
	}*/


	Inexact_Plane new_plane;
	std::vector<Inexact_Point_3> new_assigned_pts;
	std::vector<int> new_assigned_pts_ids;
	for (int pt : assigned_pts_id) {
		if (sqrt(CGAL::squared_distance(test_plane, points[pt].first)) > epsilon) {
			move_ids.push_back(pt);
		}
		else {
			new_assigned_pts.push_back(points[pt].first);
			new_assigned_pts_ids.push_back(pt);

		}
	}
	if (new_assigned_pts_ids.size() < min_points) {
		dif = -FLT_MAX;

		return false;
	}
	linear_least_squares_fitting_3(new_assigned_pts.begin(), new_assigned_pts.end(), new_plane, CGAL::Dimension_tag<0>());
	for (int p : new_assigned_pts_ids) {

		dP_norm += abs(points[p].second * new_plane.orthogonal_vector());

	}



	dif = dP_norm - dis_before;

	return true;




}


bool Shape_Detector::convert_form_merge_normal(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(i);
	o.push_back(j);
	double dif_m = 0;
	std::vector<double> s_d;
	std::vector<int> merge_moved_ids;
	//we suppose that the pair of primitives can not merged,when their normal cosine smaller than normal_threshold.
	if (abs(planes_2[i].orthogonal_vector()*planes_2[j].orthogonal_vector()) < normal_threshold) {
		return false;
	}
	bool if_good = merge_distance_changed_normal_with_epsilon(i, j, dif_m, merge_moved_ids);
	if (!if_good) { return false; }
	/*if (merge_moved_ids.size() > points.size()*0.01|| merge_moved_ids.size() > planes_to_inliers[i].size()/4|| merge_moved_ids.size() > planes_to_inliers[j].size()/4) {
		return false;
	}*/
	if ((planes_to_inliers[i].size() + planes_to_inliers[j].size() - merge_moved_ids.size()) < min_points) {
		return false;
	}
	double energy_ccc = energy_changed_normal_merge(dif_m, -1, merge_moved_ids.size(),i,j);
	if (energy_ccc > 0) {
		return false;
	}
	s_d.push_back(energy_ccc);
	o.insert(o.end(), merge_moved_ids.begin(), merge_moved_ids.end());

	s_d.push_back(dif_m);
	one_element = std::make_pair(std::make_pair(1, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_split_normal(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element) {
	std::vector<int> max_list;
	std::vector<int> min_list;
	double dif_s = 0;
	bool if_split = separate_two_out_normal(i, max_list, min_list, dif_s);

	double energy_ccc = energy_changed_normal(dif_s, 1.0,i);

	if (energy_ccc >= 0 || max_list.size() < min_points || min_list.size() < min_points || !if_split) {
		return false;
	}
	max_list_vector.push_back(max_list);
	min_list_vector.push_back(min_list);
	std::vector<int> o;
	o.push_back(i);
	o.push_back(max_list_vector.size() - 1);
	std::vector<double> s_d;
	s_d.push_back(energy_ccc);

	s_d.push_back(dif_s);
	one_element = std::make_pair(std::make_pair(2, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_exclude_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(bad_points_shape[i].first.first);

	o.insert(o.end(), bad_points_shape[i].first.second.begin(), bad_points_shape[i].first.second.end());

	one_element = std::make_pair(std::make_pair(3, o), bad_points_shape[i].second);
	return true;
}

bool Shape_Detector::convert_form_insert_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;
	o.push_back(good_points_shape[i].first.first);
	o.insert(o.end(), good_points_shape[i].first.second.begin(), good_points_shape[i].first.second.end());


	one_element = std::make_pair(std::make_pair(4, o), good_points_shape[i].second);
	return true;
}

bool Shape_Detector::update_good_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	int Nb_neigh = 10;
	std::vector<std::pair<int, double>> ids_good_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (if_insert_candidate[i] == 0) continue;
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;//avoide infinite loop

		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;

		bool if_correspond = false;
		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);
				if (inliers_to_planes[neighbor_index] == id_shape) {

					if_correspond = true;
				}
			}


		}

		if (one_neight_id.empty() || !if_correspond) { continue; }



		double mm_normal = normal_threshold;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > mm_normal) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < epsilon) {

					mm_normal = abs(this_p.second * planes_2[neight_id].orthogonal_vector());
					changed_plane_id = neight_id;

				}

			}
		}


		if (changed_plane_id == id_shape) {

			ids_good_outliers.push_back(std::make_pair(i, mm_normal));//plane_id, point_id, normal_deviation.

		}

	}


	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator_normal> p_q;

	for (std::pair<int, double> ppp : ids_good_outliers) {
		p_q.push(ppp);
	}

	std::vector<int> added_points_ids;
	double en = 0;
	double dif_e = 0;
	double n_a = 0;

	while (!p_q.empty() && n_a < number_of_insert_exclude) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		added_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_a++;
	}
	dif_e = add_changed_error_normal(id_shape, added_points_ids);
	en = energy_changed_second_normal(dif_e, -n_a, id_shape);

	if (added_points_ids.size() != 0 && en < 0) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);

		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), added_points_ids.begin(), added_points_ids.end());


		one_element = std::make_pair(std::make_pair(4, o), dd);


		return true;

	}
	else {
		return false;
	}


}


bool Shape_Detector::update_bad_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {



	std::vector<std::pair<int, double>> ids_bad_inliers;
	double bigest_normal = mean_normal_diviation;
	//double bigest_normal = 0.01;



	const Inexact_Plane & H = planes_2[id_shape];


	for (int j : planes_to_inliers[id_shape]) {
		const Point_with_normal & pt = points[j];
		if (points_if_added[j] > 2) continue;

		if (abs(H.orthogonal_vector()*pt.second) < bigest_normal) {
			ids_bad_inliers.push_back(std::make_pair(j, abs(H.orthogonal_vector()*pt.second)));
		}

	}





	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator_normal> p_q;
	for (std::pair<int, double> ppp : ids_bad_inliers) {
		p_q.push(ppp);

	}
	std::vector<int> removed_points_ids;
	double en = 0;
	double dif_e = 0;

	double n_b = 0;
	while (!p_q.empty() && n_b < number_of_insert_exclude) {
		removed_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_b++;
	}
	dif_e = remove_changed_error_normal(id_shape, removed_points_ids);
	en = energy_changed_second_normal(dif_e, n_b, id_shape);

	if (removed_points_ids.size() != 0 && (planes_to_inliers[id_shape].size() - removed_points_ids.size()) >= min_points && en < 0) {
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);
		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), removed_points_ids.begin(), removed_points_ids.end());


		one_element = std::make_pair(std::make_pair(3, o), dd);
		return true;
	}
	else { return false; }



}

bool Shape_Detector::test_if_connected(int i, int j) {



	int Nb_neigh = 10;
	for (int f = 0; f < planes_to_inliers[i].size(); ++f) {

		if ((int)spherical_neighborhood[planes_to_inliers[i][f]].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[planes_to_inliers[i][f]].size();
		}

		for (int it = 0; it < Nb_neigh; it++) {

			int neighbor_index = spherical_neighborhood[planes_to_inliers[i][f]][it];

			for (int id : planes_to_inliers[j]) {
				if (id == neighbor_index) {

					return true;
				}

			}

		}

	}
	return false;
}



void Shape_Detector::get_distance_diviation_show_merge_info(double t) {


	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;
	freedom_of_planes = 0;
	freedom_of_planes_current = 0;

	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];
		


		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}
	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;

	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_normal_current = mean_normal_diviation;

	mean_distance_current = mean_distance_diaviation;
	freedom_of_planes = get_the_recent_degrees_of_freedom();
	freedom_of_planes_current = freedom_of_planes;
	double new_coverage = double(number_of_assigned_points) / double(points.size());

    _logger->debug("Local operators: {} times" ,t);
    _logger->debug("Merge operators: {} times",t_merge);
    _logger->debug("Split operators: {} times" , t_split );
    _logger->debug("Insert operators: {} times",t_insert );
    _logger->debug("Exclude operators: {} times",t_exlude );
    _logger->debug("Regularity operators: {} times", t_regularization);

    _logger->debug("Primitives : {}" ,size_current_primitives);

    _logger->debug("Coverage  : {}" , new_coverage);
    _logger->debug("Mean error : {}" , (mean_distance_diaviation) );
    _logger->debug("Mean normal deviation : {}" , (mean_normal_diviation) );
    _logger->debug("Degree of freedom : {}" , (freedom_of_planes));

    _logger->debug("Primitives reducing : {}" , (old_size_current_primitives - size_current_primitives));

    _logger->debug("Coverage adding   : {}" , new_coverage - old_coverage );
    _logger->debug("Mean error reducing : {}" , (old_mean_distance_diaviation - mean_distance_diaviation) );
    _logger->debug("Mean normal deviation adding : {}" , (mean_normal_diviation - old_mean_normal_diaviation) );
    _logger->debug("Degree of freedom changing : {}", (freedom_of_planes-old_freedom_of_planes));

	old_size_current_primitives = size_current_primitives;
	old_mean_distance_diaviation = mean_distance_diaviation;
	old_coverage = new_coverage;
	old_mean_normal_diaviation = mean_normal_diviation;
	old_freedom_of_planes = freedom_of_planes;


}


void Shape_Detector::get_distance_diviation_show_normal_info(double t) {


	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;
	freedom_of_planes = 0;
	freedom_of_planes_current = 0;




	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];
		

		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}

	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;


	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_normal_current = mean_normal_diviation;

	mean_distance_current = mean_distance_diaviation;
	freedom_of_planes = get_the_recent_degrees_of_freedom();
	freedom_of_planes_current = freedom_of_planes;

	double new_coverage = double(number_of_assigned_points) / double(points.size());





    _logger->debug("Transfer operator: {} s." , t );
    _logger->debug("Primitives : {}" , (size_current_primitives));

    _logger->debug("Coverage  : {}" , (new_coverage) );
    _logger->debug("Mean error : {}" , (mean_distance_diaviation) );
    _logger->debug("Mean normal deviation : {}" , (mean_normal_diviation));
    _logger->debug("Degree of freedom : {}" , (freedom_of_planes));
    _logger->debug("Primitives reducing : {}" , (old_size_current_primitives - size_current_primitives));

    _logger->debug("Coverage adding   : {}" , new_coverage - old_coverage );
    _logger->debug("Mean error reducing : {}" ,(old_mean_distance_diaviation - mean_distance_diaviation));
    _logger->debug("Mean normal deviation adding : {}" , (mean_normal_diviation - old_mean_normal_diaviation) );
    _logger->debug("Degree of freedom changing : {}" , (freedom_of_planes - old_freedom_of_planes));
	old_size_current_primitives = size_current_primitives;
	old_mean_distance_diaviation = mean_distance_diaviation;
	old_coverage = new_coverage;
	old_mean_normal_diaviation = mean_normal_diviation;
	old_freedom_of_planes = freedom_of_planes;




}


void Shape_Detector::get_coverage_and_mean_error_pure()
{
	number_of_assigned_points = 0;

	mean_error = 0;
	mean_normal_diviation = 0;
	freedom_of_planes = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];
		
		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];
			mean_error += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(H.orthogonal_vector()*pt.second);
		}
	}
	freedom_of_planes = get_the_recent_degrees_of_freedom();

	mean_error /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	coverage = double(number_of_assigned_points) / double(points.size());
	primitives_number = planes_to_inliers.size();

}


void Shape_Detector::do_region_growing()
{
	int class_index = -1;
//	int ten_percent_of_points = points.size() / 10;

    _logger->debug("Region growing...");

	for (size_t i = 0; i < points.size(); i++) {
//		if (i % ten_percent_of_points == 0) std::cout << ". " << std::flush;

		if (inliers_to_planes[i] == -1) {

			// Updates the index of primitive
			inliers_to_planes[i] = ++class_index;
			int conti = 0; 	//for accelerate least_square fitting 

			// Characteristics of the seed
			Inexact_Point_3 pt_seed = points[i].first;
			Inexact_Vector_3 normal_seed = points[i].second;
			Inexact_Plane optimal_plane(pt_seed, normal_seed);

			// Initialization containers
			std::vector<int> index_container; index_container.push_back(i);
			std::vector<int> index_container_former_ring; index_container_former_ring.push_back(i);
			std::list<int> index_container_current_ring;

			// Propagation
			bool propagation = true;
			do {
				propagation = false;

				for (int k = 0; k < (int)index_container_former_ring.size(); k++) {

					int point_index = index_container_former_ring[k];


					int Nb_neigh = (int)spherical_neighborhood[point_index].size();

					if ((int)spherical_neighborhood[point_index].size() > knn) {
						Nb_neigh = knn;
					}
					for (int it = 0; it < Nb_neigh; it++) {

						int neighbor_index = spherical_neighborhood[point_index][it];

						if (inliers_to_planes[neighbor_index] == -1) {

							Inexact_Point_3 neighbor = points[neighbor_index].first;
							Inexact_Point_3 neighbor_projection = optimal_plane.projection(neighbor);
							double distance = sqrt((neighbor - neighbor_projection).squared_length());

							if (distance < epsilon && abs(points[neighbor_index].second * optimal_plane.orthogonal_vector()) > normal_threshold) {

								inliers_to_planes[neighbor_index] = class_index;
								propagation = true;
								index_container_current_ring.push_back(neighbor_index);
								conti++;

								if ((conti < 50 && conti % 10 == 0) || (conti > 50 && conti % 500 == 0)) {
									std::list<Inexact_Point_3> listp;
									for (int pm = 0; pm < (int)index_container.size(); pm++) {
										Inexact_Point_3 ptza = points[index_container[pm]].first;
										listp.push_back(ptza);
									}

									if (listp.size() >= 3) {
										Inexact_Plane reajusted_plane;
										linear_least_squares_fitting_3(listp.begin(), listp.end(), reajusted_plane, CGAL::Dimension_tag<0>());
										optimal_plane = reajusted_plane;
									}
								}
							}
						}
					}
				}

				// Updates containers
				index_container_former_ring.clear();
				for (std::list<int>::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it) {
					index_container_former_ring.push_back(*it);
					index_container.push_back(*it);
				}
				index_container_current_ring.clear();

			} while (propagation);

			// Tests the number of inliers
			if (index_container.size() >= min_points && inliers_arent_aligned(index_container)) {
				planes_to_inliers.push_back(index_container);
			}
			else {
				--class_index;
				inliers_to_planes[i] = -1;
				for (int k = 0; k < (int)index_container.size(); k++) inliers_to_planes[index_container[k]] = -1;
			}
		}
	}
}



void Shape_Detector::export_region_growing_results(const std::string & filename)
{
	std::ofstream stream(filename);

	if (stream.is_open()) {
		stream << planes_to_inliers.size() << std::endl;

		for (size_t i = 0 ; i < planes_to_inliers.size() ; ++i) {

			const Inexact_Plane & H = planes_0[i];
			stream << H.a() << " " << H.b() << " " << H.c() << " " << H.d() << std::endl;

			size_t n = planes_to_inliers[i].size();
			stream << n << std::endl;

			for (size_t j = 0 ; j < n - 1 ; ++j) {
				stream << planes_to_inliers[i][j] << " ";
			}
			stream << planes_to_inliers[i][n - 1] << std::endl;
		}
		stream.close();
	}
}



void Shape_Detector::set_path_clusters(const std::string & filename)
{
	path_clusters = filename;
}



void Shape_Detector::load_region_growing_results()
{
	planes_to_inliers.clear();
	inliers_to_planes = std::vector<int>(points.size(), -1);

	FILE* stream = fopen(path_clusters.c_str(), "r");

	if (stream != NULL) {

		int size_planes;
		fscanf(stream, "%i\n", &size_planes);
		planes_to_inliers.reserve(size_planes);

		for (int i = 0 ; i < size_planes ; ++i) {
			double a, b, c, d;
			fscanf(stream, "%lf %lf %lf %lf\n", &a, &b, &c, &d);

			int curr_plane_size;
			fscanf(stream, "%i\n", &curr_plane_size);

			std::vector<int> P;
			P.reserve(curr_plane_size);
			int index;
			for (int j = 0 ; j < curr_plane_size - 1 ; ++j) {
				fscanf(stream, "%i ", &index);
				P.push_back(index);
			}

			fscanf(stream, "%i\n", &index);
			P.push_back(index);

			for (int j = 0 ; j < curr_plane_size ; ++j) inliers_to_planes[P[j]] = i;
			planes_to_inliers.push_back(P);
		}

		fclose(stream);
	}
}



void Shape_Detector::detect_planes(bool read_clusters)
{
	clock_t t_detect_start = clock();

	// Part 1.
	// Performs a region-growing algorithm.

	// Initializes structures

	planes_0.clear();
	planes_1.clear();
	planes_2.clear();
	non_coplanar_planes = 0;
	planes_centroids.clear();
	planes_to_colors.clear();

	if (path_point_cloud_extension == ".ply") {
		planes_to_inliers.clear();
		inliers_to_planes = std::vector<int>(points.size(), -1);

		if (read_clusters) {
			load_region_growing_results();
		} else {
			if (should_compute_knn) compute_average_spacing_and_k_nearest_neighbors();
			do_region_growing();
		}
	}

	// Part 2.
	// Extracts planes

	
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(100, 225);

	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

		//Inexact_Point_3 centroid = CGAL::ORIGIN;
		double xc = 0, yc = 0, zc = 0;

		std::vector<Inexact_Point_3> inliers_i;
		inliers_i.reserve(planes_to_inliers[i].size());

		for (size_t j = 0; j < planes_to_inliers[i].size(); ++j) {
			const Inexact_Point_3 & pt = points[planes_to_inliers[i][j]].first;
			inliers_i.push_back(pt);
			//centroid = CGAL::barycenter(centroid, j, pt, 1);
			xc += pt.x(), yc += pt.y(), zc += pt.z();
		}
		
		xc /= planes_to_inliers[i].size();
		yc /= planes_to_inliers[i].size();
		zc /= planes_to_inliers[i].size();
		Inexact_Point_3 centroid (xc, yc, zc);

		Inexact_Plane plane;
		linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());

		planes_0.push_back(plane);
		planes_1.push_back(plane);

		unsigned char r = 0, g = 0, b = 0;
        r = uniform_distribution(generator);
        g = uniform_distribution(generator);
        b = uniform_distribution(generator);
		planes_to_colors.push_back(CGAL::Color(r, g, b));

		planes_centroids.push_back(centroid);
	}

	// export_region_growing_results();

	clock_t t_detect_end = clock();
    _logger->debug("Plane detection : done in {} s.",(double(t_detect_end - t_detect_start) / CLOCKS_PER_SEC));
}



bool Shape_Detector::inliers_arent_aligned(const std::vector<int> & inds)
{
	size_t n = inds.size();

	for (size_t j = 0; j < n / 3; j++) {
		size_t id_0 = j, id_1 = n / 3 + j, id_2 = 2 * n / 3 + j;
		const Inexact_Point_3 & P_0 = points[id_0].first, &P_1 = points[id_1].first, &P_2 = points[id_2].first;

		Inexact_Vector_3 N = CGAL::cross_product(P_1 - P_0, P_2 - P_0);
		if (N != CGAL::NULL_VECTOR) {
			return true;
		}
	}

	return false;
}


void Shape_Detector::regularize_planes()
{
	clock_t t_regularize_start = clock();
	planes_1 = planes_0;

    CGAL::Shape_regularization::Planes::regularize_planes(points,
			Point_map(),
			planes_1,
			CGAL::Identity_property_map<Inexact_Plane>(),
			Shape_Detector_Index_Map(inliers_to_planes),
			true, true, true, false,
			tolerance_angle,
			tolerance_coplanarity,
			Inexact_Vector_3(0, 0, 1));

	clock_t t_regularize_end = clock();
    _logger->debug("Plane regularization done in {} s.",double(t_regularize_end - t_regularize_start) / CLOCKS_PER_SEC);
}


void Shape_Detector::get_coverage_and_mean_error()
{
	number_of_assigned_points = 0;

	mean_error = 0;
	mean_normal_diviation = 0;
	freedom_of_planes = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		
		const Inexact_Plane & H = planes_2[i];
		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];
			mean_error += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}
	all_error = mean_error;
	mean_error /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);

	freedom_of_planes = get_the_recent_degrees_of_freedom();
	coverage = double(number_of_assigned_points) / double(points.size());
	primitives_number = planes_to_inliers.size();
    _logger->debug("Primitives : {}",primitives_number);
    _logger->debug("Coverage   : {}",coverage);
    _logger->debug("Mean error : {}",mean_error);
    _logger->debug("Mean normal deviation : {}",mean_normal_diviation);
    _logger->debug("Freedom of planes : {}",freedom_of_planes);

}


void Shape_Detector::coset_color() {
	

	for (int i = 0; i < orthogonal_done_to_planes.size(); ++i) {
		
		for (int id : orthogonal_done_to_planes[i]) {
			planes_to_colors[id] = planes_to_colors[orthogonal_done_to_planes[i][0]];
			

		}
	}

	for (int i = 0; i < parallel_done_to_planes.size(); ++i) {
		for (int id : parallel_done_to_planes[i]) {
			if (planes_to_orthogonal_done[id] == -1) {
				planes_to_colors[id] = planes_to_colors[parallel_done_to_planes[i][0]];
			}

		}
	}

	
}
//test tool used to est the regularized planes.
void Shape_Detector::test_regularity_planes() {
	
	for (int i = 0; i < orthogonal_done_to_planes.size(); ++i) {
		std::vector< Inexact_Vector_3> cross_normals;
		for (int id : orthogonal_done_to_planes[i]) {
		
		
			Inexact_Vector_3 last_normal = planes_2[id].orthogonal_vector();
		

			for (int kk : orthogonal_done_to_planes[i]){
				Inexact_Vector_3 that_normal = planes_2[kk].orthogonal_vector();

				double ttt = last_normal * that_normal / (sqrt(last_normal.squared_length())*sqrt(that_normal.squared_length()));
				if (abs(ttt) < 0.99&&abs(ttt)>0.001) {
					cross_normals.push_back(CGAL::cross_product(last_normal, that_normal));
				}
		}
		}
		if (cross_normals.size() > 0) {
			Inexact_Vector_3 last_normal1 = cross_normals[cross_normals.size() - 1];
			for (int kk = 0; kk < cross_normals.size(); ++kk) {
				double ttt = last_normal1 * cross_normals[kk] / (sqrt(last_normal1.squared_length())*sqrt(cross_normals[kk].squared_length()));
				if (abs(ttt) < 0.999) {
					std::cout << ttt << "problem!!!!!!!!!!!!!!!" << std::endl;
				}
				
			}
		}

	

		}
	

	double nnnumber = 0;
	double all_did = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		nnnumber += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];



		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			all_did += sqrt(CGAL::squared_distance(H, pt.first));


		}
	}

	std::cout << "mean distance::::::::::" << all_did / nnnumber << std::endl;

}

void Shape_Detector::set_primitives()
{
	//cylinder_colors.clear();
	//cylinder_alpha_shape.clear();

	alpha_shapes_pts.clear();
	best_rectangles_pts.clear();
	convex_hulls_pts.clear();
	convex_hulls_seeds.clear();
	alpha_shapes_colors.clear();
	convex_hulls_colors.clear();
	best_rectangles_colors.clear();
	clock_t t_start = clock();
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(100, 225);
	unsigned char r = 0, g = 0, b = 0;
	//change the regularizarion plane colors
	coset_color();
	
	//std::vector< Inexact_Vector_3> regular_nromals;
	for (size_t u = 0 ; u < planes_2.size() ; ++u) {
		//if (planes_to_orthogonal_done[u] != 0) continue;
		
		Inexact_Plane H = planes_2[u];

		//regular_nromals.push_back(H.orthogonal_vector());

        r = uniform_distribution(generator);
        g = uniform_distribution(generator);
        b = uniform_distribution(generator);
	
		const CGAL::Color & col = (CGAL::Color(r, g, b));
	
		
		std::vector<Inexact_Point_3> assigned_pts;
		assigned_pts.reserve(planes_to_inliers[u].size());
		for (int pt_index : planes_to_inliers[u]) {
			assigned_pts.push_back(points[pt_index].first);
		}

		// Computes alpha-shapes.

		Inexact_Vector_3 n = H.orthogonal_vector(), b1 = H.base1(), b2 = H.base2();
		n = n * 1.0 / (sqrt(n*n));
		b1 = b1 * 1.0 / (sqrt(b1*b1));
		b2 = b2 * 1.0 / (sqrt(b2*b2));

		std::vector<Inexact_Point_2> p2d;
		std::map<Inexact_Point_2, int> vMap, vMap2;
		std::vector<Inexact_Point_3> p3d;

		auto it = assigned_pts.begin();
		Inexact_Point_3 O = H.projection(*it);
		p2d.push_back(Inexact_Point_2(0, 0));
		vMap[p2d.back()] = p2d.size() - 1;
		it++;
		while (it != assigned_pts.end()) {
			Inexact_Vector_3 p = (*it) - O;
			p2d.push_back(Inexact_Point_2(p * b1, p * b2));
			vMap[p2d.back()] = p2d.size() - 1;
			it++;
		}

        Alpha_Shape as(p2d.begin(), p2d.end());
//        as.set_alpha(double(0.005 * bbox_diagonal));
        as.set_alpha(double(bbox_diagonal));

        int ind_min = alpha_shapes_pts.size();

        std::list<Inexact_Point_3> las3d;

        CDT2 cdt;
        std::map<Inexact_Point_2, CDT2::Vertex_handle> pVMap;

        auto vit = as.alpha_shape_vertices_begin();
        while (vit != as.alpha_shape_vertices_end()) {
            CDT2::Vertex_handle v = cdt.insert((*vit)->point());
            pVMap[(*vit)->point()] = v;
            vit++;
        }

        auto eit = as.alpha_shape_edges_begin();
        while (eit != as.alpha_shape_edges_end()) {
            switch (as.classify(*eit)) {
            case Alpha_Shape::SINGULAR:
                eit++;
                continue;
            default:
                break;
            }

            CDT2::Vertex_handle v1 = pVMap[as.segment(*eit).source()];
            CDT2::Vertex_handle v2 = pVMap[as.segment(*eit).target()];
            if (!v1->is_valid() || !v2->is_valid())
                std::cout << "invalid!" << std::endl;
            cdt.insert_constraint(v1, v2);
            eit++;
			
        }

        int outside = 0;
        double area = 0;
        auto fit = cdt.finite_faces_begin();
        while (fit != cdt.finite_faces_end()) {
            int res;
            Inexact_Point_2 p1 = fit->vertex(0)->point();
            Inexact_Point_2 p2 = fit->vertex(1)->point();
            Inexact_Point_2 p3 = fit->vertex(2)->point();
            if (CGAL::collinear(p1, p2, p3)) {
                std::cout << "collinear" << std::endl;
                continue;
            }
            Inexact_Point_2 center = Inexact_Point_2(0, 0) + (((p1 - Inexact_Point_2(0, 0)) + (p2 - Inexact_Point_2(0, 0)) + (p3 - Inexact_Point_2(0, 0))) * 1.0 / 3.0);
            res = as.classify(center);
            if (res == Alpha_Shape::INTERIOR) {
                area += abs(CGAL::area(p1, p2, p3));
                for (int i = 0; i < 3; i++) {
                    Inexact_Point_2 p = fit->vertex(i)->point();

                    // if point is used in cdt, copy it to
                    if (vMap2.count(p) == 0) {
                        p3d.push_back(H.projection(assigned_pts[vMap[p]]));
                        vMap2[p] = vMap[p];
                        vMap[p] = p3d.size() - 1;
                    }
                }

            } else outside++;
            fit++;
        }

        fit = cdt.finite_faces_begin();
        while (fit != cdt.finite_faces_end()) {
            int res;

            Inexact_Point_2 p1 = fit->vertex(0)->point();
            Inexact_Point_2 p2 = fit->vertex(1)->point();
            Inexact_Point_2 p3 = fit->vertex(2)->point();
            Inexact_Point_2 center = Inexact_Point_2(0, 0) + (((p1 - Inexact_Point_2(0, 0)) + (p2 - Inexact_Point_2(0, 0)) + (p3 - Inexact_Point_2(0, 0))) * 1.0 / 3.0);
            res = as.classify(center);
            if (res == Alpha_Shape::INTERIOR) {
                if (vMap.count(p1) == 0) {
                    continue;
                }
                if (vMap.count(p2) == 0) {
                    continue;
                }
                if (vMap.count(p3) == 0) {
                    continue;
                }

                std::vector<Inexact_Point_3> tr(3);
                tr[0] = p3d[vMap[p1]];
                tr[1] = p3d[vMap[p2]];
                tr[2] = p3d[vMap[p3]];
                alpha_shapes_pts.push_back(tr);
                alpha_shapes_colors.push_back(col);

                las3d.push_back(tr[0]);
                las3d.push_back(tr[1]);
                las3d.push_back(tr[2]);
            }
            fit++;
        }

        int ind_max = alpha_shapes_pts.size() - 1;

        // Computes convex hulls.

        std::vector<Inexact_Point_3> as3d(las3d.begin(), las3d.end());
        std::vector<Inexact_Point_2> as2d(as3d.size());
        for (size_t j = 0; j < as3d.size(); ++j) {
            as2d[j] = H.to_2d(as3d[j]);
        }

	

		std::vector<Inexact_Point_2> ch2d;
        CGAL::ch_graham_andrew(as2d.begin(), as2d.end(), std::back_inserter(ch2d));

		std::vector<Inexact_Point_3> ch3d;
		ch3d.reserve(ch2d.size());
		for (size_t j = 0; j < ch2d.size(); ++j) {
			ch3d.push_back(H.to_3d(ch2d[j]));
		}

		convex_hulls_seeds.push_back(std::make_pair(ind_min, ind_max));

		convex_hulls_pts.push_back(ch3d);
		convex_hulls_colors.push_back(col);
	
		// Computes best fit rectangles.
		std::vector<Inexact_Point_3> R;

		double xm = 0, ym = 0;
		for (Inexact_Point_2 pt : ch2d) {
			xm += pt.x(), ym += pt.y();
		}
		xm /= ch2d.size();
		ym /= ch2d.size();
		Inexact_Point_2 K(xm, ym);

		std::vector<Inexact_Point_2> centered_ch2d;
		for (Inexact_Point_2 pt : ch2d) {
			double x = pt.x() - xm, y = pt.y() - ym;
			centered_ch2d.push_back(Inexact_Point_2(x, y));
		}

		double min_product_dims = FLT_MAX;
		int argmin_theta_d = -1;

		for (int theta_d = 0; theta_d <= 85; theta_d += 5) {

			double x_min = FLT_MAX, x_max = -FLT_MAX;
			double y_min = FLT_MAX, y_max = -FLT_MAX;
            double cos_theta = cos(theta_d * M_PI / 180);
            double sin_theta = sin(theta_d * M_PI / 180);

			for (size_t i = 0; i < centered_ch2d.size(); ++i) {
				double x_i = cos_theta * centered_ch2d[i].x() - sin_theta * centered_ch2d[i].y();
				double y_i = sin_theta * centered_ch2d[i].x() + cos_theta * centered_ch2d[i].y();
				if (x_i < x_min) x_min = x_i;
				if (x_i > x_max) x_max = x_i;
				if (y_i < y_min) y_min = y_i;
				if (y_i > y_max) y_max = y_i;
			}

			double product_dims = (x_max - x_min) * (y_max - y_min);

			if (product_dims < min_product_dims) {
				min_product_dims = product_dims;
				argmin_theta_d = theta_d;

				Inexact_Vector_2 inv_u(cos_theta, -sin_theta);
				Inexact_Vector_2 inv_v(sin_theta, cos_theta);

				R.clear();
				R.push_back(H.to_3d(K + x_min * inv_u + y_min * inv_v));
				R.push_back(H.to_3d(K + x_min * inv_u + y_max * inv_v));
				R.push_back(H.to_3d(K + x_max * inv_u + y_max * inv_v));
				R.push_back(H.to_3d(K + x_max * inv_u + y_min * inv_v));
			}
		}
		if (R.size() == 0) {

			
			best_rectangles_pts.push_back(ch3d);
			best_rectangles_colors.push_back(col);
			continue;


		}
		best_rectangles_pts.push_back(R);
		best_rectangles_colors.push_back(col);
	}

	clock_t t_end = clock();
    _logger->debug( "Compute primitives (alpha shape, best rectangle and convex hull) : {} s." ,double(t_end - t_start) / CLOCKS_PER_SEC );
}

//void Shape_Detector::set_default_parameters(){

//    set_discretization_parameters(0.5, 0.4);
//    set_constraint(false);
//    set_weight_m(0);
//    set_metric_type(0);
//    set_max_iter(-1);

//    set_lambda_fidelity(1.0);
//    set_lambda_c(1.0);
//    set_lambda_r(1.0);
//    set_lambda_regularity(1.0);
//}

void Shape_Detector::set_stop_iteration(bool cc) {
	if_stop = cc;
}

void Shape_Detector::set_constraint(bool cc) {
	if_constraint = cc;
}


void Shape_Detector::set_lambda_r(double db) {
	lambda_r = db;
    _logger->debug("lambda_r {}",lambda_r );

}


void Shape_Detector::set_lambda_regularity(double db) {
	lambda_regular = db;
    _logger->debug("lambda_regular {}",lambda_regular );

}
void Shape_Detector::set_lambda_c(double db) {
	lambda_c = db;
    _logger->debug("lambda_c {}" ,lambda_c);

}
void Shape_Detector::set_lambda_fidelity(double db) {
	lambda_fidelity = db;
    _logger->debug("lambda_fidelity {}", lambda_fidelity );

}
void Shape_Detector::set_weight_m(int wm) {
	weight_mode = wm;
}
void Shape_Detector::set_max_iter(int mi) {
    max_iter = mi;
}

void Shape_Detector::set_metric_type(int m) {
	metric_type = m;
}

void Shape_Detector::clear_primitives()
{
	alpha_shapes_pts.clear();
	alpha_shapes_colors.clear();

	best_rectangles_pts.clear();
	best_rectangles_colors.clear();

	convex_hulls_pts.clear();
	convex_hulls_colors.clear();
	convex_hulls_seeds.clear();
}



std::tuple<Inexact_Point_3, Inexact_Vector_3, CGAL::Color, CGAL::Color> Shape_Detector::get_vertex(int point_index)
{
	const Point_with_normal & pt = points[point_index];
	const Inexact_Point_3 & position = pt.first;
	const Inexact_Vector_3 & normal = pt.second;

	CGAL::Color primitive_color(0, 0, 0);
	CGAL::Color natural_color = inliers_to_natural_colors[point_index];

	if (!inliers_to_planes.empty() && !planes_to_colors.empty()) {
		int shape_index = inliers_to_planes[point_index];
		if (shape_index != -1) {
			primitive_color = planes_to_colors[shape_index];
		}
	}
	return std::make_tuple(position, normal, primitive_color, natural_color);
}



void Shape_Detector::get_inliers(size_t primitive_index, std::vector<Point_with_normal> & inliers)
{
	const std::vector<int> & indices = planes_to_inliers[primitive_index];

	inliers.reserve(indices.size());
	for (int id : indices) inliers.push_back(points[id]);
}


void Shape_Detector::get_inliers(size_t primitive_index, std::vector<size_t> & inliers)
{
	const std::vector<int> & indices = planes_to_inliers[primitive_index];

	inliers.reserve(indices.size());
	for (int id : indices) inliers.push_back(id);
}


size_t Shape_Detector::get_total_number_of_inliers()
{
	size_t n = 0;

	for (size_t i = 0 ; i < planes_to_inliers.size() ; ++i) {
		n += planes_to_inliers[i].size();
	}

	return n;
}


int Shape_Detector::get_cloud_size() const
{
	return int(points.size());
}


double Shape_Detector::get_x_min() const
{
	return x_min;
}


double Shape_Detector::get_x_max() const
{
	return x_max;
}


double Shape_Detector::get_y_min() const
{
	return y_min;
}


double Shape_Detector::get_y_max() const
{
	return y_max;
}


double Shape_Detector::get_z_min() const
{
	return z_min;
}


double Shape_Detector::get_z_max() const
{
	return z_max;
}


double Shape_Detector::get_bbox_diagonal() const
{
	return bbox_diagonal;
}


double Shape_Detector::get_average_spacing() const
{
	return average_spacing;
}


int Shape_Detector::get_min_points() const
{
	return min_points;
}


double Shape_Detector::get_epsilon() const
{
	return epsilon;
}


int Shape_Detector::get_knn() const
{
	return knn;
}


double Shape_Detector::get_normal_threshold() const
{
	return normal_threshold;
}


int Shape_Detector::get_neighborhood_algorithm() const
{
	return neighborhood_algorithm;
}


double Shape_Detector::get_neighborhood_distance() const
{
	return neighborhood_distance;
}


double Shape_Detector::get_regularization_lambda() const
{
	return lambda;
}


double Shape_Detector::get_tolerance_angle() const
{
	return tolerance_angle;
}


double Shape_Detector::get_tolerance_coplanarity() const
{
	return tolerance_coplanarity;
}




void Shape_Detector::to_vg(const std::string & directory, std::string & basename, int X, int Y, int Z)
{

    Pwn_vector T(points_begin(),points_end());

	// Step 1.
	// For each point, we find in which subdivision it is located.

	double x_inf = FLT_MAX, x_sup = -FLT_MAX;
	double y_inf = FLT_MAX, y_sup = -FLT_MAX;
	double z_inf = FLT_MAX, z_sup = -FLT_MAX;

	for (size_t i = 0 ; i < T.size() ; ++i) {
		const Inexact_Point_3 & pt = T[i].first;
		const double &x = pt.x(), &y = pt.y(), &z = pt.z();

		if (x < x_inf) x_inf = x;
		if (x > x_sup) x_sup = x;
		if (y < y_inf) y_inf = y;
		if (y > y_sup) y_sup = y;
		if (z < z_inf) z_inf = z;
		if (z > z_sup) z_sup = z;
	}

	double dx = x_sup - x_inf;
	double dy = y_sup - y_inf;
	double dz = z_sup - z_inf;

	typedef std::tuple<int, int, int> Triplet;
	std::vector<Triplet> vertices_locations;
	vertices_locations.reserve(T.size());

	for (size_t i = 0 ; i < T.size() ; ++i) {
		const Inexact_Point_3 & pt = T[i].first;

		int a = jclamp(0, int(X * (pt.x() - x_inf) / dx), X - 1);
		int b = jclamp(0, int(Y * (pt.y() - y_inf) / dy), Y - 1);
		int c = jclamp(0, int(Z * (pt.z() - z_inf) / dz), Z - 1);
		vertices_locations.push_back(Triplet(a, b, c));
	}

	// Step 2.
	// For each planar shape we need to know in which subdivisions it can be found

	std::vector< std::set<Triplet> > shapes_locations;
	shapes_locations.reserve(planes_to_inliers.size());

	for (size_t i = 0 ; i < planes_to_inliers.size() ; ++i) {
		std::set<Triplet> L;
		for (size_t j = 0 ; j < planes_to_inliers[i].size(); ++j) {
			int ind_ij = planes_to_inliers[i][j];
			L.insert(vertices_locations[ind_ij]);
		}
		shapes_locations.push_back(L);
	}

	// Step 3.
	// Prints a .vg file for each subdivision

	for (int i = 0 ; i < X ; ++i) {
		for (int j = 0 ; j < Y ; ++j) {
			for (int k = 0 ; k < Z ; ++k) {
				Triplet tr (i, j, k);
				
				// Step 3.1
				// Collects points in T_curr

				std::vector<int> local_indices(T.size(), -1);

				Pwn_vector T_curr;
				for (size_t p = 0 ; p < T.size(); ++p) {
					if (vertices_locations[p] == tr) {
						local_indices[p] = T_curr.size();
						T_curr.push_back(T[p]);
					}
				}

				// Step 3.2
				// Gets planes associated to the current triplet

				std::vector<Inexact_Plane> S_eq_curr;
				std::vector<std::vector<int> > S_inl_curr;

				for (size_t p = 0 ; p < planes_to_inliers.size() ; ++p) {
					if (shapes_locations[p].find(tr) == shapes_locations[p].end()) continue;

					S_eq_curr.push_back(planes_2[p]);
					
					std::vector<int> I_p;
					for (size_t q = 0 ; q < planes_to_inliers[p].size() ; ++q) {
						int ind = planes_to_inliers[p][q];
						if (local_indices[ind] != -1) {
							I_p.push_back(local_indices[ind]);
						}
					}
					S_inl_curr.push_back(I_p);
				}

                _logger->debug("[{},{},{}] : {}",i , j , k , S_eq_curr.size() );

				// Step 3.3
				// Prints a vg file

                std::string filename = directory + "/" + basename + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + ".vg";
				std::ofstream stream(filename, std::ios::out);

				// Step 3.3.1
				// Provides a description of the point cloud

				size_t N = T_curr.size();
				stream << "num_points: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					const Inexact_Point_3 & pt = T_curr[i].first;
					stream << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
				}
				stream << "num_colors: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					stream << "128 128 128" << std::endl;
				}
				stream << "num_normals: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					const Inexact_Vector_3 & n = T_curr[i].second;
					stream << n.x() << " " << n.y() << " " << n.z() << std::endl;
				}

				// Step 3.3.2
				// Provides a description of the detected planes.

				size_t M = S_eq_curr.size();
				stream << "num_groups: " << M << std::endl;

				for (size_t i = 0; i < S_eq_curr.size(); ++i) {
					stream << "group_type: 0" << std::endl;
					stream << "num_group_parameters: 4" << std::endl;

					const Inexact_Plane & H = S_eq_curr[i];
					stream << "group_parameters: " << H.a() << " " << H.b() << " " << H.c() << " " << H.d() << std::endl;

					CGAL::Color C = planes_to_colors[i];
					const std::vector<int> & shape_assigned_pts = S_inl_curr[i];

                    stream << "group_label: " << i << std::endl;
					stream << "group_color: " << int(C.red()) << " " << int(C.green()) << " " << int(C.blue()) << std::endl;
					stream << "group_num_points: " << shape_assigned_pts.size() << std::endl;

					for (size_t j = 0; j < shape_assigned_pts.size(); ++j) {
						stream << shape_assigned_pts[j] << (j == shape_assigned_pts.size() - 1 ? '\n' : ' ');
					}
					stream << "num_children: 0" << std::endl;
				}

				stream.close();
			}
		}
	}
}

#include <xtensor-io/xnpz.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
void Shape_Detector::to_npz(const string& filename)
{

    Pwn_vector T(points_begin(),points_end());

    vector<double> points;
    vector<double> normals;

    for (size_t i = 0; i < T.size(); ++i){

        points.push_back(T[i].first.x());
        points.push_back(T[i].first.y());
        points.push_back(T[i].first.z());

        normals.push_back(T[i].second.x());
        normals.push_back(T[i].second.y());
        normals.push_back(T[i].second.z());

    }

    vector<double> group_parameters;
    vector<int> group_colors;
    vector<int> group_num_points;
    vector<int> group_points;

    // Part 2.
    // Provides a description of the detected planes.
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(100, 255);

    for (size_t i = 0; i < planes_2.size(); ++i) {

        group_parameters.push_back(planes_2[i].a());
        group_parameters.push_back(planes_2[i].b());
        group_parameters.push_back(planes_2[i].c());
        group_parameters.push_back(planes_2[i].d());

        const std::vector<int> & shape_assigned_pts = planes_to_inliers[i];

        group_colors.push_back(distribution(generator));
        group_colors.push_back(distribution(generator));
        group_colors.push_back(distribution(generator));

        group_num_points.push_back(shape_assigned_pts.size());

        for (size_t j = 0 ; j < shape_assigned_pts.size() ; ++j) {
            group_points.push_back(shape_assigned_pts[j]);
        }
    }

    std::vector<std::size_t> shape;

    bool compress = true;
    bool append = false;
    shape = { T.size(), 3 };
    auto xpoints = xt::adapt(points, shape);
    xt::dump_npz(filename,"points",xpoints,compress,append);
    append = true;
    auto xnormals = xt::adapt(normals, shape);
    xt::dump_npz(filename,"normals",xnormals,compress,append);

    shape = { planes_2.size(), 4 };
    auto xgroup_parameters = xt::adapt(group_parameters, shape);
    xt::dump_npz(filename,"group_parameters",xgroup_parameters,compress,append);

    shape = { planes_2.size(), 3 };
    auto xgroup_colors = xt::adapt(group_colors, shape);
    xt::dump_npz(filename,"group_colors",xgroup_colors,compress,append);

    shape = { planes_2.size(), 1 };
    auto xgroup_num_points = xt::adapt(group_num_points, shape);
    xt::dump_npz(filename,"group_num_points",xgroup_num_points,compress,append);

    shape = { group_points.size(), 1 };
    auto xgroup_points = xt::adapt(group_points, shape);
    xt::dump_npz(filename,"group_points",xgroup_points,compress,append);

}


void Shape_Detector::to_vg(const std::string& filename)
{

    Pwn_vector T(points_begin(),points_end());


	std::ofstream stream(filename, std::ios::out);

	// Part 1.
	// Provides a description of the point cloud
	size_t N = T.size();
	stream << "num_points: " << N << std::endl;
	for (size_t i = 0 ; i < N ; ++i) {
		const Inexact_Point_3 & pt = T[i].first;
        stream << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
    }
    stream << "num_colors: " << 0 << std::endl;
	stream << "num_normals: " << N << std::endl;
	for (size_t i = 0 ; i < N ; ++i) {
		const Inexact_Vector_3 & n = T[i].second;
        stream << n.x() << " " << n.y() << " " << n.z() << std::endl;
    }

	// Part 2.
	// Provides a description of the detected planes.
	size_t M = planes_2.size();
	stream << "num_groups: " << M << std::endl;

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, 255);

    for (size_t i = 0 ; i < planes_2.size() ; ++i) {
		stream << "group_type: 0" << std::endl;
		stream << "num_group_parameters: 4" << std::endl;

		const Inexact_Plane & H = planes_2[i];
		stream << "group_parameters: " << H.a() << " " << H.b() << " " << H.c() << " " << H.d() << std::endl;
		const std::vector<int> & shape_assigned_pts = planes_to_inliers[i];

        stream << "group_label: " << i << std::endl;
        stream << "group_color: " << distribution(generator) << " " << distribution(generator) << " " << distribution(generator) << std::endl;
        stream << "group_num_points: " << shape_assigned_pts.size() << std::endl;

		for (size_t j = 0 ; j < shape_assigned_pts.size() ; ++j) {
			stream << shape_assigned_pts[j] << (j == shape_assigned_pts.size() - 1 ? '\n' : ' ');
		}
		stream << "num_children: 0" << std::endl;
	}
	stream.close();
}


std::pair<double,double> Shape_Detector::get_gauss_sphere_value( const Inexact_Vector_3 & v)
{
	


	double subdivs_90 = 90.0;

    //const double M_PI = 3.1415926535897832384626;

	// We receive a vector v.
	// We assume that it is oriented towards the top of the unit sphere.

	// We compute the longitude and latitude values associated to n.
	// These values belong to the [0, 360] and [0, 90] intervals,
	// which are discretized in (4 * n) and (n + 1) values.

	double f_latitude, f_longitude;

	if (v.x() == 0 && v.y() == 0) {
		std::make_pair(0, 0);
	}
	else {
		double x = CGAL::to_double(v.x()), y = CGAL::to_double(v.y()), z = CGAL::to_double(v.z());
		double r = sqrt(x * x + y * y + z * z);

		f_latitude = acos(z / r);
		if (fabs(f_latitude) < 1e-10) {
			f_longitude = 0;
		}
		else {
			f_longitude = atan2(y / (r * sin(f_latitude)), x / (r * sin(f_latitude)));
            if (f_longitude < 0) f_longitude += 2 * M_PI;
		}
	}

    f_latitude = f_latitude * 180 / M_PI;   // in [0, 90]
    f_longitude = f_longitude * 180 / M_PI; // in [0, 360]

	
	if (round(f_longitude) == round(abs(4 * subdivs_90))) f_longitude = 0;
	return std::make_pair(f_longitude, f_latitude);
}

void Shape_Detector::plane_parellel_group_mean_shift() {
	//get the guass sphere value for each plane
	planes_to_unit_sphere_degree.clear();
	planes_to_inliernumber.clear();
	
	for (int i = 0; i < planes_2.size(); ++i) {
		const Inexact_Plane & P = planes_2[i];
		

		Inexact_Vector_3 u = Inexact_Vector_3(P.a(), P.b(), P.c());
		if (u.z() < 0 || (u.z() == 0 && u.y() < 0) || (u.z() == 0 && u.y() == 0 && u.x() < 0)) u = -u;
		
		planes_to_unit_sphere_degree.push_back(get_gauss_sphere_value(u));
		
		planes_to_inliernumber.push_back(double(planes_to_inliers[i].size())/double(points.size()));
	}



	planes_to_parallel_clusters.clear();
	parallel_clusters_to_planes.clear();
	planes_to_parallel_clusters = std::vector<int>(planes_to_unit_sphere_degree.size(),-1);
	int cluster_id = 0;
	//mean shift to cluster
	MeanShift this_mean_shift;

	for (int i = 0; i < planes_to_unit_sphere_degree.size(); ++i) {
		if (planes_to_parallel_clusters[i] == -1) {
			
			std::pair<double,double> this_center = this_mean_shift.meanshift(planes_to_unit_sphere_degree[i], planes_to_unit_sphere_degree, planes_to_inliernumber, planes_to_parallel_clusters, tolerance_angle, 1);


		
			std::vector<int> this_parallel_planes;
			for (int j = 0; j < planes_to_unit_sphere_degree.size(); ++j) {
				if (planes_to_parallel_clusters[j] != -1) continue;
				double distance_latitude= (planes_to_unit_sphere_degree[j].first - this_center.first)*(planes_to_unit_sphere_degree[j].first - this_center.first) + (planes_to_unit_sphere_degree[j].second - this_center.second)*(planes_to_unit_sphere_degree[j].second - this_center.second);
				
				if (std::sqrt(distance_latitude) < tolerance_angle) {
					planes_to_parallel_clusters[j] = cluster_id;
					
					this_parallel_planes.push_back(j);
				
				}
			}
			if(this_parallel_planes.size()!=0){
			
				
				parallel_clusters_to_planes.push_back(this_parallel_planes);
				cluster_id++;
			}
			
		}
	}

	for (int i = 0; i < planes_to_parallel_clusters.size(); ++i) {

		if (planes_to_parallel_clusters[i] == -1) {
			planes_to_parallel_clusters[i] = parallel_clusters_to_planes.size();
			parallel_clusters_to_planes.push_back(std::vector<int>(1, i));

		}
	}

	
}



void Shape_Detector::initialize_parallel_clusters() {

	plane_parellel_group_mean_shift();
	
	normal_parallel_clusters.clear();
	normal_parallel_clusters.resize(parallel_clusters_to_planes.size());
	


	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		std::vector<int> this_parallel_cluster = parallel_clusters_to_planes[i];
		if (this_parallel_cluster.size() == 1) {

			normal_parallel_clusters[i] = planes_2[this_parallel_cluster[0]].orthogonal_vector();

		}
		else {
			//calculate the new orientation
			std::vector<Inexact_Point_3> this_translated_points;//all the points of a parallel cluster translated to original
			for (int j = 0; j < this_parallel_cluster.size(); ++j) {
				Inexact_Vector_3 this_tranlate_vector = Inexact_Point_3(0, 0, 0) - planes_centroids[this_parallel_cluster[j]];
				for (int cluster_point_id : planes_to_inliers[this_parallel_cluster[j]]) {
						this_translated_points.push_back(points[cluster_point_id].first + this_tranlate_vector);
					}

				}
				Inexact_Plane translated_original_plane;
				linear_least_squares_fitting_3(this_translated_points.begin(), this_translated_points.end(), translated_original_plane, CGAL::Dimension_tag<0>());
				Inexact_Vector_3 this_parallel_cluster_normal = translated_original_plane.orthogonal_vector();

				normal_parallel_clusters[i] = (this_parallel_cluster_normal);

			
		}


	}
}

void Shape_Detector::initialize_orthogonal_clusters_after_parallel() {
	//for each parallel cluster, we check its parallel clusters.
	parallel_id_changed_normal_after_orthogonal.clear();

	parallel_id_changed_normal_after_orthogonal.resize(parallel_clusters_to_planes.size());

	
	
	for (int i = 0; i < parallel_clusters_to_planes.size(); ++i) {
		Inexact_Vector_3 this_normal = normal_parallel_clusters[i];
		
		std::vector<int> this_orthogonal_cluster_ids;
		
		this_orthogonal_cluster_ids.push_back(i);
		for (int j = 0; j < parallel_clusters_to_planes.size(); ++j) {

			Inexact_Vector_3 that_normal = normal_parallel_clusters[j];


			double cos_normal = abs(this_normal*that_normal) / (sqrt(this_normal.squared_length())*sqrt(that_normal.squared_length()));
            if (cos_normal < sin(tolerance_angle * M_PI / 180.0)) {
				
				this_orthogonal_cluster_ids.push_back(j);
			}

		}
		if (this_orthogonal_cluster_ids.size() > 1) {
			

			Inexact_Vector_3 center_normal = normal_parallel_clusters[this_orthogonal_cluster_ids[0]];
			std::vector<std::pair<int, Inexact_Vector_3>> this_prallel_ids_normals;
			this_prallel_ids_normals.push_back(std::make_pair(this_orthogonal_cluster_ids[0], center_normal));
			for (int j = 1; j < this_orthogonal_cluster_ids.size(); ++j) {
				Inexact_Vector_3 this_normal = normal_parallel_clusters[this_orthogonal_cluster_ids[j]];
				Inexact_Vector_3 mid_normal = CGAL::cross_product(center_normal, this_normal);
				Inexact_Vector_3 aligned_this_normal = CGAL::cross_product(center_normal, mid_normal);
				this_prallel_ids_normals.push_back(std::make_pair(this_orthogonal_cluster_ids[j], aligned_this_normal));
			}
			if(this_prallel_ids_normals.size()>1){
				parallel_id_changed_normal_after_orthogonal[i]=this_prallel_ids_normals;//id of parallel cluster.
			}

		}


	}

	




}

void Shape_Detector::initialize_coplanar_clusters_after_parallel() {
	planes_centroids_coplanar.clear();
	planes_centroids_coplanar = planes_centroids;
	parallel_cluster_to_coplanr_cases.clear();
	parallel_cluster_to_coplanr_cases.resize(parallel_clusters_to_planes.size());
	//For each parallel cluster, using mean shift to cluster their primitives according to their signed distance to origin.
	for (size_t j = 0; j < parallel_clusters_to_planes.size(); ++j) {
		//int co_n = 0;
		if (parallel_clusters_to_planes[j].size() == 1) continue;
		Inexact_Point_3 original_p = Inexact_Point_3(0, 0, 0);
		std::vector<double> distances_to_original;
		//signed distance
		std::vector<double> this_group_number_points;
		for (int l = 0; l < parallel_clusters_to_planes[j].size(); ++l) {
			Inexact_Point_3 O_one = planes_centroids[parallel_clusters_to_planes[j][l]];
			const Inexact_Plane & P = Inexact_Plane(O_one, normal_parallel_clusters[j]);
			double dP_norm = 0;
			Inexact_Vector_3 dP(P.projection(original_p) - original_p);
			if (P.has_on_negative_side(original_p)) {
				dP_norm = sqrt(dP.squared_length());

			}
			else {
				dP_norm = -sqrt(dP.squared_length());

			}
			
			distances_to_original.push_back(dP_norm);
			this_group_number_points.push_back(double(planes_to_inliers[parallel_clusters_to_planes[j][l]].size())/double(points.size()));
		}
		//mean shift clustering
		MeanShift this_mean_shift;
		std::vector<int> if_coplanar(parallel_clusters_to_planes[j].size(), -1);
		std::vector<std::vector<int>> coplanar_cased_this_parallel_cluster;

		for (int i = 0; i < parallel_clusters_to_planes[j].size(); ++i) {

			if (if_coplanar[i] == 0)
			{
				continue;
			}


			double this_center = this_mean_shift.meanshift_oned(distances_to_original[i], distances_to_original, this_group_number_points, if_coplanar, tolerance_coplanarity, 1);

			std::vector<int> this_coplanar_id;
			for (int e = 0; e < distances_to_original.size(); ++e) {
				if (if_coplanar[e] != -1) continue;
				double distance_latitude = (distances_to_original[e] - this_center)*(distances_to_original[e] - this_center);
				if (std::sqrt(distance_latitude) < tolerance_coplanarity) {
					if_coplanar[e] = 0;

					this_coplanar_id.push_back(parallel_clusters_to_planes[j][e]);

				}
			}
			//save the new centers if the coplanar operation is conducted.
			if (this_coplanar_id.size() > 1) {
				
				
				coplanar_cased_this_parallel_cluster.push_back(this_coplanar_id);
				double xc = 0, yc = 0, zc = 0;
				double number_p = 0;

				for (int this_id : this_coplanar_id) {
				
					for (size_t j = 0; j < planes_to_inliers[this_id].size(); ++j) {
						const Inexact_Point_3 & pt = points[planes_to_inliers[this_id][j]].first;

						xc += pt.x(), yc += pt.y(), zc += pt.z();
					}
					number_p += planes_to_inliers[this_id].size();
				}
				
				xc /= number_p;
				yc /= number_p;
				zc /= number_p;
				Inexact_Point_3 centroid_common(xc, yc, zc);
				for (int plane_id : this_coplanar_id) {
					planes_centroids_coplanar[plane_id] = centroid_common;

				}
			}


		}
		
		parallel_cluster_to_coplanr_cases[j] = coplanar_cased_this_parallel_cluster;
	}
	

}


Inexact_Point_3 Shape_Detector::get_specific_centroids(int i) {
	double xc = 0, yc = 0, zc = 0;
	double number_p = 0;

		for (int p_id: planes_to_inliers[i]) {
			const Inexact_Point_3 & pt = points[p_id].first;

			xc += pt.x(), yc += pt.y(), zc += pt.z();
		}
		number_p += double(planes_to_inliers[i].size());
	
	xc /= number_p;
	yc /= number_p;
	zc /= number_p;
	Inexact_Point_3 centroid_common(xc, yc, zc);
	return centroid_common;

}


//calculate the L2 changing for a parallel or a coplanar operator.

double Shape_Detector::euclidean_distance_change_after_a_coplanar_or_parallel(std::vector<int> this_parallel_cluster, int id_cluster) {
	
	double dis_original = 0;
	double dis_change = 0;

	for (int id : this_parallel_cluster) {
		Inexact_Plane original_plane = planes_2[id];
		Inexact_Plane changed_plane = Inexact_Plane(planes_centroids_coplanar[id], normal_parallel_clusters[id_cluster]);
		for (int id_point : planes_to_inliers[id]) {
			dis_original += sqrt(CGAL::squared_distance(original_plane, points[id_point].first));
			dis_change += sqrt(CGAL::squared_distance(changed_plane, points[id_point].first));

		}

	}

	return (dis_change - dis_original);



}

double Shape_Detector::normal_deviation_change_after_a_coplanar_or_parallel(std::vector<int> this_parallel_cluster, int id_cluster) {

	double dis_original = 0;
	double dis_change = 0;

	for (int id : this_parallel_cluster) {
		Inexact_Plane original_plane = planes_2[id];
		Inexact_Plane changed_plane = Inexact_Plane(planes_centroids_coplanar[id], normal_parallel_clusters[id_cluster]);
		for (int id_point : planes_to_inliers[id]) {
			dis_original += abs(points[id_point].second * original_plane.orthogonal_vector());
			dis_change += abs(points[id_point].second * changed_plane.orthogonal_vector());

		}

	}

	return (dis_change - dis_original);



}



// euclidean distance changing after a orthogonal operation.
double Shape_Detector::euclidean_distance_change_after_orthogonal(std::vector<std::pair<int, Inexact_Vector_3>> this_orthogonal_cluster) {

	double dis_original = 0;
	double dis_change = 0;

	for (int i = 0; i < this_orthogonal_cluster.size();++i) {
		Inexact_Vector_3 this_normal = this_orthogonal_cluster[i].second;
		int parallel_cluster_id = this_orthogonal_cluster[i].first;
		std::vector<int> planes_ids = parallel_clusters_to_planes[parallel_cluster_id];
		for (int id : planes_ids) {
			Inexact_Plane original_plane = planes_2[id];
			Inexact_Plane changed_plane = Inexact_Plane(planes_centroids_coplanar[id], this_normal);
			for (int id_point : planes_to_inliers[id]) {
				dis_original += sqrt(CGAL::squared_distance(original_plane, points[id_point].first));
				dis_change += sqrt(CGAL::squared_distance(changed_plane, points[id_point].first));

			}
		}
	}

	return (dis_change - dis_original);



}

double Shape_Detector::normal_deviation_change_after_orthogonal(std::vector<std::pair<int, Inexact_Vector_3>> this_orthogonal_cluster) {

	double dis_original = 0;
	double dis_change = 0;

	for (int i = 0; i < this_orthogonal_cluster.size(); ++i) {
		Inexact_Vector_3 this_normal = this_orthogonal_cluster[i].second;
		int parallel_cluster_id = this_orthogonal_cluster[i].first;
		std::vector<int> planes_ids = parallel_clusters_to_planes[parallel_cluster_id];
		for (int id : planes_ids) {
			Inexact_Plane original_plane = planes_2[id];
			Inexact_Plane changed_plane = Inexact_Plane(planes_centroids_coplanar[id], this_normal);
			for (int id_point : planes_to_inliers[id]) {
				dis_original += abs(points[id_point].second * original_plane.orthogonal_vector());
				dis_change += abs(points[id_point].second * changed_plane.orthogonal_vector());

			}
		}
	}

	return (dis_change - dis_original);



}



int Shape_Detector::get_changed_degrees_of_freedom_after_a_regular_operation(int i) {

	int before_dof = get_the_recent_degrees_of_freedom();// dof before this regularity operation conducting.
	std::vector<int> planes_to_parallel_done_after = planes_to_parallel_done;
	std::vector<int> planes_to_orthogonal_done_after = planes_to_orthogonal_done;
	std::vector<int> planes_to_coplanar_done_after = planes_to_coplanar_done;
	std::vector<std::vector<int>> parallel_done_to_planes_after = parallel_done_to_planes;
	std::vector<std::vector<int>> orthogonal_done_to_planes_after = orthogonal_done_to_planes;
	std::vector<std::vector<int>> coplanar_done_to_planes_after = coplanar_done_to_planes;
	//different case form this regularity operation, which is a combination of making exact parallel, orthogonal and coplanar.
	if (parallel_clusters_to_planes[i].size() == 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 && parallel_cluster_to_coplanr_cases[i].size() == 0) {
		//nothing to do.
		
		
	}
	else if (parallel_clusters_to_planes[i].size() > 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 && parallel_cluster_to_coplanr_cases[i].size() == 0) {
		//only parallel operator.
		//recorde the affected primitives in planes_to_parallel_done_after, parallel_done_to_planes_after.

		int id_para_group = parallel_done_to_planes_after.size();
		for (int id : parallel_clusters_to_planes[i]) {
			planes_to_parallel_done_after[id] = id_para_group;

		}
		parallel_done_to_planes_after.push_back(parallel_clusters_to_planes[i]);
		


	}

	
	else if (parallel_clusters_to_planes[i].size() > 1 && parallel_id_changed_normal_after_orthogonal[i].size() == 0 && parallel_cluster_to_coplanr_cases[i].size() > 0) {
		//parallel and coplanar operator.
		//recorde the affected primitives in planes_to_parallel_done_after, parallel_done_to_planes_after, planes_to_coplanar_done_after and coplanar_done_to_planes_after.
		int id_para_group = parallel_done_to_planes_after.size();
		for (int id : parallel_clusters_to_planes[i]) {
			planes_to_parallel_done_after[id] = id_para_group;

		}
		parallel_done_to_planes_after.push_back(parallel_clusters_to_planes[i]);


		
		
		for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[i]) {
			int id_copl_group = coplanar_done_to_planes_after.size();
			coplanar_done_to_planes_after.push_back(a_co);
			for (int id_c : a_co) {
				
				planes_to_coplanar_done_after[id_c] = id_copl_group;

			}
		}
		

	}

	else if (parallel_id_changed_normal_after_orthogonal[i].size() > 0) {
		//orthogonal operator, we should also check each prarllel cluster, if they have coplanar cases.
		//recorde the affected primitives in planes_to_orthogonal_done_after, orthogonal_done_to_planes_after, planes_to_coplanar_done_after and coplanar_done_to_planes_after.

		int id_ortho_group = orthogonal_done_to_planes_after.size();
		std::vector<int> this_ors;
		for (std::pair<int, Inexact_Vector_3> a_or : parallel_id_changed_normal_after_orthogonal[i]) {

			for (int or_id : parallel_clusters_to_planes[a_or.first]) {
				planes_to_orthogonal_done_after[or_id] = id_ortho_group;
				this_ors.push_back(or_id);
			}

		}
		orthogonal_done_to_planes_after.push_back(this_ors);


		for (std::pair<int, Inexact_Vector_3> a_or : parallel_id_changed_normal_after_orthogonal[i]) {
			if (parallel_cluster_to_coplanr_cases[a_or.first].size() > 0) {//coplanar cases exist
				for (std::vector<int> a_co : parallel_cluster_to_coplanr_cases[a_or.first]) {
					int id_copl_group = coplanar_done_to_planes_after.size();
					coplanar_done_to_planes_after.push_back(a_co);
					for (int id_c : a_co) {

						planes_to_coplanar_done_after[id_c] = id_copl_group;

					}
				}
			}
		}

	}
	

	else {
		std::cout << "wrong situation!!!!!!!!!!!" << std::endl;
		getchar();

	}

	//update the regular by planes, some last regular operator will cover the previous one.

	int number_parallel_group = parallel_done_to_planes_after.size();
	int number_orthogonal_group = orthogonal_done_to_planes_after.size();
	int number_coplanar_group = coplanar_done_to_planes_after.size();

	parallel_done_to_planes_after.clear();
	orthogonal_done_to_planes_after.clear();
	coplanar_done_to_planes_after.clear();
	parallel_done_to_planes_after.resize(number_parallel_group);
	orthogonal_done_to_planes_after.resize(number_orthogonal_group);
	coplanar_done_to_planes_after.resize(number_coplanar_group);

	for (int i = 0; i < planes_to_parallel_done_after.size(); ++i) {
		if (planes_to_parallel_done_after[i] != -1) {
			parallel_done_to_planes_after[planes_to_parallel_done_after[i]].push_back(i);
		}
		if (planes_to_orthogonal_done_after[i] != -1) {
			orthogonal_done_to_planes_after[planes_to_orthogonal_done_after[i]].push_back(i);
		}
		if (planes_to_coplanar_done_after[i] != -1) {

			coplanar_done_to_planes_after[planes_to_coplanar_done_after[i]].push_back(i);
		}
	}
	std::vector<std::vector<int>> parallel_done_to_planes_new;
	for (std::vector<int> this_parallel_done : parallel_done_to_planes_after) {
		if (this_parallel_done.size() > 1) {

			parallel_done_to_planes_new.push_back(this_parallel_done);
		}

	}

	if (parallel_done_to_planes_new.size() != parallel_done_to_planes_after.size()) {
		parallel_done_to_planes_after = parallel_done_to_planes_new;
		planes_to_parallel_done_after.clear();
		planes_to_parallel_done_after = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < parallel_done_to_planes_after.size(); ++i) {
			for (int id_p : parallel_done_to_planes_after[i]) {
				planes_to_parallel_done_after[id_p] = i;
			}

		}

	}


	std::vector<std::vector<int>> orthogonal_done_to_planes_new;
	for (std::vector<int> this_orthogonal_done : orthogonal_done_to_planes_after) {
		if (this_orthogonal_done.size() > 1) {

			orthogonal_done_to_planes_new.push_back(this_orthogonal_done);
		}

	}

	if (orthogonal_done_to_planes_new.size() != orthogonal_done_to_planes_after.size()) {
		orthogonal_done_to_planes_after = orthogonal_done_to_planes_new;
		planes_to_orthogonal_done_after.clear();
		planes_to_orthogonal_done_after = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < orthogonal_done_to_planes_after.size(); ++i) {
			for (int id_p : orthogonal_done_to_planes_after[i]) {
				planes_to_orthogonal_done_after[id_p] = i;
			}

		}

	}

	std::vector<std::vector<int>> coplanar_done_to_planes_new;
	for (std::vector<int> this_coplanar_done : coplanar_done_to_planes_after) {
		if (this_coplanar_done.size() > 1) {

			coplanar_done_to_planes_new.push_back(this_coplanar_done);
		}

	}

	if (coplanar_done_to_planes_new.size() != coplanar_done_to_planes_after.size()) {
		coplanar_done_to_planes_after = coplanar_done_to_planes_new;
		planes_to_coplanar_done_after.clear();
		planes_to_coplanar_done_after = std::vector<int>(planes_2.size(), -1);

		for (int i = 0; i < coplanar_done_to_planes_after.size(); ++i) {
			for (int id_p : coplanar_done_to_planes_after[i]) {
				planes_to_coplanar_done_after[id_p] = i;
			}

		}

	}

	// calculate the after dof
	int after_dof = 0;
	for (int i = 0; i < planes_to_coplanar_done_after.size(); ++i) {
		if (planes_to_coplanar_done_after[i] == -1 && planes_to_orthogonal_done_after[i] == -1 && planes_to_parallel_done_after[i] == -1) {
			after_dof += 3;

		}
		else if (planes_to_coplanar_done_after[i] == -1 && planes_to_orthogonal_done_after[i] == -1 && planes_to_parallel_done_after[i] != -1) {
			after_dof += 1;

		}
		else if (planes_to_coplanar_done_after[i] == -1 && planes_to_orthogonal_done_after[i] != -1 && planes_to_parallel_done_after[i] == -1) {
			after_dof += 1;

		}

	}
	after_dof += 3 * orthogonal_done_to_planes_after.size();
	after_dof += 2 * parallel_done_to_planes_after.size();
	after_dof += coplanar_done_to_planes_after.size();

	return after_dof - before_dof;


}

void Shape_Detector::discretize_planes()
{
	
	double dist_min_threshold = discretization_distance;


    // this line seems to be problematic, as it is only triggered when
    // epsilon is not multiplied by the bounding box diagonal
    // (because discretization distance is multiplied)!
    // to solve the problem for now I will keep on multiplying epsilon with the bounding
    // box diagonal
	if (dist_min_threshold > epsilon / 2.0) {

		dist_min_threshold = epsilon / 2.0;

	}
    _logger->debug("dist_min_threshold = {}", dist_min_threshold );

	non_coplanar_planes = 0;


	size_t n = planes_0.size();
	planes_2.clear();
	planes_2.reserve(n);

	coplanarity_table.clear();
	coplanarity_table = std::vector<int>(n, -1);

	// Discretize the orientations of the extracted planes

	std::map<std::pair<int, int>, std::list<int> > atlas;

	double ds = discretization_angle;
	//int subdivs = 90 / discretization_angle;
	int subdivs_90 = int(round(90 / ds));

	for (size_t i = 0 ; i < n ; ++i) {

		const Inexact_Plane & P = planes_1[i];
		Inexact_Point_3 O = planes_centroids[i]; // P.projection(CGAL::ORIGIN);

		Inexact_Vector_3 u = Inexact_Vector_3(P.a(), P.b(), P.c());
		if (u.z() < 0 || (u.z() == 0 && u.y() < 0) || (u.z() == 0 && u.y() == 0 && u.x() < 0)) u = -u;

		int k_longitude, k_latitude;
		discretize_vector(ds, u, k_longitude, k_latitude);

		double a = 0, b = 0, c = 0;
		bool compute_plane_coefficients = true;
		int plane_index = -1;

		// Now, two possible solutions :
		// either there already exists a parallel plane sufficiently close to O in the list of planes,
		// and if so we consider that the current polygon actually belongs to that plane,
		// or such a plane doesn't exist and we register P.
		
		std::pair<int, int> key = std::make_pair(k_longitude, k_latitude);
		std::map<std::pair<int, int>, std::list<int> >::iterator it_key = atlas.find(key);
		
		if (it_key != atlas.end()) {

			const std::list<int> & parallel_planes = it_key->second;
			const Inexact_Plane & P_parallel = planes_2[parallel_planes.front()];

			a = P_parallel.a(), b = P_parallel.b(), c = P_parallel.c();
			compute_plane_coefficients = false;

			double dist_min = dist_min_threshold;

			for (std::list<int>::const_iterator it_p = parallel_planes.begin(); it_p != parallel_planes.end(); it_p++) {
				const Inexact_Plane & Q = planes_2[*it_p];
				Inexact_Vector_3 dP (Q.projection(O) - O);
				double dP_norm = sqrt(dP.squared_length());

				if (dP_norm < dist_min) {
					dist_min = dP_norm;
					plane_index = (*it_p);
				}
			}
		}

		if (plane_index != -1) {
			// Assimilates the polygon to the plane pointed by the previous variable
			planes_2.push_back(planes_2[plane_index]);
			coplanarity_table[i] = plane_index;

		} else {
			// Registers a new plane
			if (compute_plane_coefficients) {

				if (k_latitude == 0) {
					a = 0, b = 0, c = 1;
				} else if (k_latitude == subdivs_90) {
					if (k_longitude == 0 || k_longitude == 4 * subdivs_90) {
						a = 1;
						b = 0;
					} else if (k_longitude == subdivs_90) {
						a = 0;
						b = 1;
					} else if (k_longitude == 2 * subdivs_90) {
						a = -1;
						b = 0;
					} else if (k_longitude == 3 * subdivs_90) {
						a = 0;
						b = -1;
					} else {
                        a = cos(k_longitude * ds * M_PI / 180);
                        b = sin(k_longitude * ds * M_PI / 180);
					}
					c = 0;
				} else if (k_longitude % subdivs_90 == 0) {
					if (k_longitude == 0 || k_longitude == 4 * subdivs_90) {
                        a = sin(k_latitude * ds * M_PI / 180);
						b = 0;
					} else if (k_longitude == subdivs_90) {
						a = 0;
                        b = sin(k_latitude * ds * M_PI / 180);
					} else if (k_longitude == 2 * subdivs_90) {
                        a = -sin(k_latitude * ds * M_PI / 180);
						b = 0;
					} else if (k_longitude == 3 * subdivs_90) {
						a = 0;
                        b = -sin(k_latitude * ds * M_PI / 180);
					}
                    c = cos(k_latitude * ds * M_PI / 180);
				} else {
                    a = sin(k_latitude * ds * M_PI / 180) * cos(k_longitude * ds * M_PI / 180);
                    b = sin(k_latitude * ds * M_PI / 180) * sin(k_longitude * ds * M_PI / 180);
                    c = cos(k_latitude * ds * M_PI / 180);
				}

				if (c < 0 || (c == 0 && b < 0) || (c == 0 && b == 0 && a < 0)) {
					if (a != 0) a = -a;
					if (b != 0) b = -b;
					if (c != 0) c = -c;
				}
			}

			double d = -(a * O.x() + b * O.y() + c * O.z());
			Inexact_Plane P(a, b, c, d);
			
			plane_index = int(planes_2.size());
			++non_coplanar_planes;

			planes_2.push_back(P);
			atlas[key].push_back(plane_index);
		}
	}


}

void Shape_Detector::make_3d_histogram(std::map<std::pair<int, int>, std::list<int> > & A, const double ds, std::vector<std::pair<Inexact_Vector_3, double> > & H)
{
	for (auto it_a = A.begin() ; it_a != A.end() ; ++it_a) 
	{
		// Gets a key
		int k_longitude = it_a->first.first, k_latitude = it_a->first.second;
        double a = sin(k_latitude * ds * M_PI / 180) * cos(k_longitude * ds * M_PI / 180);
        double b = sin(k_latitude * ds * M_PI / 180) * sin(k_longitude * ds * M_PI / 180);
        double c = cos(k_latitude * ds * M_PI / 180);
		Inexact_Vector_3 u(a, b, c);

		// Gets a list of planes associated to a certain orientation
		// Computes a sum of areas or any quantity of similar nature
		const std::list<int> & P = it_a->second;
		double s = 0;
		for (int p : P) {			
			s += planes_to_inliers[p].size();
		}

		// Inserts the final result
		H.push_back(std::make_pair(u, s));		
	}
}

void Shape_Detector::print_3d_histogram(const std::string & filename, std::vector<std::pair<Inexact_Vector_3, double> > & H)
{
	std::ofstream stream(filename);
	if (!stream.is_open()) return;

	stream << H.size() << std::endl;
	for (const std::pair<Inexact_Vector_3, double> & h : H) {
		const Inexact_Vector_3 & u = h.first;
		const double & s = h.second;
		stream << u.x() << " " << u.y() << " " << u.z() << " " << s << std::endl;
	}

	stream.close();
}

void Shape_Detector::discretize_vector(const double ds, const Inexact_Vector_3 & v, int & k_longitude, int & k_latitude)
{
	int subdivs_90 = int(round(90 / ds));

	// We receive a vector v.
	// We assume that it is oriented towards the top of the unit sphere.

	// We compute the longitude and latitude values associated to n.
	// These values belong to the [0, 360] and [0, 90] intervals,
	// which are discretized in (4 * n) and (n + 1) values.

	double f_latitude, f_longitude;

	if (v.x() == 0 && v.y() == 0) {
		k_latitude = k_longitude = 0;
		return;
	} else {
		double x = CGAL::to_double(v.x()), y = CGAL::to_double(v.y()), z = CGAL::to_double(v.z());
		double r = sqrt(x * x + y * y + z * z);

		f_latitude = acos(z / r);
		if (fabs(f_latitude) < 1e-10) {
			f_longitude = 0;
		} else {
			f_longitude = atan2(y / (r * sin(f_latitude)), x / (r * sin(f_latitude)));
            if (f_longitude < 0) f_longitude += 2 * M_PI;
		}
	}

    f_latitude = f_latitude * 180 / M_PI;   // in [0, 90]
    f_longitude = f_longitude * 180 / M_PI; // in [0, 360]

	// Discretizes

	int i_inf = floor(f_latitude / ds), i_sup = ceil(f_latitude / ds);
	double dist_i_inf = fabs(i_inf * ds - f_latitude), dist_i_sup = fabs(i_sup * ds - f_latitude);
	k_latitude = (dist_i_inf < dist_i_sup ? i_inf : i_sup);
	if (k_latitude == 0) {
		k_longitude = 0;
	} else {
		int j_inf = floor(f_longitude / ds), j_sup = ceil(f_longitude / ds);
		double dist_j_inf = fabs(j_inf * ds - f_longitude), dist_j_sup = fabs(j_sup * ds - f_longitude);
		k_longitude = (dist_j_inf < dist_j_sup ? j_inf : j_sup);

		if (k_longitude == 4 * subdivs_90) k_longitude = 0;
	}
}

int Shape_Detector::get_number_of_non_coplanar_planes()
{
	return non_coplanar_planes;
}

Pwn_vector::const_iterator Shape_Detector::points_begin() const
{
	return points.cbegin();
}

Pwn_vector::const_iterator Shape_Detector::points_end() const
{
	return points.cend();
}
