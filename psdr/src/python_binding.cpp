//#include <iostream>
#include <python_binding.h>
//#include <random>
//#include <sys/stat.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>


#include "shape_container.h"
#include "shape_detector.h"

using namespace std;
namespace nb = nanobind;
using namespace nb::literals;


pyPSDR::pyPSDR(int verbosity){

    _SD = Shape_Detector();
    _SC = Shape_Container(&_SD);

    if(verbosity == 0)
        _SD._logger->set_level(spdlog::level::warn);
    else if(verbosity == 1)
        _SD._logger->set_level(spdlog::level::info);
    else if(verbosity == 2)
        _SD._logger->set_level(spdlog::level::debug);
    else
        _SD._logger->set_level(spdlog::level::off);

}


int pyPSDR::load_points(const string filename){
    return _SD.load_points(filename);
}


int pyPSDR::load_points(const array3& points, const array3& normals, const array1& classes){

    _SD.points.clear();
    _SD.point_classes.clear();
    _SD.points_if_added.clear();
    _SD.planes.clear();
    _SD.planes_0.clear();
    _SD.planes_1.clear();
    _SD.planes_2.clear();

    assert(points.shape(0) == normals.shape(0));
    assert(points.shape(0) == classes.shape(0));

    // this is later used to determine if planes should be detected or simply read from existing file, here we want them to be detected
    _SD.path_point_cloud_extension = ".ply";


    for(size_t i = 0; i < points.shape(0); i++){
        Point_with_normal pwn;
        pwn.first = Inexact_Point_3(points(i,0),points(i,1),points(i,2));
        pwn.second = Inexact_Vector_3(normals(i,0),normals(i,1),normals(i,2));
        _SD.points.push_back(pwn);
        _SD.point_classes.push_back(classes(i));
    }

    _SD._logger->debug("Loaded {} points",points.shape(0));

    _SD.if_oriented_normal = true;

    _SD.inliers_to_natural_colors = std::vector<CGAL::Color>(points.size(), CGAL::black());
    _SD.spacing_is_known = false;
    _SD.set_extrema();

    return 0;

}

int pyPSDR::load_points(const array3& points, const array3& normals){

    _SD.points.clear();
    _SD.points_if_added.clear();
    _SD.planes.clear();
    _SD.planes_0.clear();
    _SD.planes_1.clear();
    _SD.planes_2.clear();

    assert(points.shape(0) == normals.shape(0));

    // this is later used to determine if planes should be detected or simply read from existing file, here we want them to be detected
    _SD.path_point_cloud_extension = ".ply";


    for(size_t i = 0; i < points.shape(0); i++){
        Point_with_normal pwn;
        pwn.first = Inexact_Point_3(points(i,0),points(i,1),points(i,2));
        pwn.second = Inexact_Vector_3(normals(i,0),normals(i,1),normals(i,2));
        _SD.points.push_back(pwn);
    }

    _SD._logger->debug("Loaded {} points",points.shape(0));

    _SD.if_oriented_normal = true;

    _SD.inliers_to_natural_colors = std::vector<CGAL::Color>(points.size(), CGAL::black());
    _SD.spacing_is_known = false;
    _SD.set_extrema();



    return 0;

}


int pyPSDR::load_points(const array3& points){

    _SD.points.clear();

//    _logger->warn("input .ply file does not contain valid 'normals'. 'normals' will be estimated.");

//    _logger->debug("estimate global k neighbor scale");
//    // first estimate k for knn
//    vector<Inexact_Point_3> tpoints;
//    for(const auto pnv : points)
//        tpoints.push_back(pnv.first);
//    knn = CGAL::estimate_global_k_neighbor_scale(tpoints);
//    should_compute_knn = false;

//    _logger->debug("estimate normals using jets");

//    CGAL::jet_estimate_normals<Concurrency_tag>
//        (points,
//            knn, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
//            CGAL::parameters::point_map(Point_map())
//            .normal_map(Normal_map()));
//    if_oriented_normal = false;

    return 0;

}




int pyPSDR::detect(int rg_min_points = 25, double rg_epsilon = 0.2, double rg_normal_threshold = 0.85, int knn = 10){

    _SD.set_detection_parameters(rg_min_points, rg_epsilon, knn, rg_normal_threshold);

    _SD.set_regularization_parameters(0.8, 5, rg_epsilon/2.0);
    _SD.set_discretization_parameters(0.5, 0.4);
    _SD.set_constraint(false);
    _SD.set_weight_m(0);
    _SD.set_metric_type(0);
    _SD.set_lambda_fidelity(1.0);
    _SD.set_lambda_c(1.0);
    _SD.set_lambda_r(1.0);
    _SD.set_lambda_regularity(1.0);

    return _SC.detect();
}

int pyPSDR::refine(int max_iter = -1){
    return _SC.refine(max_iter);
}

int pyPSDR::save(const string filename, const string type = "convex"){
    return _SC.save(filename, type);
}


NB_MODULE(pypsdr_ext, m) {
    nb::class_<pyPSDR>(m, "psdr")
            .def(nb::init<int>(),"verbosity"_a = 0)
            .def("load_points", nb::overload_cast<const string>(&pyPSDR::load_points), "file"_a ,"Load a point cloud with normals or a point group file.")
            .def("load_points", nb::overload_cast<const array3&, const array3&, const array1&>(&pyPSDR::load_points), "points"_a, "normals"_a, "classes"_a, "Load points and normals from numpy arrays.")
            .def("load_points", nb::overload_cast<const array3&, const array3&>(&pyPSDR::load_points), "points"_a, "normals"_a, "Load points and normals from numpy arrays.")
            .def("load_points", nb::overload_cast<const array3&>(&pyPSDR::load_points), "points"_a ,"Load points from a numpy array.")
            .def("detect", &pyPSDR::detect, "min_inliers"_a = 25, "epsilon"_a = 0.2, "normal_th"_a = 0.85, "knn"_a = 10, "Detect planar shapes.")
            .def("refine", &pyPSDR::refine, "max_iter"_a = -1, "Refine planar shapes.")
            .def("save", &pyPSDR::save, "file"_a, "primitive_type"_a = "convex", "Save planar shapes to file.");
}


