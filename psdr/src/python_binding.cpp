//#include <iostream>
#include <python_binding.h>
//#include <random>
//#include <sys/stat.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>


#include "shape_container.h"
#include "shape_detector.h"

using namespace std;
namespace nb = nanobind;
using namespace nb::literals;

pyPSDR::pyPSDR(int verbosity = 0){

    _SD = Shape_Detector();
    _SC = Shape_Container(&_SD);

    if(verbosity == 0)
        _SD._logger->set_level(spdlog::level::off);
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
            .def("load_points", &pyPSDR::load_points, "file"_a ,"Load a point cloud with normals or a point group file.")
            .def("detect", &pyPSDR::detect, "min_inliers"_a = 25, "epsilon"_a = 0.2, "normal_th"_a = 0.85, "knn"_a = 10, "Detect planar shapes.")
            .def("refine", &pyPSDR::refine, "max_iter"_a = -1, "Refine planar shapes.")
            .def("save", &pyPSDR::save, "file"_a, "primitive_type"_a = "convex", "Save planar shapes to file.");
}


