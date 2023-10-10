#pragma once

using namespace std;

#include <random>

#include "shape_container.h"
#include "shape_detector.h"


#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>
namespace nb = nanobind;
using namespace nb::literals;

typedef nb::ndarray<double, nb::shape<nb::any, 3>> array3;
typedef nb::ndarray<int, nb::shape<nb::any>> array1;

class pyPSDR
{
public:
//    std::shared_ptr<spdlog::logger> _logger;
    pyPSDR(int verbosity = 0);
//    ~pyPSDR();

    Shape_Detector _SD;
    Shape_Container _SC;

    int load_points(const string filename);
    int load_points(const array3& points, const array3& normals, const array1& classes);
    int load_points(const array3& points, const array3& normals);
    int load_points(const array3& points);

    int detect(int rg_min_points, double rg_epsilon, double rg_normal_threshold, int knn);
    int refine(int max_iter);
    int save(const string filename, const string type);
};


