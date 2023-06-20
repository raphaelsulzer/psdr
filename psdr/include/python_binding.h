#pragma once

#include <random>

#include "shape_container.h"
#include "shape_detector.h"

using namespace std;

class pyPSDR
{
public:
//    std::shared_ptr<spdlog::logger> _logger;
    pyPSDR(int verbosity);
//    ~pyPSDR();

    Shape_Detector _SD;
    Shape_Container _SC;

    int load_points(const string filename);
    int detect(int rg_min_points, double rg_epsilon, double rg_normal_threshold, int knn);
    int refine(int max_iter);
    int save(const string filename, const string type);
};


