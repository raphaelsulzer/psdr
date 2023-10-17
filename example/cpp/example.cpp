#include <iostream>
#include <chrono>
#include <thread>

// EXTERNAL
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "yaml-cpp/yaml.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

//// INTERNAL
#include "shape_container.h"
#include "shape_detector.h"


//using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;



int main(int argc, char const *argv[]){

    auto logger = spdlog::stdout_color_mt("PSDR");
    #ifdef NDEBUG
     // nondebug
    #else
    spdlog::set_level(spdlog::level::debug);
    #endif
    spdlog::set_pattern("[%H:%M:%S] [%n] [%l] %v");


    string pointcloud_file;
    if(argc > 1)
        pointcloud_file = argv[1];
    else
        pointcloud_file = "../../../anchor.ply";

    logger->info("Detect planar shapes...");
    logger->debug("Load pointcloud from {}",pointcloud_file);


    auto SD = Shape_Detector();
    if (SD.load_points(pointcloud_file)){
        logger->error("Could not load input pointcloud!");
        return 1;
    }



    double pd_epsilon = 0.02*SD.get_bbox_diagonal();
    int min_inliers = 50;
    double normal_th = 0.85;
    int knn = 10;

    logger->info("Detect planes with epsilon {}, min inliers {} and normal threshold {}",pd_epsilon,min_inliers,normal_th);


    SD.set_detection_parameters(min_inliers, pd_epsilon, knn, normal_th);


    auto SC = Shape_Container(&SD);
    if(SC.detect()){
        logger->error("Could not detect planes");
        return 1;
    }

    string plane_file = "/home/rsulzer/cpp/psdr/example/data/anchor/convexes_detected/file.ply";
    logger->debug("save planes to {}",plane_file);
    SC.save(plane_file, "convex");
    plane_file = "/home/rsulzer/cpp/psdr/example/data/anchor/convexes_detected/file.npz";
    SC.save(plane_file);

    // refinement params
    SD.set_constraint(false);
    SD.set_weight_m(0);
    SD.set_metric_type(0);
    SD.set_lambda_c(1.0);
    SD.set_lambda_r(1.0);
    SD.set_lambda_regularity(1.0);
    SD.set_lambda_fidelity(1.0);

    SD.set_regularization_parameters(5.0, pd_epsilon/2.0); // not called anywhere

    SD.set_discretization_parameters(0.5, pd_epsilon/2.0);

    if(SC.refine(2)){
        logger->error("Could not refine planes!");
        return 1;
    }


    plane_file = "/home/rsulzer/cpp/psdr/example/data/anchor/convexes_refined/file.ply";
    logger->debug("save planes to {}",plane_file);
    SC.save(plane_file, "convex");
    plane_file = "/home/rsulzer/cpp/psdr/example/data/anchor/convexes_refined/file.npz";
    SC.save(plane_file);

    return 0;
}
