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

    double pd_epsilon = 0.2;
    int min_inliers = 20;
    double normal_th = 0.85;
    int knn = 10;

    logger->info("Detect planes with epsilon {}, min inliers {} and normal threshold {}",pd_epsilon,min_inliers,normal_th);

    SD.set_detection_parameters(min_inliers, pd_epsilon, knn, normal_th);

    auto SC = Shape_Container(&SD);
    if(SC.detect()){
        logger->error("Could not detect planes");
        return 1;
    }

    if(SC.refine()){
        logger->error("Could not refine planes!");
        return 1;
    }


    string plane_file = "/home/rsulzer/cpp/psdr/example/data/bunny/convexes_detected/file.ply";
    logger->debug("save planes to {}",plane_file);
    SC.save(plane_file, "convex");
    plane_file = "/home/rsulzer/cpp/psdr/example/data/bunny/convexes_detected/file.npz";
    SC.save(plane_file);

    return 0;
}
