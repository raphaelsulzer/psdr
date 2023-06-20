#include <iostream>
#include <chrono>
#include <thread>

// EXTERNAL
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "yaml-cpp/yaml.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

// INTERNAL
#include "input_parser.h"
#include "shape_container.h"
#include "shape_detector.h"


using namespace std;
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

//    Parser ip;
//    ip.parse(argc,argv);

    string yaml_file;
    if(argc > 1)
        yaml_file = argv[1];
    else
        yaml_file = "/home/rsulzer/data/reconbench/confs/anchor_1.yaml";
    YAML::Node config = YAML::LoadFile(yaml_file);

    logger->info("Detect planar shapes...");
    logger->debug("Load pointcloud from {}",config["data"]["pointcloud"].as<string>());

    auto tsd = Shape_Detector();

    auto SD = Shape_Detector();
    if (SD.load_points(config["data"]["pointcloud"].as<string>())){
        logger->error("Could not load input pointcloud!");
        return 1;
    }

    double pd_epsilon = config["plane_detection"]["epsilon"].as<double>();
    int min_inliers = config["plane_detection"]["min_inliers"].as<int>();
    double normal_th = config["plane_detection"]["normal_th"].as<double>();
    int knn = config["plane_detection"]["knn"].as<int>();

    logger->info("Detect planes with epsilon {}, min inliers {} and normal threshold {}",pd_epsilon,min_inliers,normal_th);

    SD.set_detection_parameters(min_inliers, pd_epsilon, knn, normal_th);

    auto SC = Shape_Container(&SD);
    if(SC.detect()){
        logger->error("Could not detect planes");
        return 1;
    }
    if(config["plane_detection"]["refine"].as<int>()){
        if(SC.refine()){
            logger->error("Could not refine planes!");
            return 1;
        }
    }

    auto plane_root = fs::path(config["data"]["planes"].as<string>()).parent_path();
    if(!fs::is_directory(plane_root)) fs::create_directories(plane_root);
    logger->debug("save planes to {}",config["data"]["planes"].as<string>());
    string path_vg = fs::path(config["data"]["planes"].as<string>()).replace_extension(".vg").string();
    SC.save(path_vg);
    string path_npz = fs::path(config["data"]["planes"].as<string>()).replace_extension(".npz").string();
    logger->debug("save planes to {}",path_npz);
    SC.save(path_npz);
    string path_ply = fs::path(config["data"]["planes"].as<string>()).replace_extension(".ply").string();
    logger->debug("save planes to {}",path_ply);
    SC.save(path_ply, "convex");

    return 0;
}
