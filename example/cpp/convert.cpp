#include <iostream>
#include <chrono>
#include <thread>

// EXTERNAL
#include <boost/filesystem.hpp>
// #include <boost/program_options.hpp>
// #include "yaml-cpp/yaml.h"
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
    auto SC = Shape_Container(&SD);

    auto plane_file = fs::path(pointcloud_file).replace_extension(".ply").string();
    logger->debug("save planes to {}",plane_file);
    SC.save(plane_file, "convex");
    return 0;
}
