cmake_minimum_required(VERSION 3.14.0)
project(psdr_example LANGUAGES CXX)

######################################################
###################### SETTINGS ######################
######################################################
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


#######################################################
####################### EXTERNAL ######################
#######################################################
find_package(CGAL QUIET REQUIRED)
include(${CGAL_USE_FILE})
### issues with eigen? Do this: https://stackoverflow.com/a/25537967


find_package(Boost REQUIRED COMPONENTS system filesystem thread serialization program_options)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(spdlog REQUIRED)
find_package(xtl REQUIRED) # for xtensor
find_package(ZLIB REQUIRED) # for xtensor
find_package(xtensor REQUIRED)
find_package(xtensor-io REQUIRED)


include_directories(../../psdr/include)
# make sure to include psdr build which was generated with the same environment
link_directories(../../psdr/build/scikit)
#link_directories(../../psdr/build/release)

##############################
############ psdr ############
##############################
add_executable(example
    "example.cpp")

target_link_libraries(example ${Boost_LIBRARIES})
target_link_libraries(example ${CGAL_LIBRARIES})
target_link_libraries(example spdlog::spdlog)
target_link_libraries(example ZLIB::ZLIB xtensor xtensor-io)
target_link_libraries(example ShapeDetection)

