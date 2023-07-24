cmake_minimum_required(VERSION 3.14.0)
project(planar_shape_detection_and_refinement LANGUAGES CXX)

######################################################
###################### SETTINGS ######################
######################################################
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

######################################################
###################### EXTERNAL ######################
######################################################
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
## issues with eigen? Do this: https://stackoverflow.com/a/25537967

if(NOT SKBUILD)
    set(CONDA_ENV "ksr3")
    set(Python_ROOT_DIR "/opt/miniconda3/envs/${CONDA_ENV}/bin")
endif()
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
include_directories(${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS system filesystem thread serialization program_options)

# Run `nanobind.cmake_dir()` from Python to detect where nanobind is installed
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

find_package(Eigen3 REQUIRED PATHS $CONDA_PREFIX/include/eigen3)
message("Eigen include dir: " ${EIGEN3_INCLUDE_DIR})
find_package(yaml-cpp REQUIRED)
find_package(spdlog REQUIRED)
find_package(xtl REQUIRED) # for xtensor
find_package(ZLIB REQUIRED) # for xtensor
find_package(xtensor REQUIRED)
find_package(xtensor-io REQUIRED)


######################################################
###################### INTERNAL ######################
######################################################
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")


###########################################
############ libShapeDetection ############
###########################################
add_library(ShapeDetection STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shape_container.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shape_detector.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shape_detector_index_map.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/MeanShift.cpp"
)
target_link_libraries(ShapeDetection ${Boost_LIBRARIES})
target_link_libraries(ShapeDetection ${CGAL_LIBRARIES})
target_link_libraries(ShapeDetection spdlog::spdlog)
target_link_libraries(ShapeDetection ZLIB::ZLIB xtensor xtensor-io)

##########################################
############# python bindings ############
##########################################
nanobind_add_module(pypsdr_ext "${CMAKE_CURRENT_SOURCE_DIR}/src/python_binding.cpp")
if(SKBUILD)
    include_directories("${CMAKE_CURRENT_BINARY_DIR}")
    target_link_libraries(pypsdr_ext PUBLIC ShapeDetection spdlog::spdlog)
    install(TARGETS pypsdr_ext LIBRARY DESTINATION pypsdr)
endif()
