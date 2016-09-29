#
# CMakeLists.txt for cellsim_34_vcl
#

cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(cellsim CXX)

# enable testing (disable with -DBUILD_TESTING=OFF)
include(CTest)

# add cmake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# add options (builds serial 3 and 4 variable version by default)
option(BUILD_SERIAL "Build serial version" ON)
option(BUILD_CUDA "Build CUDA ViennaCL version" OFF)
option(BUILD_MKL "Build MKL version" OFF)
option(BUILD_OPENMP "Build OpenMP ViennaCL version" OFF)
option(THREE_VARIABLES "Build three variable version" ON)
option(FOUR_VARIABLES "Build four variable version" ON)
option(TEST_SKIP_COMPARE "Skip the comparison part of the tests" OFF)
option(TEST_PYTHON_REDUCE "Add tests for the python version of reduce" ON)

# set CXX standard
add_definitions("-std=c++11")

# detect system for naming executable
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # we are on a mac
    set(platform_suffix "mac")
elseif(${UNIX})
    # differentiate between pan and generic linux
    if(IS_DIRECTORY /gpfs1m/apps AND IS_DIRECTORY /gpfs1m/projects)
        set(platform_suffix "pan")
    else()
        set(platform_suffix "linux")
    endif(IS_DIRECTORY /gpfs1m/apps AND IS_DIRECTORY /gpfs1m/projects)
else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(platform_suffix "unknown")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
message(STATUS "Platform suffix is: '${platform_suffix}'")

# detect compiler for naming executable
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" compiler_suffix)
message(STATUS "Compiler suffix is: '${compiler_suffix}'")

# Set RELEASEHOST build:
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type detected, defaulting to RELEASE")
    set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

# hint for Boost location on PAN
if(NOT BOOST_ROOT)
    set(BOOST_ROOT "/share/easybuild/RHEL6.3/sandybridge/software/Boost/1.58.0-intel-2015a-Python-2.7.9")
endif(NOT BOOST_ROOT)

# find Boost
find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# find Eigen3
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# find ViennaCL
if(${BUILD_SERIAL} OR ${BUILD_OPENMP} OR ${BUILD_CUDA})
    find_package(ViennaCL REQUIRED)
    include_directories(SYSTEM ${VIENNACL_INCLUDE_DIRS})
endif(${BUILD_SERIAL} OR ${BUILD_OPENMP} OR ${BUILD_CUDA})

# find OpenMP
if(${BUILD_OPENMP})
    find_package(OpenMP REQUIRED)
endif(${BUILD_OPENMP})

# find CUDA
if(${BUILD_CUDA})
    find_package(CUDA REQUIRED)
endif(${BUILD_CUDA})

# find MKL
if(${BUILD_MKL})
    set(MKL_SDL ON)
    find_package(MKL REQUIRED)
    include_directories(SYSTEM ${MKL_INCLUDE_DIRS})
endif(${BUILD_MKL})

# if intel compiler, modify flags
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    # force precise floating point calculations
    add_definitions(-fp-model precise)
    
    # add host optimisation to RELEASEHOST build
    set(CMAKE_CXX_FLAGS_RELEASEHOST "-O3 -xHost -DNDEBUG")
    
    # set flags for custom release on pan (optimal on wm and avx)
    set(CMAKE_CXX_FLAGS_RELEASEPAN "-O3 -xSSE4.2 -axAVX -DNDEBUG")
# modify flags on GNU compiler too
elseif(${CMAKE_COMPILER_IS_GNUCXX})
    # add host optimisation to RELEASEHOST build
    set(CMAKE_CXX_FLAGS_RELEASEHOST "-O3 -march=native -mtune=native -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASEPAN "-O3 -march=westmere -mtune=sandybridge -DNDEBUG")
# set defaults
else(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    # set defaults
    set(CMAKE_CXX_FLAGS_RELEASEPAN "${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELEASEHOST "${CMAKE_CXX_FLAGS_RELEASE}")
endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")

# add src directory
add_subdirectory(src)

# add test directory
add_subdirectory(test)
