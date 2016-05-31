#
# CMakeLists.txt for cellsim_34_vcl
#

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(cellsim)

# add cmake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# add options (builds serial 3 and 4 variable version by default)
option(BUILD_SERIAL "Build serial version" ON)
option(BUILD_CUDA "Build CUDA ViennaCL version" OFF)
option(BUILD_MKL "Build MKL version" OFF)
option(BUILD_OPENMP "Build OpenMP ViennaCL version" OFF)
option(THREE_VARIABLES "Build three variable version" ON)
option(FOUR_VARIABLES "Build four variable version" ON)

# set CXX standard
add_definitions("-std=c++11")

# Set release build:
if(NOT CMAKE_BUILD_TYPE)
    message("-- No build type detected, defaulting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# source
file(GLOB cpu_src ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB all_headers ${PROJECT_SOURCE_DIR}/src/*.h)
file(GLOB cuda_src ${PROJECT_SOURCE_DIR}/src/*.cu)

# hint for Boost location on PAN
if(NOT BOOST_ROOT)
    set(BOOST_ROOT "/share/easybuild/RHEL6.3/sandybridge/software/Boost/1.58.0-intel-2015a-Python-2.7.9")
endif(NOT BOOST_ROOT)

# find Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# find Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# find ViennaCL
if(${BUILD_SERIAL} OR ${BUILD_OPENMP} OR ${BUILD_CUDA})
    find_package(ViennaCL REQUIRED)
    include_directories(${VIENNACL_INCLUDE_DIRS})
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
    include_directories(${MKL_INCLUDE_DIRS})
endif(${BUILD_MKL})

# if intel compiler, add host optimisation to release build
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    # force precise floating point calculations
    add_definitions(-fp-model precise)
    # add host optimisation
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -DNDEBUG")
endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")

# serial version (ViennaCL)
if(${BUILD_SERIAL})
    # three variable
    if(${THREE_VARIABLES})
        message("-- Building ViennaCL serial three variable version")
        add_executable(generic3d_03 ${cpu_src} ${all_headers})
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        message("-- Building ViennaCL serial four variable version")
        add_executable(generic3d_04 ${cpu_src} ${all_headers})
        set_property(TARGET generic3d_04 APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_SERIAL})

# OpenMP version (ViennaCL)
if(${BUILD_OPENMP})
    # three variable
    if(${THREE_VARIABLES})
        message("-- Building ViennaCL OpenMP three variable version")
        add_executable(generic3d_03_openmp ${cpu_src} ${all_headers})
        set_property(TARGET generic3d_03_openmp APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET generic3d_03_openmp APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET generic3d_03_openmp APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        message("-- Building ViennaCL OpenMP four variable version")
        add_executable(generic3d_04_openmp ${cpu_src} ${all_headers})
        set_property(TARGET generic3d_04_openmp APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET generic3d_04_openmp APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET generic3d_04_openmp APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
        set_property(TARGET generic3d_04_openmp APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_OPENMP})

# CUDA version (ViennaCL)
if(${BUILD_CUDA})
    # if intel is selected we force gcc instead
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        message("-- WARNING: trying to compile CUDA version with Intel compilers (this probably won't work!)")
    endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    
    # three variable
    if(${THREE_VARIABLES})
        message("-- Building ViennaCL CUDA three variable version")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA")
        cuda_add_executable(generic3d_03_cuda ${cuda_src} ${all_headers}
            OPTIONS
                DEBUG ${DEFINITIONS_STR}
                RELEASE ${DEFINITIONS_STR}
                MINSIZEREL ${DEFINITIONS_STR}
                RELWITHDEBINFO ${DEFINITIONS_STR}
        )
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        message("-- Building ViennaCL CUDA four variable version")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA -DFOUR_VARIABLES")
        cuda_add_executable(generic3d_04_cuda ${cuda_src} ${all_headers}
            OPTIONS
                DEBUG ${DEFINITIONS_STR}
                RELEASE ${DEFINITIONS_STR}
                MINSIZEREL ${DEFINITIONS_STR}
                RELWITHDEBINFO ${DEFINITIONS_STR}
        )
    endif(${FOUR_VARIABLES})
endif(${BUILD_CUDA})

# MKL version
if(${BUILD_MKL})
    # three variable
    if(${THREE_VARIABLES})
        message("-- Building MKL three variable version")
        add_executable(generic3d_03_mkl ${cpu_src} ${all_headers})
        target_link_libraries(generic3d_03_mkl ${MKL_LIBRARIES})
        set_property(TARGET generic3d_03_mkl APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        message("-- Building MKL four variable version")
        add_executable(generic3d_04_mkl ${cpu_src} ${all_headers})
        target_link_libraries(generic3d_04_mkl ${MKL_LIBRARIES})
        set_property(TARGET generic3d_04_mkl APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
        set_property(TARGET generic3d_04_mkl APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_MKL})
