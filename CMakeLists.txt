#
# CMakeLists.txt for cellsim_34_vcl
#

cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(cellsim CXX)

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

# detect system for naming executable
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # we are on a mac
    set(platform_suffix "_mac")
elseif(${UNIX})
    # differentiate between pan and generic linux
    if(IS_DIRECTORY /gpfs1m/apps AND IS_DIRECTORY /gpfs1m/projects)
        set(platform_suffix "_pan")
    else()
        set(platform_suffix "_linux")
    endif(IS_DIRECTORY /gpfs1m/apps AND IS_DIRECTORY /gpfs1m/projects)
else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(platform_suffix "_unknown")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
message(STATUS "Platform suffix is: '${platform_suffix}'")

# detect compiler for naming executable
string(TOLOWER "_${CMAKE_CXX_COMPILER_ID}" compiler_suffix)
message(STATUS "Compiler suffix is: '${compiler_suffix}'")

# Set release build:
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type detected, defaulting to Release")
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
    # add host optimisation
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -DNDEBUG")
    # set flags for custom release on pan (optimal on wm and avx)
    set(CMAKE_CXX_FLAGS_RELPAN "-O3 -xSSE4.2 -axAVX -DNDEBUG")
endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")

# serial version (ViennaCL)
if(${BUILD_SERIAL})
    # three variable
    if(${THREE_VARIABLES})
        set(target_name "generic3d_03${platform_suffix}_vclserial${compiler_suffix}")
        message(STATUS "Building ViennaCL serial three variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_name "generic3d_04${platform_suffix}_vclserial${compiler_suffix}")
        message(STATUS "Building ViennaCL serial four variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_SERIAL})

# OpenMP version (ViennaCL)
if(${BUILD_OPENMP})
    # three variable
    if(${THREE_VARIABLES})
        set(target_name "generic3d_03${platform_suffix}_vclopenmp${compiler_suffix}")
        message(STATUS "Building ViennaCL OpenMP three variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_name "generic3d_04${platform_suffix}_vclopenmp${compiler_suffix}")
        message(STATUS "Building ViennaCL OpenMP four variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_OPENMP})

# CUDA version (ViennaCL)
if(${BUILD_CUDA})
    # if intel is selected we exit
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        message(FATAL_ERROR "Cannot currently build CUDA version with Intel compilers (use GCC instead)!")
    endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    
    # three variable
    if(${THREE_VARIABLES})
        set(target_name "generic3d_03${platform_suffix}_vclcuda${compiler_suffix}")
        message(STATUS "Building ViennaCL CUDA three variable version: ${target_name}")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA")
        cuda_add_executable(${target_name} ${cuda_src} ${all_headers}
            OPTIONS
                DEBUG ${DEFINITIONS_STR}
                RELEASE ${DEFINITIONS_STR}
                MINSIZEREL ${DEFINITIONS_STR}
                RELWITHDEBINFO ${DEFINITIONS_STR}
        )
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_name "generic3d_04${platform_suffix}_vclcuda${compiler_suffix}")
        message(STATUS "Building ViennaCL CUDA four variable version: ${target_name}")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA -DFOUR_VARIABLES")
        cuda_add_executable(${target_name} ${cuda_src} ${all_headers}
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
        set(target_name "generic3d_03${platform_suffix}_mkl${compiler_suffix}")
        message(STATUS "Building MKL three variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
        target_link_libraries(${target_name} ${MKL_LIBRARIES})
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_name "generic3d_04${platform_suffix}_mkl${compiler_suffix}")
        message(STATUS "Building MKL four variable version: ${target_name}")
        add_executable(${target_name} ${cpu_src} ${all_headers})
        target_link_libraries(${target_name} ${MKL_LIBRARIES})
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
        set_property(TARGET ${target_name} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
    endif(${FOUR_VARIABLES})
endif(${BUILD_MKL})
