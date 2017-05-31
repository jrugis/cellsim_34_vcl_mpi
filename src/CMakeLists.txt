
# source
file(GLOB cpu_src ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB all_headers ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB cuda_src ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)

# serial version (ViennaCL)
if(${BUILD_SERIAL})
    # three variable
    if(${THREE_VARIABLES})
        set(target_vcl_serial_03 "generic3d_03_${platform_suffix}_vcl-serial-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_serial_03 "generic3d_03_${platform_suffix}_vcl-serial-${compiler_suffix}")
        message(STATUS "Building ViennaCL serial three variable version: ${target_vcl_serial_03}")
        add_executable(${target_vcl_serial_03} ${cpu_src} ${all_headers})
        install(TARGETS ${target_vcl_serial_03} DESTINATION bin)
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_vcl_serial_04 "generic3d_04_${platform_suffix}_vcl-serial-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_serial_04 "generic3d_04_${platform_suffix}_vcl-serial-${compiler_suffix}")
        message(STATUS "Building ViennaCL serial four variable version: ${target_vcl_serial_04}")
        add_executable(${target_vcl_serial_04} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_vcl_serial_04} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
        install(TARGETS ${target_vcl_serial_04} DESTINATION bin)
    endif(${FOUR_VARIABLES})
endif(${BUILD_SERIAL})

# OpenMP version (ViennaCL)
if(${BUILD_OPENMP})
    # three variable
    if(${THREE_VARIABLES})
        set(target_vcl_omp_03 "generic3d_03_${platform_suffix}_vcl-openmp-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_omp_03 "generic3d_03_${platform_suffix}_vcl-openmp-${compiler_suffix}")
        message(STATUS "Building ViennaCL OpenMP three variable version: ${target_vcl_omp_03}")
        add_executable(${target_vcl_omp_03} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_vcl_omp_03} APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_vcl_omp_03} APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_vcl_omp_03} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
        install(TARGETS ${target_vcl_omp_03} DESTINATION bin)
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_vcl_omp_04 "generic3d_04_${platform_suffix}_vcl-openmp-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_omp_04 "generic3d_04_${platform_suffix}_vcl-openmp-${compiler_suffix}")
        message(STATUS "Building ViennaCL OpenMP four variable version: ${target_vcl_omp_04}")
        add_executable(${target_vcl_omp_04} ${cpu_src} ${all_headers})
        set_property(TARGET ${target_vcl_omp_04} APPEND_STRING PROPERTY COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_vcl_omp_04} APPEND_STRING PROPERTY LINK_FLAGS "${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
        set_property(TARGET ${target_vcl_omp_04} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "VIENNACL_WITH_OPENMP")
        set_property(TARGET ${target_vcl_omp_04} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
        install(TARGETS ${target_vcl_omp_04} DESTINATION bin)
    endif(${FOUR_VARIABLES})
endif(${BUILD_OPENMP})

# CUDA version (ViennaCL)
if(${BUILD_CUDA})
    # three variable
    if(${THREE_VARIABLES})
        set(target_vcl_cuda_03 "generic3d_03_${platform_suffix}_vcl-cuda-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_cuda_03 "generic3d_03_${platform_suffix}_vcl-cuda-${compiler_suffix}")
        message(STATUS "Building ViennaCL CUDA three variable version: ${target_vcl_cuda_03}")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA")
        cuda_add_executable(${target_vcl_cuda_03} ${cuda_src} ${all_headers}
            OPTIONS
                DEBUG ${DEFINITIONS_STR}
                RELEASE ${DEFINITIONS_STR}
                MINSIZEREL ${DEFINITIONS_STR}
                RELWITHDEBINFO ${DEFINITIONS_STR}
        )
        install(TARGETS ${target_vcl_cuda_03} DESTINATION bin)
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_vcl_cuda_04 "generic3d_04_${platform_suffix}_vcl-cuda-${compiler_suffix}" PARENT_SCOPE)
        set(target_vcl_cuda_04 "generic3d_04_${platform_suffix}_vcl-cuda-${compiler_suffix}")
        message(STATUS "Building ViennaCL CUDA four variable version: ${target_vcl_cuda_04}")
        set(DEFINITIONS_STR "-DVIENNACL_WITH_CUDA -DFOUR_VARIABLES")
        cuda_add_executable(${target_vcl_cuda_04} ${cuda_src} ${all_headers}
            OPTIONS
                DEBUG ${DEFINITIONS_STR}
                RELEASE ${DEFINITIONS_STR}
                MINSIZEREL ${DEFINITIONS_STR}
                RELWITHDEBINFO ${DEFINITIONS_STR}
        )
        install(TARGETS ${target_vcl_cuda_04} DESTINATION bin)
    endif(${FOUR_VARIABLES})
endif(${BUILD_CUDA})

# MKL version
if(${BUILD_MKL})
    # three variable
    if(${THREE_VARIABLES})
        set(target_mkl_03 "generic3d_03_${platform_suffix}_mkl-${compiler_suffix}" PARENT_SCOPE)
        set(target_mkl_03 "generic3d_03_${platform_suffix}_mkl-${compiler_suffix}")
        message(STATUS "Building MKL three variable version: ${target_mkl_03}")
        add_executable(${target_mkl_03} ${cpu_src} ${all_headers})
        target_link_libraries(${target_mkl_03} ${MKL_LIBRARIES})
        set_property(TARGET ${target_mkl_03} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
        install(TARGETS ${target_mkl_03} DESTINATION bin)
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        set(target_mkl_04 "generic3d_04_${platform_suffix}_mkl-${compiler_suffix}" PARENT_SCOPE)
        set(target_mkl_04 "generic3d_04_${platform_suffix}_mkl-${compiler_suffix}")
        message(STATUS "Building MKL four variable version: ${target_mkl_04}")
        add_executable(${target_mkl_04} ${cpu_src} ${all_headers})
        target_link_libraries(${target_mkl_04} ${MKL_LIBRARIES})
        set_property(TARGET ${target_mkl_04} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "MKL_SOLVER")
        set_property(TARGET ${target_mkl_04} APPEND_STRING PROPERTY COMPILE_DEFINITIONS "FOUR_VARIABLES")
        install(TARGETS ${target_mkl_04} DESTINATION bin)
    endif(${FOUR_VARIABLES})
endif(${BUILD_MKL})
