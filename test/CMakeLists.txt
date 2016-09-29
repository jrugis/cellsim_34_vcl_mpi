
#
# macro for adding a test
#
macro(add_cellsim_test TEST_NAME)
    message(STATUS "Adding test: '${TEST_NAME}'")
    set(TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DTEST_NAME=${TEST_NAME}
            -DTEST_RUN_DIR=${TEST_RUN_DIR}
            -DTEST_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DTEST_BINARY=$<TARGET_FILE:${TEST_NAME}>
            -DSKIP_COMPARE=${TEST_SKIP_COMPARE}
            -DTEST_PY_REDUCE=OFF
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
    )
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 3000)
    
    if (TEST_PYTHON_REDUCE)
        message(STATUS "Adding test: '${TEST_NAME}_pyreduce'")
        set(TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}_pyreduce)
        add_test(
            NAME ${TEST_NAME}_pyreduce
            COMMAND ${CMAKE_COMMAND}
                -DTEST_NAME=${TEST_NAME}_pyreduce
                -DTEST_RUN_DIR=${TEST_RUN_DIR}
                -DTEST_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DTEST_BINARY=$<TARGET_FILE:${TEST_NAME}>
                -DSKIP_COMPARE=${TEST_SKIP_COMPARE}
                -DTEST_PY_REDUCE=ON
                -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
        )
        set_tests_properties(${TEST_NAME}_pyreduce PROPERTIES TIMEOUT 3000)
    endif ()
endmacro(add_cellsim_test)

# serial version (ViennaCL)
if(${BUILD_SERIAL})
    # three variable
    if(${THREE_VARIABLES})
        add_cellsim_test(${target_vcl_serial_03})
    endif(${THREE_VARIABLES})
    # four variable
    if(${FOUR_VARIABLES})
        add_cellsim_test(${target_vcl_serial_04})
    endif(${FOUR_VARIABLES})
endif(${BUILD_SERIAL})
