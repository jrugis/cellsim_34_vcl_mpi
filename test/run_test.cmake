
#
# CMake script to run a single test case
#
message(STATUS "Running test: ${TEST_NAME}")
message(STATUS "  Test run directory: ${TEST_RUN_DIR}")
message(STATUS "  Test src directory: ${TEST_SRC_DIR}")
message(STATUS "  CellSim binary: ${TEST_BINARY}")
message(STATUS "  Skip comparison: ${SKIP_COMPARE}")
message(STATUS "  Python reduce: ${TEST_PY_REDUCE}")

#
# make the test directory
#
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_RUN_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${TEST_RUN_DIR})

#
# copy mesh file
#
file(COPY ${TEST_SRC_DIR}/cell01m_HARMONIC_100p.msh DESTINATION ${TEST_RUN_DIR})
file(RENAME ${TEST_RUN_DIR}/cell01m_HARMONIC_100p.msh ${TEST_RUN_DIR}/cs.msh)

#
# copy parameter file
#
string(REGEX MATCH "03" IS_3VAR ${TEST_BINARY})
if (IS_3VAR)
    file(COPY ${TEST_SRC_DIR}/generic3d_03-cs.dat DESTINATION ${TEST_RUN_DIR})
    file(RENAME ${TEST_RUN_DIR}/generic3d_03-cs.dat ${TEST_RUN_DIR}/cs.dat)
else()
    file(COPY ${TEST_SRC_DIR}/generic3d_04-cs.dat DESTINATION ${TEST_RUN_DIR})
    file(RENAME ${TEST_RUN_DIR}/generic3d_04-cs.dat ${TEST_RUN_DIR}/cs.dat)
endif()

#
# If testing python reduce, disable C++ reduce
#
if (TEST_PY_REDUCE)
    # change variable in parameter file
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
        sed -e "s/0.05 120 1/0.05 120 0/g" cs.dat
        RESULT_VARIABLE STATUS
        OUTPUT_FILE ${TEST_RUN_DIR}/cs.dat.new
    )
    if (STATUS)
        message(FATAL_ERROR "sed command failed: ${STATUS}")
    endif ()
    file(RENAME ${TEST_RUN_DIR}/cs.dat.new ${TEST_RUN_DIR}/cs.dat)
endif()

#
# run cellsim
#
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR} ${TEST_BINARY}
    RESULT_VARIABLE STATUS
)
if (STATUS)
    message(FATAL_ERROR "Error running cellsim: '${STATUS}'")
endif (STATUS)

#
# python is required for the comparison
#
find_package(PythonInterp)  # we require python for this
if (NOT PYTHONINTERP_FOUND)
    message(WARNING "Skipping comparison as python not available")
    set(SKIP_COMPARE ON)
endif()

#
# If testing python reduce, run the reduce now
#
if (TEST_PY_REDUCE)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
        ${PYTHON_EXECUTABLE} ${TEST_SRC_DIR}/cs_reduce_min-max.py
        RESULT_VARIABLE STATUS
    )
    if (STATUS)
        message(FATAL_ERROR "Error running reduce: '${STATUS}'")
    endif (STATUS)
endif()

#
# compare output files
#
if (NOT SKIP_COMPARE)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
        ${PYTHON_EXECUTABLE} ${TEST_SRC_DIR}/cs_compare_peaks.py cR.bin ${TEST_SRC_DIR}/generic3d-cR.bin
        RESULT_VARIABLE STATUS
    )
    if (STATUS)
        message(FATAL_ERROR "Error during comparison: '${STATUS}'")
    endif (STATUS)
endif ()
