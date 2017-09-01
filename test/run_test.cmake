
#
# CMake script to run a single test case
#
message(STATUS "Running test: ${TEST_NAME}")
message(STATUS "  Test run directory: ${TEST_RUN_DIR}")
message(STATUS "  Test src directory: ${TEST_SRC_DIR}")
message(STATUS "  CellSim binary: ${TEST_BINARY}")
message(STATUS "  Skip comparison: ${SKIP_COMPARE}")
message(STATUS "  Python reduce: ${TEST_PY_REDUCE}")
message(STATUS "  Python found: ${PYTHONINTERP_FOUND}")
message(STATUS "  Python executable: ${PYTHON_EXECUTABLE}")
message(STATUS "  Enable gprof: ${ENABLE_GPROF}")
message(STATUS "  MPIEXEC: ${MPIEXEC}")
message(STATUS "  MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
message(STATUS "  MPIEXEC_PREFLAGS: ${MPIEXEC_PREFLAGS}")
message(STATUS "  MPIEXEC_POSTFLAGS: ${MPIEXEC_POSTFLAGS}")
if (ENABLE_GPROF)
    message(STATUS "    gprof: ${GPROF_PROGRAM}")
    message(STATUS "    gprof2dot: ${GPROF2DOT_PROGRAM}")
    message(STATUS "    dot: ${DOT_PROGRAM}")
endif()

#
# make the test directory
#
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_RUN_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${TEST_RUN_DIR})

#
# copy mesh files
#
foreach (index RANGE 1 7)
    file(COPY ${TEST_SRC_DIR}/out_N4_p3-p2-p4-${index}tet.msh DESTINATION ${TEST_RUN_DIR})
    file(RENAME ${TEST_RUN_DIR}/out_N4_p3-p2-p4-${index}tet.msh ${TEST_RUN_DIR}/cs${index}.msh)
endforeach()

#
# copy parameter file
#
foreach (index RANGE 1 7)
    string(REGEX MATCH "03" IS_3VAR ${TEST_BINARY})
    if (IS_3VAR)
        file(COPY ${TEST_SRC_DIR}/generic3d_03-cs.dat DESTINATION ${TEST_RUN_DIR})
        file(RENAME ${TEST_RUN_DIR}/generic3d_03-cs.dat ${TEST_RUN_DIR}/cs${index}.dat)
    else()
        file(COPY ${TEST_SRC_DIR}/generic3d_04-cs.dat DESTINATION ${TEST_RUN_DIR})
        file(RENAME ${TEST_RUN_DIR}/generic3d_04-cs.dat ${TEST_RUN_DIR}/cs${index}.dat)
    endif()

    #
    # If testing python reduce, disable C++ reduce
    #
    if (TEST_PY_REDUCE)
        # change variable in parameter file
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
            sed -e "s/0.05 120 1/0.05 120 0/g" cs${index}.dat
            RESULT_VARIABLE STATUS
            OUTPUT_FILE ${TEST_RUN_DIR}/cs${index}.dat.new
        )
        if (STATUS)
            message(FATAL_ERROR "sed command failed: ${STATUS}")
        endif ()
        file(RENAME ${TEST_RUN_DIR}/cs${index}.dat.new ${TEST_RUN_DIR}/cs${index}.dat)
    endif()
endforeach()

#
# run cellsim
#
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 7
            ${MPIEXEC_PREFLAGS} ${TEST_BINARY} ${MPIEXEC_POSTFLAGS}
    RESULT_VARIABLE STATUS
)
if (STATUS)
    message(FATAL_ERROR "Error running cellsim: '${STATUS}'")
endif (STATUS)

#
# Handle the profiling result
#
if (ENABLE_GPROF AND GPROF_PROGRAM)
    message(STATUS "Running gprof")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
        ${GPROF_PROGRAM} ${TEST_BINARY}
        OUTPUT_FILE ${TEST_RUN_DIR}/profile_${TEST_NAME}.txt
        RESULT_VARIABLE STATUS
    )
    message(STATUS "Status of gprof: ${STATUS}")

    if (NOT STATUS)
        # also try to locate gprof2dot and dot
        if (GPROF2DOT_PROGRAM AND DOT_PROGRAM)
            message(STATUS "Running gprof2dot")
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
                ${GPROF2DOT_PROGRAM} --colour-nodes-by-selftime --strip profile_${TEST_NAME}.txt
                OUTPUT_FILE ${TEST_RUN_DIR}/gprof2dot_output
                RESULT_VARIABLE STATUS
            )
            message(STATUS "Status of gprof2dot: ${STATUS}")

            if (NOT STATUS)
                message(STATUS "Running dot")
                execute_process(
                    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
                    ${DOT_PROGRAM} -Tpng -o ${TEST_RUN_DIR}/profile_${TEST_NAME}.png
                        ${TEST_RUN_DIR}/gprof2dot_output
                    RESULT_VARIABLE STATUS
                )
                message(STATUS "Status of dot: ${STATUS}")

                if (NOT STATUS)
                    execute_process(
                        COMMAND ${CMAKE_COMMAND} -E copy
                        ${TEST_RUN_DIR}/profile_${TEST_NAME}.png
                        ../profile_${TEST_NAME}.png
                    )
                endif()
            endif()
        endif()
    endif()
elseif (ENABLE_GPROF AND NOT GPROF_PROGRAM)
    message(WARNING "Skipping gprof as could not locate executable")
endif ()

#
# skip the comparison if python wasn't found
#
if (NOT PYTHONINTERP_FOUND)
    message(WARNING "Python not available: skipping reduction/comparison")
    set(SKIP_COMPARE ON)
    set(TEST_PY_REDUCE OFF)
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
        ${PYTHON_EXECUTABLE} ${TEST_SRC_DIR}/cs_compare_peaks.py cR1.bin ${TEST_SRC_DIR}/generic3d-cR.bin
        RESULT_VARIABLE STATUS
    )
    if (STATUS)
        message(FATAL_ERROR "Error during comparison: '${STATUS}'")
    endif (STATUS)
endif ()
