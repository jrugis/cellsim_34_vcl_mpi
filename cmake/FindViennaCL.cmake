
find_path(
    VIENNACL_INCLUDE_DIR
    NAMES
    viennacl/vector.hpp
    PATHS
    /projects/nesi00119/code/CS_dev/viennacl/ViennaCL-1.7.1/
    /projects/nesi00119/code/viennacl/ViennaCL-1.7.1
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ViennaCL DEFAULT_MSG VIENNACL_INCLUDE_DIR)

set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR})
mark_as_advanced(VIENNACL_INCLUDE_DIRS VIENNACL_INCLUDE_DIR)
