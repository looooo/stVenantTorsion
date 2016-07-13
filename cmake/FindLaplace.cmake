
find_path(LAPLACE_INCLUDE_DIR laplace/laplaceKern2D.h)

find_library(LAPLACE_LIB_2D NAMES laplaceKern2D)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LaplaceKernels  DEFAULT_MSG
                                  LAPLACE_INCLUDE_DIR LAPLACE_LIB_2D)