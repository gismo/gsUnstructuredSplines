######################################################################
## CMakeLists.txt --- gsUnstructuredSplines/unittests
## This file is part of the G+Smo library.
##
######################################################################

# add example files
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR} unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")