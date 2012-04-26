IF(WIN32)
include_directories( ${include_directories} E:\\src\\boost_1_48_0)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_SCL_SECURE_NO_WARNINGS") # as long as there is no support for std::array, these warnings are plain stupid!
ELSE()
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -march=native")

ENDIF()

add_library( ublas_jama CholeskyDecomposition.cpp  EigenvalueDecomposition.cpp  LUDecomposition.cpp  QRDecomposition.cpp )
