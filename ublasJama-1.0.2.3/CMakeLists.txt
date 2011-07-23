IF(WIN32)
include_directories( ${include_directories} E:\\src\\boost_1_46_1\\boost_1_46_1)

ELSE()
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -march=native")

ENDIF()

add_library( ublas_jama CholeskyDecomposition.cpp  EigenvalueDecomposition.cpp  LUDecomposition.cpp  QRDecomposition.cpp )
