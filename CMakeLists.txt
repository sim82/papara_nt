project(papara_nt)
add_subdirectory(ivy_mike)
include_directories( ivy_mike/src . C:\\src\\boost_1_46_1)
LINK_DIRECTORIES( ${LINK_DIRECTORIES} C:\\src\\boost_1_46_1\\stage\\lib)

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
add_executable(papara_nt main.cpp pars_align_vec.cpp)
add_executable(smith_waterman smith_waterman.cpp )
add_executable(stepwise_addition stepwise_addition.cpp pars_align_vec.cpp pairwise_seq_distance.cpp )


target_link_libraries(papara_nt ivymike)
#target_link_libraries(stepwise_addition ivymike boost_thread)
target_link_libraries(stepwise_addition ivymike)