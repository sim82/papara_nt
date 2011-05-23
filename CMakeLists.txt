cmake_minimum_required(VERSION 2.8)

project(papara_nt)
add_subdirectory(ivy_mike)

include_directories( ${include_directories} ivy_mike/src . )

# handle the stupid autolinking shit of visualc++
IF(WIN32)
include_directories( ${include_directories} C:\\src\\boost_1_46_1)
LINK_DIRECTORIES( ${LINK_DIRECTORIES} C:\\src\\boost_1_46_1\\stage\\lib)
set( BOOST_LIBS )
ELSE()
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
set( BOOST_LIBS boost_thread boost_program_options)
ENDIF()


add_executable(papara_nt main.cpp pars_align_vec.cpp)
add_executable(smith_waterman smith_waterman.cpp )
add_executable(stepwise_addition stepwise_addition.cpp pars_align_vec.cpp pairwise_seq_distance.cpp )
add_executable(pw_dist pw_dist.cpp pairwise_seq_distance.cpp )

target_link_libraries(papara_nt ivymike pthread)
#target_link_libraries(stepwise_addition ivymike)
target_link_libraries(stepwise_addition ivymike ${BOOST_LIBS})
target_link_libraries(pw_dist ivymike ${BOOST_LIBS})


