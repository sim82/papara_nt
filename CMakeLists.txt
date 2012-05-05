cmake_minimum_required(VERSION 2.8)

project(papara_nt)
set( CMAKE_NO_VERBOSE yes )
#add_subdirectory(ivy_mike)
add_subdirectory(ublasJama-1.0.2.3)


IF( NOT IVYMIKE_INCLUDE )
SET( IVYMIKE_INCLUDE ./ivy_mike )
ENDIF()

include_directories( ${include_directories} ${IVYMIKE_INCLUDE} . ublasJama-1.0.2.3 )

# handle the stupid autolinking shit of visualc++
IF(WIN32)
include_directories( ${include_directories} E:\\src\\boost_1_48_0)
LINK_DIRECTORIES( ${LINK_DIRECTORIES} E:\\src\\boost_1_48_0\\stage\\lib)
set( BOOST_LIBS )
set( SYSDEP_LIBS )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_SCL_SECURE_NO_WARNINGS") # as long as there is no support for std::array, these warnings are plain stupid!

file(GLOB ALL_HEADERS *.h)

ELSE()
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -march=native")
set( BOOST_LIBS boost_thread boost_program_options)
set(SYSDEP_LIBS pthread)
#LINK_DIRECTORIES( ${LINK_DIRECTORIES} /usr/lib64/atlas-sse2 )

set( ALL_HEADERS )
ENDIF()


ADD_LIBRARY( papara_core STATIC papara.cpp pvec.cpp pars_align_seq.cpp pars_align_gapp_seq.cpp parsimony.cpp sequence_model.cpp )

# add_executable(papara_nt main.cpp pvec.cpp pars_align_seq.cpp pars_align_gapp_seq.cpp parsimony.cpp ${ALL_HEADERS})
add_executable(papara papara2_main.cpp  ${ALL_HEADERS})
add_executable(fasta_random_sample2 fasta_random_sample2.cpp  ${ALL_HEADERS})
add_executable(fasta_to_phy fasta_to_phy.cpp  ${ALL_HEADERS})
# add_executable(pw_dist pw_dist.cpp pairwise_seq_distance.cpp )

# target_link_libraries(papara_nt ivymike ublas_jama ${SYSDEP_LIBS} )
target_link_libraries(papara papara_core ivymike ublas_jama ${SYSDEP_LIBS} )
# target_link_libraries(pw_dist ivymike ${BOOST_LIBS} ${SYSDEP_LIBS})


# more experimental stuff that I don't care to force through microsoft's interpretation of c++...
# also that stuff externally only depends on boost headers, which removes most windows related P.I.T.A.
# if( NOT WIN32 ) 
# add_executable(tree_similarity tree_similarity.cpp call_main2.cpp)
# add_executable(propara propara.cpp pvec.cpp raxml_interface.cpp tree_similarity.cpp parsimony.cpp )
# add_executable(smith_waterman smith_waterman.cpp )
# add_executable(stepwise_addition stepwise_addition.cpp pars_align_seq.cpp pairwise_seq_distance.cpp tree_similarity.cpp parsimony.cpp )
# add_executable(stepwise_addition_gapp stepwise_addition_gapp.cpp pvec.cpp pars_align_seq.cpp pairwise_seq_distance.cpp tree_similarity.cpp parsimony.cpp raxml_interface.cpp )
# add_executable(stepwise_addition_pro stepwise_addition_pro.cpp pvec.cpp pairwise_seq_distance.cpp tree_similarity.cpp parsimony.cpp raxml_interface.cpp sequence_model.cpp )
# add_executable(epa_extract_qs_covered epa_extract_qs_covered.cpp )
# 
# 
# add_executable(tbb_test tbb.cpp )
# add_executable(dtw dtw.cpp )
# add_executable(inherit_test inherit_test.cpp )
# add_executable(testbench testbench.cpp pars_align_seq.cpp )
# # add_executable(atlas_test atlas.cpp )
# 
# add_executable(test_bitset test_bitset.cpp)
# #add_executable(endian_benchmark endian_benchmark.cpp)
# 
# target_link_libraries(tree_similarity ivymike)
# target_link_libraries(propara ivymike ${SYSDEP_LIBS} PocoFoundation ublas_jama ${BOOST_LIBS})
# #target_link_libraries(stepwise_addition ivymike)
# target_link_libraries(stepwise_addition ivymike ${BOOST_LIBS} ${SYSDEP_LIBS})
# target_link_libraries(stepwise_addition_gapp ivymike ublas_jama ${BOOST_LIBS} ${SYSDEP_LIBS} PocoFoundation)
# target_link_libraries(stepwise_addition_pro ivymike ${BOOST_LIBS} ${SYSDEP_LIBS} PocoFoundation)
# 
# #target_link_libraries( atlas_test ivymike ${BOOST_LIBS} atlas cblas )
# target_link_libraries(tbb_test ublas_jama ivymike )
# target_link_libraries(testbench ${BOOST_LIBS} )
# target_link_libraries(epa_extract_qs_covered ivymike)
# 
# target_link_libraries(dtw ivymike)
# target_link_libraries(inherit_test ivymike)
# 
# target_link_libraries(test_bitset ivymike)
# #target_link_libraries(endian_benchmark ivymike)
# 
# endif()
