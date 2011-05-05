project(papara_nt)
add_subdirectory(ivy_mike)
include_directories( ivy_mike/src . )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
add_executable(papara_nt main.cpp)

target_link_libraries(papara_nt ivymike)

