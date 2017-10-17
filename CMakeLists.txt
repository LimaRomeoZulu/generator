cmake_minimum_required (VERSION 2.6)
project (Generator)

include_directories ("${PROJECT_SOURCE_DIR}")
add_subdirectory (standard-RaxML) 

set(CMAKE_BUILD_TYPE DEBUG)
 
# add the executable
add_executable(Generator Generator.c)
target_link_libraries (Generator standard-RaxML)
