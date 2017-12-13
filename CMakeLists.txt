	cmake_minimum_required (VERSION 3.1)
project (Generator)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_BUILD_TYPE DEBUG)
#set(CMAKE_CXX_STANDARD 11) # C++11...
#set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
#set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

set (WARN_FLAGS              "-Wall -Wextra -Wno-unknown-pragmas")
set (CMAKE_C_FLAGS         "${CMAKE_C_FLAGS} ${WARN_FLAGS}")

set (CMAKE_C_FLAGS_DEBUG   "-O0 -g -lm -ggdb3 -DDEBUG")

SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lm ")

include_directories ("${PROJECT_SOURCE_DIR}")
add_subdirectory (RAxML) 

 
# add the executable
add_executable(Generator generator.c)
target_link_libraries (Generator RAxML)
