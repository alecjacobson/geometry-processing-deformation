cmake_minimum_required(VERSION 2.6)
project(deformation)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED)

if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
  add_compile_options(-mavx2)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(FATAL_ERROR "Intel not supported")
elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
  add_compile_options(/arch:AVX2)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_definitions(${CMAKE_CXX_FLAGS} "-mfpmath=sse")
add_definitions(${CMAKE_CXX_FLAGS} "-march=native")
add_definitions(${CMAKE_CXX_FLAGS} "-g")
add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

# Add your project files
include_directories("include/")
if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES src/*.cpp)
endif()

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::opengl igl::opengl_glfw )

add_executable(deformation "main.cpp")
target_link_libraries(deformation core igl::core igl::opengl igl::opengl_glfw )
