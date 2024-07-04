cmake_minimum_required(VERSION 3.10)
project(zoomPM)

find_package(PNG REQUIRED)
add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})

find_package(OpenMP)

add_executable(zoomPM zoomPM.cpp cmdLine.h io_png.c io_png.h xmtime.h)
target_link_libraries(zoomPM ${PNG_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_link_libraries(zoomPM OpenMP::OpenMP_CXX)
endif()
