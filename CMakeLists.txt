cmake_minimum_required(VERSION 3.10)
project(zoom)

find_package(PNG REQUIRED)
add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})

find_package(OpenMP)

add_executable(zoomAB zoomAB.cpp cmdLine.h io_png.c io_png.h xmtime.h)
target_link_libraries(zoomAB ${PNG_LIBRARIES})

add_executable(zmorpho zmorpho.cpp cmdLine.h xmtime.h io_png.c io_png.h
  boundaries/flst.cpp boundaries/shape.cpp boundaries/shape.h boundaries/tree.cpp boundaries/tree.h
  microCurv/image.cpp microCurv/image.h microCurv/fill_curve.cpp microCurv/fill_curve.h microCurv/gass.cpp microCurv/gass.cpp microCurv/levelLine.cpp microCurv/levelLine.h microCurv/lltree.cpp microCurv/lltree.h)

target_compile_definitions(zmorpho PUBLIC BOUNDARY)
target_link_libraries(zmorpho ${PNG_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_link_libraries(zoomAB OpenMP::OpenMP_CXX)
  target_link_libraries(zmorpho OpenMP::OpenMP_CXX)
endif()
