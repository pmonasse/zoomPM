cmake_minimum_required(VERSION 2.6)
project(zoom)

find_package(PNG REQUIRED)
add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})

add_executable(zoomAB main_zoomAB.cpp cmdLine.h io_png.c io_png.h)
target_link_libraries(zoomAB ${PNG_LIBRARIES})

add_executable(zmorpho zmorpho.cpp cmdLine.h xmtime.h io_png.c io_png.h
  boundaries/flst.cpp boundaries/shape.cpp boundaries/shape.h boundaries/tree.cpp boundaries/tree.h
  microCurv/image.cpp microCurv/image.h microCurv/fill_curve.cpp microCurv/fill_curve.h microCurv/gass.cpp microCurv/gass.cpp microCurv/levelLine.cpp microCurv/levelLine.h microCurv/lltree.cpp microCurv/lltree.h)

target_compile_definitions(zmorpho PUBLIC BOUNDARY)
target_link_libraries(zmorpho ${PNG_LIBRARIES})
