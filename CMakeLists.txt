cmake_minimum_required(VERSION 2.8)
 
PROJECT(nlmClassic)

SET(SRC_LIST main.cpp process.cpp process.h nlm.h nlm.cpp)
add_executable(nlmClassic ${SRC_LIST})

set(VTK_DIR "/home/vinllen/opt/vtk/VTK-build-6.1.0")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(nlmClassic ${VTK_LIBRARIES})
