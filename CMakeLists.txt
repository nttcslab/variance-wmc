cmake_minimum_required(VERSION 3.8)
project(variance)

set(CMAKE_CXX_FLAGS "-Wno-unused-result")

add_library(variancelib STATIC src/vtree.cpp src/sdd.cpp src/variancecalc.cpp)
target_compile_features(variancelib PUBLIC cxx_std_11)

add_executable(main src/main.cpp)
target_link_libraries(main variancelib)

add_executable(enc2lmap src/enc2lmap.cpp)
target_compile_features(enc2lmap PUBLIC cxx_std_11)