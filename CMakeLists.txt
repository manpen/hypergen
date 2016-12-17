cmake_minimum_required(VERSION 3.0)
project(mhyperbolic)

include_directories(libs/NetworKit/include/)
link_directories(libs/NetworKit)

#set(Vc_DIR "libs/Vc/")
find_package(Vc)
add_compile_options(${Vc_ALL_FLAGS})
#include_directories(libs/Vc/include/)
#link_directories(libs/Vc/build/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp -g -march=native -ffast-math")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fsanitize=address")


add_executable(mhyperbolic main.cpp)

add_executable(benchmark_networkit benchmark_networkit.cpp)
target_link_libraries(benchmark_networkit NetworKit)

add_executable(main_hyper 
    main_hyper.cpp
    include/Geometry.cpp
    include/PointGenerator.cpp
    include/RandomHelper.cpp
    include/BandSegment.cpp
    include/Segment.cpp
    include/Generator.cpp
    include/Histogram.cpp
)

target_link_libraries(main_hyper NetworKit Vc)


#target_link_libraries(main_hyper NetworKit)
#install(TARGETS mhyperbolic RUNTIME DESTINATION bin)