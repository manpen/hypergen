cmake_minimum_required(VERSION 3.0)
project(mhyperbolic)

include_directories(libs/NetworKit/include/)
link_directories(libs/NetworKit)

#set(Vc_DIR "libs/Vc/")
find_package(Vc REQUIRED)
add_compile_options(${Vc_ALL_FLAGS})
message("VC Dir: ${Vc_INCLUDE_DIR}")
message("VC Flags: ${Vc_ALL_FLAGS}")
include_directories(${Vc_INCLUDE_DIR})
link_directories("${Vc_INCLUDE_DIR}/../lib/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g -O0")

if (${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -ffast-math -DNDEBUG")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # enable -xHOST -O3 -ipo -no-prec-div -static -fp-model fast=2
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 -march=core-avx2 -fno-alias ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -use-gen -prof-dir=./profiler")
endif()


add_executable(benchmark_networkit benchmark_networkit.cpp)
target_link_libraries(benchmark_networkit NetworKitNoEdges)

set(HYPER_FILES include/Geometry.cpp
    include/PointGenerator.cpp
    include/RandomHelper.cpp
    include/BandSegment.cpp
    include/Segment.cpp
    include/Generator.cpp
    include/Histogram.cpp
    include/Configuration.cpp
    include/ActiveManager.cpp
)

# Main Binary
add_executable(main_hyper 
    main_hyper.cpp
    ${HYPER_FILES}
)
target_link_libraries(main_hyper NetworKit Vc)

# Binary to measure time to compute distances
add_executable(main_hyper_skipdist
        main_hyper.cpp
        ${HYPER_FILES}
)
target_link_libraries(main_hyper_skipdist NetworKit Vc)
set_target_properties(main_hyper_skipdist PROPERTIES COMPILE_FLAGS "-DSKIP_DIST_COMP=1")

