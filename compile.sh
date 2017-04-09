#!/bin/bash
export CXX="g++"
export CC="gcc"
ARCH=broadwell

git submodule update --init

CXX_VERSION=$($CXX -dumpversion | cut -f1 -d.)
if [ "$CXX_VERSION" -lt "5" ]; then
   echo "Only tested with g++5 (or higher); may work anyways"
   exit
fi

if ! grep -q "avx2" /proc/cpuinfo; then
   echo "It seems that your CPU does not support AVX2"
   echo "Please adopt the ARCH parameter in this file and comment"
   echo "out this check; you may experience a reduced performance"
   echo "see libs/Vc/cmake/OptimizeForArchitecture.cmake"
   exit
fi

CORES=$(cat /proc/cpuinfo | grep "processor" | wc -l)
echo "Use $CORES cores for building"
   
echo "Build libvc"
BDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VCDIR="$BDIR/libs/Vc/"
VCDIRInst="$VCDIR/install_$CC"
export Vc_DIR="$VCDIRInst/lib/cmake/Vc/"

# BUILD VC
cd $VCDIR
    mkdir build_$CC
    cd build_$CC
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$VCDIRInst -DTARGET_ARCHITECTURE=$ARCH -DBUILD_TESTING=OFF ..
    make -j $CORES
    make install


echo "Build NetworKIT"
cd $BDIR/libs/NetworKit
    perl -pi -e "s/cpp\s=.*/cpp = $CXX/g" build.conf
    rm *.a
 
    echo "#define HYPERBOLIC_NO_EDGES" > networkit/cpp/generators/HyperbolicBuildConfig.h
    scons --optimize=Opt --target=Lib -j$CORES
    mv libNetworKit-Core-Opt.a libNetworKitNoEdges.a

    echo "#define HYPERBOLIC_SKIP_DIST" >> networkit/cpp/generators/HyperbolicBuildConfig.h
    scons --optimize=Opt --target=Lib -j$CORES
    mv libNetworKit-Core-Opt.a libNetworKitSkipDist.a

    echo "" > networkit/cpp/generators/HyperbolicBuildConfig.h
    scons --optimize=Opt --target=Lib -j$CORES
    rm libNetworKit.a
    mv libNetworKit-Core-Opt.a libNetworKit.a
    
cd $BDIR

echo "Build HyperGen"
    mkdir debug_$CC
    cd debug_$CC
    cmake -DCMAKE_BUILD_TYPE=Debug -DTARGET_ARCHITECTURE=$ARCH ..
    make -j $CORES
    cd ..

    mkdir build_$CC
    cd build_$CC
    cmake -DCMAKE_BUILD_TYPE=Release -DTARGET_ARCHITECTURE=$ARCH ..
    make -j $CORES
