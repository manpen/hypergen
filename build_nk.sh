#!/usr/bin/env bash
cd libs/NetworKit

rm *.a

echo "#define HYPERBOLIC_NO_EDGES" > networkit/cpp/generators/HyperbolicBuildConfig.h
scons --optimize=Opt --target=Lib -j8
mv libNetworKit-Core-Opt.a libNetworKitNoEdges.a

echo "" > networkit/cpp/generators/HyperbolicBuildConfig.h
scons --optimize=Opt --target=Lib -j8
rm libNetworKit.a
mv libNetworKit-Core-Opt.a libNetworKit.a