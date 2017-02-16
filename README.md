# Scripts to reproduce Experiments will be added before 22. Feb!

# Building all binaries
In order to build HyperGen you need a recent GCC compiler with OpenMP support,
cmake and scons (for NetworKit). On Ubuntu 16.04 you can install the dependencies
using the following command:

```
sudo apt-get install g++-5 git scons libgomp1 cmake
```

In order to retrieve the code, fetch libVc and to compile NetworKit, libVc and HyperGen
you can use the following commands:

```
sudo apt-get install g++-5 git scons libgomp1 cmake
git clone https://github.com/manpen/hypergen.git
cd hypergen
./compile.sh
cd build_gcc
./main_hyper -n 1000000 -d 1000
```

All experiments are conducted on a system with AVX2. The code builds without,
but is untested and may yield different results. Therefore, there is a check in
the compile script.

In case your system is not compatible, you may set ARCH=auto and remove the "exit"
after the error message (line ~19). Then simply reexecute ./compile.sh

The implementation itself is a prototype and requires at least one streaming band.
For very small networks, you may therefore get an error message. This is not an
issue of the algorithm itself, but rather out-of-scope of this application. It may
get fixed in the future.
