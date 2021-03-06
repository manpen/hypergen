# HyperGen
This repository contains a prototype implementation of
"Generating practical random hyperbolic graphs in near-linear time
and with sub-linear memory" [M. Penschuck, SEA 2017].
Please refer to the branch SEA17 to reproduce the results.

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

All experiments are conducted on a system with AVX2. The code builds without it,
but is untested and may yield different results. Therefore, there is a check in
the compile script.

In case your system is not compatible, you may set ARCH=auto and remove the "exit"
after the error message (line ~19). Then simply reexecute ./compile.sh

The implementation itself is a prototype and requires at least one streaming band.
For very small networks, you may therefore get an error message. This is not an
issue of the algorithm itself, but rather out-of-scope of this application. It may
get fixed in the future.

# Experiments
## Building Hyperbolic Embedder
For completeness our experiments optionally include an implementation of the GIRG
model provided by [https://bitbucket.org/HaiZhung/hyperbolic-embedder/]. It is
included as a GIT submodule in related_work. As a few manual modifications are necessary,
it is however not automatically build by the compile script. The following steps allow
building it on a Ubuntu 16.04 system:

Add "#include <ctime>" related_work/hyp_emb/main.cpp.
Replace "-lcblas" by "-lgslcblas" in related_work/hyp_emp/makefile

```
sudo apt-get install libgflags-dev libgoogle-glog-dev libgsl-dev
cd related_work/hyp_emb/
make
```

## Running the Experiments
In order to compare HyperGen, NkGen and NkGenOpt you can run our invokation script.
The analysis scripts require Python3 and a couple of modules. On Ubuntu they may be
installed using the following command:

```
sudo apt-get install python3-scipy python3-pandas python3-matplotlib
```

To run the experiments change into the build directory (not the debug directory).
We assume that you executed the ./compile.sh script (see above) and used gcc.

```
cd build_gcc
python3 ../experiments/run_experiments.py
```

In the current setting, the module will record five repetitions for each data point
and store the results in a "data"-folder in the root of this repository.

WARNING: We limited the number of nodes produced by NkGen/NkGen-Opt using the
MAX_NO_NODES flag in the header of the script to enable experiments on machines
with (comparably) small main memory.

WARNING: Expect a runtime of a few days (and be glad, if it runs faster)

If you want faster results, you may reduce the number of repetitions (iterations variable)
or reduce the max number of edges generated (MAX_NO_EDGES) in the head of the script.

Since the script will skip runs that already have a log file, you can start small,
increase the problem size and rerun the script. In case you terminate a run, please
delete the corresponding partial log file manually.

# Analysing the runs
After running the experiments (or during a run), you can analyse the log files using the
Jupyter Notebook contained in experiments/analyse_logs.ipynb.
