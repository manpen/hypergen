#!/usr/bin/env python3

# @file
# @brief Perform Experiments
#
# This script invokes the generators for all settings investigated.
# The executation will require a machine with at least 64 GB RAM
# (decrease MAX_NO_EDGES if your machine has less memory). Expect
# a runtime of several days.
#
# @author Manuel Penschuck
# @copyright
# Copyright (C) 2017 Manuel Penschuck
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# @copyright
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# @copyright
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import subprocess, os, sys
from math import exp, pi
import multiprocessing
import socket
import itertools
import random

# decrease these constants to speed-up simulation
MAX_NO_EDGES = int(1e12)
MAX_NO_NODES = int(1e8)
iterations = 5

def invokeGenerator(n, avgDeg, expv, outf, seed=-1, segments=-1, worker=-1, bandSpacing=2):
    if not os.path.isfile("./main_hyper"):
        print("Did not find binary ./main_hyper -- please run this script")
        print("from within the build directory.")
        sys.exit(-1)

    if seed < 0:
        seed = random.randint(0, 0xfffffff)

    args = ["./main_hyper", "-n", n, "-d", avgDeg, "-e", expv, "-s", seed, "-b", bandSpacing]

    if (segments > 0):
        args += ["-w", segments]

    my_env = os.environ.copy()
    if "OMP_NUM_THREADS" in my_env:
        del(my_env["OMP_NUM_THREADS"])
    if (worker > 0):
        my_env["OMP_NUM_THREADS"] = str(worker)

    args = [str(x) for x in args]
    print("Execute " + " ".join(args))
    args = ["/usr/bin/time", "-av"] + args
    subprocess.call(args, stdout=outf, stderr=outf, env=my_env)

def invokeNkGenerator(n, avgDeg, expv, mode, outf, seed=-1, segments=-1, worker=-1):
    if not os.path.isfile("./benchmark_networkit"):
        print("Did not find binary ./benchmark_networkit")
        print("Please run this script from within the build directory.")
        sys.exit(-1)

    if seed < 0:
        seed = random.randint(0, 0xfffffff)

    args = ["./benchmark_networkit", "-n", n, "-d", avgDeg, "-e", expv, "-s", seed, "-g", mode]

    if (segments > 0):
        args += ["-w", segments]

    my_env = os.environ.copy()
    if "OMP_NUM_THREADS" in my_env:
        del(my_env["OMP_NUM_THREADS"])
    if (worker > 0):
        my_env["OMP_NUM_THREADS"] = str(worker)

    args = [str(x) for x in args]
    print("Execute " + " ".join(args))
    args = ["/usr/bin/time", "-av"] + args
    subprocess.call(args, stdout=outf, stderr=outf, env=my_env)

def exp10Series(n0, decs, pointsPerDec=3, forceInt = True):
    series = [n0 * 10 ** (x / pointsPerDec) for x in range(decs*pointsPerDec + 1)]
    if forceInt:
        series = map(int, series)
    return series

def runSeries(prefix, iterations, confs, exps, algos, cpus, datadir="../data/"):
    for it in range(1, iterations+1):
        for algo in algos:
            try: os.makedirs(datadir + "/" + prefix + algo + "/" + str(it))
            except: pass

        for n,d in confs:
            for ex in exps:
                for algo in algos:
                    fn = "%s/%s/%d/a%s_n%d_d%d_e%.1f.log" % (datadir, prefix + algo, it, algo, n, d, ex)
                    print('Process ' + fn)

                    if (algo != "mh" and n > MAX_NO_NODES):
                        print(" ... skip due to MAX_NO_NODES constraint")
                        continue

                    if (n * d > MAX_NO_EDGES):
                        print(" ... skip due to MAX_NO_EDGES constraint")
                        continue

                    if os.path.isfile(fn):
                        print(" ... skip since file already exists")
                        continue

                    with open(fn, 'w') as outf:
                        outf.write("Hostname: %s" % socket.gethostname())

                        if ("mh" == algo):
                            invokeGenerator(n, d, ex, outf, worker=cpus, segments=(1 if d < 100 else 2) * cpus)
                        elif("nkorg" == algo):
                            invokeNkGenerator(n, d, ex, 0, outf)
                        elif("nkopt" == algo):
                            invokeNkGenerator(n, d, ex, 1, outf)
                        else:
                            assert(False)

cpus = multiprocessing.cpu_count()
print("Number of threads: %d" % cpus)
algos = ['nkorg', 'nkopt', 'mh']

# measure runtime and memory consumption as function of number of nodes
nodes = exp10Series(1e5, 6, 1)
degrees = [10, 1000]
exps = [3.0, 2.1]
confs = list(itertools.product(nodes, degrees))
confs.sort(key=lambda x: x[0]*x[1])
prefix = ""
runSeries("", iterations, confs, exps, algos, cpus)

# measure runtime and memory consumption as function of number of nodes
nodes = [int(1e8)]
degrees = exp10Series(10, 2, 3)
exps = [3.0]
confs = list(itertools.product(nodes, degrees))
confs.sort(key=lambda x: x[0]*x[1])
prefix = ""
runSeries("degs_", iterations, confs, exps, algos, cpus)
