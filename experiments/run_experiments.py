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
from math import exp, pi, log
import multiprocessing
import socket
import itertools
import random

# decrease these constants to speed-up simulation
MAX_NO_EDGES = int(1e14)
MAX_NO_NODES = int(1e8)
iterations = 9

def getExpectedDegree(n, R, alpha):
    gamma = 2*alpha+1
    xi = (gamma-1)/(gamma-2)
    firstSumTerm = exp(-R/2)
    secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((pi/4)*pow((1/alpha),2)-(pi-1)*(1/alpha)+(pi-2))-1)
    return (2 / pi) * n * xi * xi *(firstSumTerm + secondSumTerm)

def getTargetRadius(n, alpha, avgDeg):
    gamma = 2*alpha+1
    xi = (gamma-1)/(gamma-2)
    xiInv = ((gamma-2)/(gamma-1))
    v = avgDeg * (pi/2)*xiInv*xiInv

    epsilon = 1e-10

    currentR = 2.0*log(n / v)
    lowerBound = currentR/2
    upperBound = currentR*2

    assert(getExpectedDegree(n, lowerBound, alpha) > avgDeg)
    assert(getExpectedDegree(n, upperBound, alpha) < avgDeg)

    currentDev = 2*epsilon

    for iteration in range(100):
        currentR = (lowerBound + upperBound)/2
        currentK = getExpectedDegree(n, currentR, alpha)
        currentDev = abs(getExpectedDegree(n, currentR, alpha) / avgDeg - 1.0)

        if currentDev < epsilon:
            return currentR

        if currentK < avgDeg:
            upperBound = currentR
        else:
            lowerBound = currentR

    print("[WARNING] ComputeTargetRadius seems not to converge; use current value")
    return currentR


def invokeGenerator(n, avgDeg, expv, outf, seed=-1, segments=-1, worker=-1, bandSpacing=2, skipDistComp=False):
    if not os.path.isfile("./main_hyper"):
        print("Did not find binary ./main_hyper -- please run this script")
        print("from within the build directory.")
        sys.exit(-1)

    if seed < 0:
        seed = random.randint(0, 0xfffffff)

    args = ["./main_hyper" + ("_skipdist" if skipDistComp else ""),
            "-n", n, "-d", avgDeg, "-e", expv, "-s", seed, "-b", bandSpacing]

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

def invokeNkGenerator(n, avgDeg, expv, mode, outf, seed=-1, segments=-1, worker=-1, skipDistComp=False):
    if not os.path.isfile("./benchmark_networkit"):
        print("Did not find binary ./benchmark_networkit")
        print("Please run this script from within the build directory.")
        sys.exit(-1)

    if seed < 0:
        seed = random.randint(0, 0xfffffff)

    args = ["./benchmark_networkit" + ("_skipdist" if skipDistComp else ""),
            "-n", n, "-d", avgDeg, "-e", expv, "-s", seed, "-g", mode]

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

def invokeEmbedderGenerator(n, avgDeg, expv, outf, seed=-1):
    if not os.path.isfile("../related_work/girg/embedder"):
        print("Did not find binary ./embedder")
        print("Please run this script from within the build directory.")
        sys.exit(-1)

    if seed < 0:
        seed = random.randint(0, 0xfffffff)

    R = getTargetRadius(n, expv, avgDeg)
    C = R - 2*log(n)

    args = ["../related_work/girg/embedder",
            "--generate", "dummy",
            "-T", 0,
            "-n", n,
            "-C", C,
            "-alpha", expv,
            "-seed", seed]

    args = [str(x) for x in args]
    print("Execute " + " ".join(args))
    args = ["/usr/bin/time", "-av"] + args
    subprocess.call(args, stdout=outf, stderr=outf)

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

                    if (algo != "mh" and algo != "mh_skipdist" and n > MAX_NO_NODES):
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

                        if ("mh" == algo or "mh_skipdist" == algo):
                            invokeGenerator(n, d, ex, outf, worker=cpus,
                                            segments=(1 if d < 100 else 2) * cpus,
                                            skipDistComp=("mh_skipdist" == algo))
                        elif("nkorg" == algo or "nkorg_skipdist" == algo):
                            invokeNkGenerator(n, d, ex, 0, outf,
                                              skipDistComp=("nkorg_skipdist" == algo))
                        elif("nkopt" == algo):
                            invokeNkGenerator(n, d, ex, 1, outf)
                        elif("emb" == algo):
                            invokeEmbedderGenerator(n, d, ex, outf)
                        else:
                            assert(False)

cpus = multiprocessing.cpu_count()
print("Number of threads: %d" % cpus)
algos = ['nkorg', 'nkopt', 'mh', "emb"]

# measure runtime and memory consumption as function of number of nodes
nodes = exp10Series(1e5, 6, 1)
degrees = [10, 1000]
exps = [3.0, 2.1]
confs = list(itertools.product(nodes, degrees))
confs.sort(key=lambda x: x[0]*x[1])
prefix = ""
#runSeries("", iterations, confs, exps, algos, cpus)

# measure runtime and memory consumption as function of number of nodes
algos = ['nkorg', 'nkopt', 'mh', "emb", 'mh_skipdist', 'nkorg_skipdist']
nodes = [int(2**24)]
degrees = exp10Series(10, 2, 3)
exps = [3.0]
confs = list(itertools.product(nodes, degrees))
confs.sort(key=lambda x: x[0]*x[1])
prefix = ""
runSeries("degs_", iterations, confs, exps, algos, cpus)
