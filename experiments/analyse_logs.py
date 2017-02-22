#!/usr/bin/env python3

# @file
# @brief Analyse Experiments: Generate Plots and Tables
#
# This script analyses log file collected with ./run_experiments.py
# It expects the log files at ../data/*, so please execute the script
# from WITHIN the experiments folder. The plots and tables are written
# into the folder ./results
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

import glob, os, sys
import re
import itertools

import pandas as pa
import numpy as np
import pylab as plt

MIN_ITERATIONS = 1
PLOT_TYPE="png" # pgf

# configure visuals
markers=['o', 's', 'd', 'D']
linestyles=['-', '--', ':', '-.']
colors = ["#ff0000", "#008000", "#000080", "#000000"]

def extractWTime(line):
    value = line.replace("Elapsed (wall clock) time (h:mm:ss or m:ss): ", "")
    res = re.match(r"\s*((\d+)h )?(\d+)m (\d+(\.\d+)?)s", value)
    if (res):
        hours, minutes, seconds, _ = [float(x) for x in res.groups(0)[1:]]
        return (3600.0 * hours + 60.0 * minutes + seconds)

    fields = [float(x) for x in value.split(":")]
    assert(len(fields) <= 3)
    while (len(fields) < 3):
        fields = [0.0] + fields
    return fields[0] * 3600 + fields[1] * 60 + fields[2]

def extractNum(str, t="int"):
    if (t == "float"):
        res=re.search(r"\s(\d.\d+e[+\-]\d+)", str)
        if (res):
            return float(res.groups()[0])


    nums = [s for s in str.replace("ms", "").split() if s.replace('.','').isdigit()]
    if not nums:
        return None

    assert(len(nums) == 1)

    num = nums[0]

    if "." in num or t == "float":
        return float(num)
    else:
        return int(num)

def insertUnique(di, key, value):
    assert(not key in di)
    di[key] = value


# parse log-files in ../data directory
def loadLogs(path="./runtime/"):
    files = glob.glob(path + "*/*/*.log")

    logs = []
    for file in files:
        fn = os.path.basename(file)

        row = {'group': file.split('/')[-3]}

        with open(file, "r") as fh:
            for line in fh:
                if "conf.nodes:" in line or "-n No. Nodes" in line:
                    insertUnique(row, "nodes", extractNum(line))
                    row['algo'] = 'HyperGen' + (" (Phi)" if "phi" in row["group"] else "")
                elif "conf.avgDegree:" in line or "-d Avg. Degree" in line:
                    insertUnique(row, "deg", extractNum(line))
                elif "conf.degreeExp:" in line or "-e Distr. Exp" in line:
                    insertUnique(row, "exp", extractNum(line))
                    row["exp"] = int(10.0*row["exp"] + 0.5)
                elif "conf.noSegments:" in line:
                    insertUnique(row, "segs", extractNum(line))
                elif "conf.noWorker:" in line:
                    insertUnique(row, "worker", extractNum(line))
                elif "TargetRadius: " in line:
                    insertUnique(row, "R", extractNum(line, "float"))
                elif "Number of edges: " in line:
                    insertUnique(row, "edges", extractNum(line))
                elif "Generator time Time elapsed" in line:
                    insertUnique(row, "gtime", extractNum(line, "float"))

                elif "[INFO ]: Produced " in line and " edges" in line:
                    insertUnique(row, "edges", extractNum(line))
                elif "[INFO ]: Required " in line and " compares" in line:
                    insertUnique(row, "compares", extractNum(line.split("compares")[0]))
                elif "Timer Time elapsed" in line:
                    insertUnique(row, "gtime", extractNum(line, "float"))

                elif "Maximum resident set size (kbytes):" in line:
                    insertUnique(row, "ressize", extractNum(line))
                elif "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                    insertUnique(row, "wtime", extractWTime(line))

                elif "[INFO ]: Start Original" in line:
                    row['algo'] = 'NkGen-Org'
                elif "[INFO ]: Start Opt" in line:
                    row['algo'] = 'NkGen-Opt'

        if not "gtime" in row:
            print("Skip '%s' as it seems incomplete" % file)
            continue

        logs.append(row)

    if not logs:
        return []

    data = pa.DataFrame(logs)
    print("Gathered %d files from %d groups" % (len(row), len(data.group.unique())) )
    return data

def latexFriendlyPlot(fig_width_pt = 246.0, ratio=(5**0.5-1.0)/2.0):
    # Get fig_width_pt from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*ratio      # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'font.size': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'figure.figsize': fig_size,
              "font.family": "serif",
              "font.serif": [],
              "pgf.preamble": [
                  r"\usepackage[utf8x]{inputenc}",
                  r"\usepackage[T1]{fontenc}",
              ]
              }
    plt.rcParams.update(params)

print("Load Logs")
data = loadLogs("../data/")
if (0 == len(data)):
    print("No valid data gathered. Did you execute this script from WITHIN")
    print("the folder experiments?")
    sys.exit(-1)

try: os.mkdir("results")
except: pass

def sortAlgos(algos):
    res = sorted(list(algos), key=lambda x: "A"+x if "Phi" in x else x, reverse=True)
    return res

# Plot runtime/memory consumption as function of number of nodes/edges
# Also extract information for speed-up table
def comparePlots(data, relative, enable_plots = True):
    latexFriendlyPlot(300)
    exps = data["exp"].unique()
    degs = data["deg"].unique()
    exps = [30, 21]
    degs = [10, 1000]

    if relative:
        objectives = {
            "ressize": ("Max. Res. Set Size [B] per node", 2**10, "nodes")
            ,"gtime": ("Walltime [ns] per edge ", 1e6, "edges")
                      }
    else:
        objectives = {
            "gtime": ("$t(n)$ Walltime [s]", 1e-3, "edges")
            ,"ressize": ("Max. Resident Set Size [MiB]", 2.0**-10, "nodes")}

    avg_rows = []
    for obj, (ylabel, yscale, xkey) in objectives.items():
        for (deg, dexp) in itertools.product(degs, exps):
            plot_data = data.loc[(data["deg"]==deg)&(data["exp"] == dexp)]

            print('-' * 80)

            if enable_plots:
                fig = plt.figure(1)
                plt.clf()

                plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])
                ax = plt.subplot()
                ax.set_xscale("log", nonposx='clip')
                ax.set_yscale("log", nonposy='clip')

            minx=1e100
            maxx=1

            for (algo, ls, col, mark) in zip(sortAlgos(plot_data["algo"].unique()), linestyles, colors, markers):
                if "ressize" == obj and "Phi" in algo:
                    continue

                try:
                    grp = plot_data[plot_data["algo"] == algo].groupby("nodes")
                    mask = (grp.size() >= MIN_ITERATIONS)
                    if (not np.all(mask)):
                        print("Skip %d of %d groups as they include too few entries" % (np.sum(np.invert(mask)), len(mask)))
                    if not np.any(mask):
                        continue


                    if xkey == "nodes":
                        xparam = grp["nodes"].median()[mask]
                    else:
                        xparam = grp.median()[xkey][mask]

                    minx = min(minx, np.min(xparam))
                    maxx = max(maxx, np.max(xparam[mask]))

                    obj_mean = grp.median()[obj][mask] * yscale
                    obj_err  = np.array(list(grp.std()[obj][mask] * yscale))

                    if relative:
                        obj_mean /= xparam
                        obj_err /= xparam

                        mean_mask = (xparam >= (9e7 * deg))

                        if (obj == "gtime"):
                            avg_rows.append ({
                                'algo': algo,
                                'xkey': xkey,
                                'exp': dexp,
                                'deg': deg,
                                'edgetime_mean': np.mean(obj_mean[mean_mask]),
                                'edgetime_std': np.std(obj_mean[mean_mask])
                            })

                    if enable_plots:
                        plt.errorbar(xparam, obj_mean, yerr=obj_err, label=algo, ls=ls, color=col, marker=mark, markersize=3)
                except:
                    pass

            title = "Deg: %d, Exp: %.1f" % (deg, dexp/10.0)
            print(title, xkey, obj)

            if enable_plots:
                if (False and obj == "ressize"):
                    plt.axhline(2**10 * 64)
                plt.title(title)
                plt.xlabel("Number $n$ of nodes" if xkey == "nodes" else "Number $m$ of edges")
                plt.ylabel(ylabel)
                ax.set_xlim([minx/1.2, maxx * 1.2 ])

            fn = 'results/plot_%s_d%d_e%d_%s.%s' % (obj, deg, dexp, "rel" if relative else "abs", PLOT_TYPE)
            if not fn in ["plot_gtime_d1000_e30_rel.pgf", "plot_ressize_d1000_e30_rel.pgf"]:
                plt.legend(loc="best")

            print(fn)

            if enable_plots:
                plt.savefig(fn, transparnt=True, bbox_inches="tight")

    return pa.DataFrame(avg_rows)

comparePlots(data.loc[data["group"].map(lambda x: x in ["nkorg", "nkopt", "mh", "phi"])], False)
avgs = comparePlots(data.loc[data["group"].map(lambda x: x in ["nkorg", "nkopt", "mh", "phi"])], True, True)

def generateSpeedupTable(avgs):
    latex = r"\begin{tabular}{|l|" + ("|r|r|" * 4)  + "}\n"
    latex += "\t" + r"\hline" + "\n"
    latex += "\t" + r"Implementation "
    latex +=    r"& \multicolumn{2}{c||}{$\bar d {=} 10^1,\ \alpha{=}0.55$}"
    latex +=    r"& \multicolumn{2}{c||}{$\bar d {=} 10^1,\ \alpha{=}1$}"
    latex +=    r"& \multicolumn{2}{c||}{$\bar d {=} 10^3,\ \alpha{=}0.55$}"
    latex +=    r"& \multicolumn{2}{c||}{$\bar d {=} 10^3,\ \alpha{=}1$}"
    latex +=    r"\\\hline" + "\n"
    latex += "\t" + (r"  & {\footnotesize Time/edge} & {\footnotesize Speed-up} " * 4) + r"\\\hline" + "\n"

    refAlgo = "NkGen-Org"
    for algo in ["NkGen-Org", "NkGen-Opt", "HyperGen"]:
        latex += "\t" + algo + " \n"
        for deg in [10, 1000]:
            for dexp in [21, 30]:
                row = avgs[(avgs.algo == algo) & (avgs.exp == dexp) & (avgs.deg == deg)]
                time = row.edgetime_mean.values[0]
                std  = row.edgetime_std.values[0]
                latex += "\t\t & $%.1f \pm %.1f$" % (time, std)

                speedup = "n/a"
                if (algo != refAlgo):
                    refrow = avgs[(avgs.algo == refAlgo) & (avgs.exp == dexp) & (avgs.deg == deg)]
                    reftime = refrow.edgetime_mean.values[0]
                    refstd  = row.edgetime_std.values[0]

                    speedup = "$%.1f \pm %.1f$" % (reftime / time, refstd / time + std * reftime / time / time)

                latex += " & " + speedup + "\n"


        latex += "\t" + r"\\\hline" "\n"

    latex += r"\end{tabular}" + "\n"
    return latex

with open("results/speedup.tex", "w") as f:
    f.write(generateSpeedupTable(avgs))

# Print runtime as function of average degree
def compareDegPlots(data, relative):
    latexFriendlyPlot(300)
    exps = data["exp"].unique()
    unodes = data["nodes"].unique()

    if relative:
        objectives = {"gtime": ("Walltime [ns] per edge ", 1e6, "edges")}
    else:
        objectives = {"gtime": ("$t(n)$ Walltime [s]", 1e-3, "deg")}

    for obj, (ylabel, yscale, xkey) in objectives.items():
        for (nodes, dexp) in itertools.product(unodes, exps):
            plot_data = data.loc[(data["nodes"]==nodes)&(data["exp"] == dexp)]

            fig = plt.figure(1)
            plt.cla()
            plt.clf()

            plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])
            ax = plt.subplot()
            ax.set_xscale("log", nonposx='clip')
            ax.set_yscale("log", nonposy='clip')

            for (algo, ls, col, mark) in zip(sorted(plot_data["algo"].unique())[::-1], linestyles, colors, markers):
                try:
                    grp = plot_data[plot_data["algo"] == algo].groupby("deg")
                    mask = grp.size() > MIN_ITERATIONS
                    if (not np.all(mask)):
                        print("Skip %d of %d groups as they include too few entries" % (np.sum(np.invert(mask)), len(mask)))

                    if not np.any(mask):
                        continue



                    if xkey == "deg":
                        xparam = grp["deg"].median()[mask]
                    else:
                        xparam = grp.median()[xkey][mask]


                    obj_mean = grp.mean()[obj][mask] * yscale
                    obj_err  = np.array(list(grp.std()[obj][mask] * yscale))

                    if relative:
                        obj_mean /= xparam
                        obj_err /= xparam
                        print(algo, np.mean(obj_mean))

                    obj_err = 0

                    plt.errorbar(xparam, obj_mean, yerr=obj_err, label=algo, ls=ls, color=col, marker=mark, markersize=3)
                except:
                    pass

            if (False and obj == "ressize"):
                plt.axhline(2**10 * 64)

            plt.title("Nodes: $10^{%d}$, Exp: %.1f" % ( np.log10(nodes), dexp/10.0))
            plt.xlabel(r"Average Degree $\bar d$")
            plt.ylabel(ylabel)

            fn = 'results/plot_deg_%s_n%d_e%d_%s.%s' % (obj, nodes, dexp, "rel" if relative else "abs", PLOT_TYPE)
            plt.legend(loc="best")

            plt.savefig(fn, transparnt=True, bbox_inches="tight")

compareDegPlots(data.loc[data["group"].map(lambda x: "degs_" in x)], False)