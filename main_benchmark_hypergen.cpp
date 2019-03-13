/**
 * @file
*
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2019 Manuel Penschuck
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <algorithm>
#include <numeric>
#include <iostream>
#include <omp.h>
#include <cmath>

#include "include/Generator.hpp"
#include "include/Configuration.hpp"

#ifdef __unix__
#include "unistd.h"

static std::string hostname() {
    char tmp[128];
    if (gethostname(tmp, 128))
        return "n/a";
    return {tmp};
}

#else
static std::string hostname() {
    return "n/a";
}
#endif



double benchmark(std::ostream& os, const std::string& host, unsigned iter, unsigned int n, unsigned int avgDeg, double alpha, unsigned int seed = 0) {
    double time_total;
    EdgeId num_edges;
    {
        Configuration config;
        config.noSegments = config.noWorker = omp_get_max_threads();
        config.nodes = n;
        config.avgDegree = avgDeg;
        config.alpha = alpha;
        config.degreeExp = 2.0 * alpha + 1;
        config.seed = seed;
        config.R = -1; // trigger computation of R
        config.dump();

        ScopedTimer tot_timer(time_total);

        Generator gen(config);
        {
            ScopedTimer timer("Generator time");

            auto addPoint = [&](const Point &pt, unsigned int s) {};
            auto addEdge = [&](Edge /*e*/, unsigned int segmentId) {};
            gen.generate(addEdge, addPoint);
        }
        num_edges = gen.num_edges();
    }

    // Logging
    {
        std::stringstream ss;
        ss << "[CSV]"
           << host << ","
           << iter << ","
           << "HyperGen,"
           << n << ","
           << avgDeg << ","
           << alpha << ","
           << "0.0,"
           << time_total << ","
           << num_edges << ","
           << (2.0 * num_edges / n);

        os << ss.str() << std::endl;
    }

    return time_total;
}

int main(int argc, char* argv[]) {
    std::cout <<
              "Build configuration\n"
              " sizeof(Coord): "   << sizeof(Coord) << "\n"
              " sizeof(Node): "    << sizeof(Node) << "\n"
              " sizeof(Point): "   << sizeof(Point) << "\n"
              " sizeof(Request): " << sizeof(Request) << "\n"
              " sizeof(DefaultPrng):" << sizeof(DefaultPrng) <<
              std::endl;

    std::cout << " SIMD: NodePacking=" << NodePacking << " CoordPacking=" << CoordPacking << std::endl;
    {
        std::cout << " Log-Transform: "
                     #ifdef LOG_TRANSFORM
                     "ENABLED"
                     #else
                     "DISABLED"
                  #endif
                  << std::endl;
    }

    std::cout << "\n\n";



    // Print Header
    std::cerr << "[CSV]"
       "host,"
       "iter,"
       "algo,"
       "n,"
       "avgDeg,"
       "alpha,"
       "T,"
       "TimeTotal,"
       "GenNumEdge,"
       "GenAvgDeg\n";

    unsigned seed = 0;
    const auto host = hostname();

    const unsigned n0 = 1e4;
    const unsigned nMax = 1e8;
    const unsigned steps_per_dec = 3;
    const double timeout = 100 * 1e3; // ms

    for(int iter = 0; iter < 5; iter++) {
        for (double ple : {2.2, 3.0}) {
            const auto alpha = (ple - 1.0) / 2.0;
            unsigned int skip_n = nMax + 1;
            for (const auto avgDeg : {10, 100, 1000}) {

                int ni = 0;
                for (auto n = n0; n <= nMax; n = n0 * std::pow(10.0, 1.0 * ni / steps_per_dec), ++ni) {
                    std::cout << "\033[31miter=" << iter << ", PLE=" << ple << ", n=" << n << ", avgDeg=" << avgDeg
                              << "\033[0m\n";

                    if (avgDeg * 20 > n) continue;

                    std::cout << "iter=" << iter << ", n=" << n << ", avgDeg=" << avgDeg << "\n";

                    double time;
                    if (n < skip_n) {
                        time = benchmark(std::cerr, host, iter, n, avgDeg, alpha, seed);
                        if (time > timeout) {
                            skip_n = n;
                            std::cout << " took too long\n";
                        }
                    } else {
                        std::cout << " skip_n = " << skip_n << "\n";
                    }
                    seed += 10;
                }
            }
        }
    }

    return 0;
}
