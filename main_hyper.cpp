/**
 * @file
 * @brief Summation
 *
 * Sample implementation of HyperGen. If CROSS_REFERENCE is defined,
 * the points received are fed to the NetworKIT generator and results
 * are compared. SLOW!
 *
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2017 Manuel Penschuck
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
#include <iostream>
#include <omp.h>

//#define COUNT_NEIGHBORS
//#define DEGREE_DIST


#include "include/Generator.hpp"
#include "include/Configuration.hpp"


#include <iostream>
#ifdef CROSS_REFERENCE
    #include <NetworKit/generators/HyperbolicGenerator.h>
    #include <NetworKit/auxiliary/Timer.h>
    #include <NetworKit/graph/GraphBuilder.h>
    #include <NetworKit/auxiliary/Parallel.h>
    #include "include/Histogram.hpp"
#endif

#include <cmath>

#include <parallel/algorithm>
#include <string>



int main(int argc, char* argv[]) {
#if defined(CROSS_REFERENCE) && defined(SKIP_DIST_COMP)
    static_assert(false, "Cannot build CROSS_REFERENCE with SKIP_DIST_COMP");
#endif

    std::cout <<
        "main_hyper  Copyright (C) 2017 Manuel Penschuck\n"
        "This program comes with ABSOLUTELY NO WARRANTY;\n"
        "This is free software, and you are welcome to redistribute it\n"
        "under certain conditions\n"
    << std::endl;

    Configuration config(argc, argv);

    const auto threadsBefore = omp_get_max_threads();
    omp_set_num_threads(config.noWorker);
    std::cout <<
        "sizeof(Coord): "   << sizeof(Coord) << "\n"
        "sizeof(Node): "    << sizeof(Node) << "\n"
        "sizeof(Point): "   << sizeof(Point) << "\n"
        "sizeof(Request): " << sizeof(Request) << "\n"
        "sizeof(DefaultPrng):" << sizeof(DefaultPrng) <<
    std::endl;

    std::cout << "SIMD: NodePacking=" << NodePacking << " CoordPacking=" << CoordPacking << std::endl;
    {
        std::cout << "Cross-Referencing: "
        #ifdef CROSS_REFERENCE
            "ENABLED"
        #else
            "DISABLED"
        #endif
        << std::endl;
    }
    {
        std::cout << "Log-Transform: "
        #ifdef LOG_TRANSFORM
            "ENABLED"
        #else
            "DISABLED"
        #endif
        << std::endl;
    }


// run new generator
    std::vector<Point> points;
#ifdef CROSS_REFERENCE
    points.resize(config.nodes, Point(0, -100, 1));
#endif

#ifdef DEGREE_DIST
    std::vector<std::vector<uint32_t>> degrees(config.noSegments);
    for(auto& vec : degrees)
        vec.resize(config.nodes, 0);
#endif

    std::vector<EdgeId> nodeCounters(config.noSegments);
    std::vector<EdgeId> edgeCounters(config.noSegments);
    std::vector<std::vector<Edge>> edges(config.noSegments);
    std::vector<Node> nodeAccum(config.noSegments*8);

#ifdef COUNT_NEIGHBORS
    std::vector<Node> neighborhood(config.nodes, 0);
#endif

    Generator gen(config);
    {
        ScopedTimer timer("Generator time");

        auto addPoint = [&](const Point &pt, unsigned int s) {
#ifdef CROSS_REFERENCE
            points.at(pt.id) = pt;
            nodeCounters.at(s)++;
#endif
        };

        auto addEdge = [&](Edge e, unsigned int segmentId) {
#ifdef COUNT_NEIGHBORS
            neighborhood[e.first]++;
#endif

#ifdef DEGREE_DIST
            degrees[segmentId][e.first]++;
            degrees[segmentId][e.second]++;
#endif

#ifdef CROSS_REFERENCE
            edgeCounters.at(segmentId)++;

            if (e.first > e.second) {
                std::swap(e.first, e.second);
            }

            assert(e.first != e.second);

            edges.at(segmentId).push_back(e);
#else
            nodeAccum[8*segmentId] += (e.first + e.second);
#endif
        };

        gen.generate(addEdge, addPoint);
    }

#ifdef DEGREE_DIST
    {
        // sum all partial counts
        for(unsigned int i=1; i < config.noSegments; ++i) {
            auto dest = degrees[0].begin();
            for(auto src = degrees[i].cbegin(); src != degrees[i].cend(); ++src, ++dest)
                *dest += *src;
        }

        Histogram<true> hist;
        for(const auto & d : degrees[0]) {
            hist.addPoint(d);
        }

        hist.toStream(std::cout, "DEGREE-DIST");
    };
#endif


#ifdef CROSS_REFERENCE
    const auto noPoints = std::accumulate(nodeCounters.cbegin(), nodeCounters.cend(), 0);
    std::cout << "Points generated: " << noPoints << std::endl;
    
    for(unsigned int i=0; i<points.size(); i++) {
        if (points[i].phi < -50) {
            std::cout << "Missing point " << i << std::endl;
        }
    }
#endif


#ifdef COUNT_NEIGHBORS
    {
        Histogram<true> hist;
        for (const auto &d : neighborhood)
            hist.addPoint(d);
        hist.toStream(std::cout, "NBR-DIST-EXT");
    }
#endif

#ifndef CROSS_REFERENCE
    std::cout << "Accum key is " << std::accumulate(nodeAccum.cbegin(), nodeAccum.cend(), 0) << std::endl;
#endif

    std::cout << "Set number of threads to: " << threadsBefore << std::endl;
    omp_set_num_threads(threadsBefore);

// combine edges
#ifdef CROSS_REFERENCE
    std::vector<Edge> myEdges;
    {
        for (const auto& e : edges)
            myEdges.insert(myEdges.end(), e.cbegin(), e.cend());
        __gnu_parallel::sort(myEdges.begin(), myEdges.end());
    }


    // compute degree distribution
    if (0) {
        std::vector<Node> degrees(config.nodes, 0);
        for (const auto &e : myEdges) {
            degrees.at(e.first)++;
            degrees.at(e.second)++;
        }

        Histogram<true> hist;
        for (const auto &d : degrees)
            hist.addPoint(d);
        hist.toStream(std::cout, "DEG-DIST");
    }


// compare with networkit
    std::vector<Edge> nkEdges;
    {
        // map points to networkit rep
        std::vector<double> angles;
        std::vector<double> rads;
        angles.reserve(points.size());
        rads.reserve(points.size());
        
        for(const Point& pt : points) {
            angles.push_back(fmod(pt.phi, 2*M_PI));
            rads.push_back(pt.r.r);
        }

        NetworKit::HyperbolicGenerator nkGen(config.nodes, config.avgDegree, config.degreeExp, 0.0);
        nkGen.R = gen.getGeometry().R;
        auto graph = nkGen.generateColdOrig(angles, rads, gen.getGeometry().R);

        for(auto e : graph.edges()) {
            if (e.first > e.second) {
                std::swap(e.first, e.second);
            }

            nkEdges.push_back(e);
        }
        __gnu_parallel::sort(nkEdges.begin(), nkEdges.end());
    }

// print edges
    // detect duplicates
    Count duplicates = 0;
    {
        auto mb = myEdges.cbegin();
        auto mn = mb + 1;

        for(; mn != myEdges.cend(); ++mb, ++mn)
            duplicates += (*mb == *mn);
    }


    {
        auto mi = myEdges.cbegin();
        auto ni = nkEdges.cbegin();
        EdgeId matches = 0;
        EdgeId missing = 0;
        EdgeId missing_num = 0;
        EdgeId wrong  = 0;
        EdgeId wrong_num = 0;

        while(mi != myEdges.cend() || ni != nkEdges.cend()) {
            const bool mdone = (mi == myEdges.cend());
            const bool ndone = (ni == nkEdges.cend());

            auto print = [&] (bool my, bool nk) {
                bool numerical_issue = false;

                if (config.verbosity > 2 && my != nk) {
                    Edge e = my ? *mi : *ni;
                    std::cout << "["
                              << std::setw(10) << std::right << e.first << ", "
                              << std::setw(10) << std::right << e.second << "] "
                              << (my ? "My " : "   ")
                              << (nk ? "Nk " : "   ")
                              << "\n";
                }

                if (!my || !nk) {
                    const Edge e = my ? *mi : *ni;
                    numerical_issue = (std::abs(points[e.first].coshDistanceToPoincare(points[e.second]) - gen.getGeometry().coshR) /  gen.getGeometry().coshR < 1e-2);

                    if (config.verbosity > 2 && !numerical_issue)
                        std::cout << points[e.first] << "\n "
                                  << points[e.second] << ", "
                                  << "distance(P): " << std::setw(10) << points[e.first].distanceToPoincare(points[e.second]) << ", "
                                  << "distance(H): " << std::setw(10) << points[e.first].distanceToHyper(points[e.second]) << ", "
                                  << "R: " << std::setw(10) <<  gen.getGeometry().R
                                  << std::endl;
                }


                if (nk) ni++;
                if (my) mi++;

                return numerical_issue;
            };

            if (mdone || (!ndone && *ni < *mi)){
                missing_num += print(false, true);
                missing++;
            } else if (ndone || *mi < *ni) {
                if (mi != myEdges.cbegin()) {
                    wrong_num += print(true, false);
                    wrong++;
                }
            } else {
                assert(!mdone && !ndone && *mi == *ni);
                print(true, true);
                matches++;
            }
        }

        std::cout << "Edges produced\n"
            " My: " << myEdges.size() << "\n"
            " Nk: " << nkEdges.size() << "\n"
            " Ex: " << static_cast<EdgeId>(config.avgDegree*config.nodes/2) << "\n"
            " Matches:    " << matches << "\n"
            " Missing:    " << missing << "\n"
            " Missing(numerical issues): " << missing_num << "\n"
            " Wrong:   " << wrong << "\n"
            " Wrong(issues issues):   " << wrong_num << "\n"
            "Duplicates: " << duplicates << "\n"
        << std::endl;
    }
#endif

    return 0;
}
