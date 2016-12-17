#include <iostream>
#include <omp.h>

#include "include/Generator.hpp"

#ifndef NDEBUG
#define CROSS_REFERENCE
#endif

#include <iostream>
#ifdef CROSS_REFERENCE
    #include <NetworKit/generators/HyperbolicGenerator.h>
    #include <NetworKit/auxiliary/Timer.h>
    #include <NetworKit/graph/GraphBuilder.h>
    #include <NetworKit/auxiliary/Parallel.h>
#endif

#include <cmath>

#include <parallel/algorithm>
#include <string>


int main(int argc, char* argv[]) {
    unsigned int noWorker = omp_get_max_threads();
    
    unsigned int confNoPoints = 1000;
    double confAvgDeg = 5;
    double confAlpha = 2.1;
    Seed confSeed = 1234;

    for(unsigned int i=1; i+1 < argc; i+=2) {
        std::string key = argv[i];
        std::string value = argv[i+1];

        if      (key == "-n") {confNoPoints = stoll(value);}
        else if (key == "-d") {confAvgDeg = stod(value);}
        else if (key == "-a") {confAlpha = stod(value);}
        else if (key == "-s") {confSeed = stoi(value);}
        else if (key == "-w") {noWorker = stoi(value);}
        else {
            std::cerr << "Unknown argument: " << key << std::endl;
            abort();
        }
    }

    std::cout << "Parameters:\n"
              "-n No. Nodes   " << confNoPoints << "\n"
              "-d Avg. Degree " << confAvgDeg << "\n"
              "-a Alpha       " << confAlpha << "\n"
              "-w No. Worker  " << noWorker << "\n"
              "-s Seed        " << confSeed
    << std::endl;

    const auto threadsBefore = omp_get_max_threads();

    std::cout << "SIMD: NodePacking=" << NodePacking << " CoordPacking=" << CoordPacking << std::endl;
#ifdef CROSS_REFERENCE
    std::cout << "Cross-Referencing: Enabled" << std::endl;
#else
    std::cout << "Cross-Referencing: Disabled" << std::endl;
#endif


// run new generator
    std::vector<Point> points;
#ifdef CROSS_REFERENCE
    points.resize(confNoPoints, Point(0, -100, 1));
#endif

    std::vector<EdgeId> nodeCounters(noWorker);
    std::vector<EdgeId> edgeCounters(noWorker);
    std::vector<std::vector<Edge>> edges(noWorker);

    Generator gen(confNoPoints, confAvgDeg, confAlpha, confSeed, noWorker);
    {
        ScopedTimer timer("Generator time");

        auto addPoint = [&](const Point &pt, unsigned int s) {
#ifdef CROSS_REFERENCE
            points.at(pt.id) = pt;
            nodeCounters.at(s)++;
#endif
        };

        auto addEdge = [&](Edge e, unsigned int segmentId) {
#ifdef CROSS_REFERENCE
            edgeCounters.at(segmentId)++;

            if (e.first > e.second) {
                std::swap(e.first, e.second);
            }

            assert(e.first != e.second);

            edges.at(segmentId).push_back(e);
#endif
        };

        gen.generate(addEdge, addPoint);
    }

#ifdef CROSS_REFERENCE
    const auto noPoints = std::accumulate(nodeCounters.cbegin(), nodeCounters.cend(), 0);
    std::cout << "Points generated: " << noPoints << std::endl;
    
    for(unsigned int i=0; i<points.size(); i++) {
        if (points[i].phi < -50) {
            std::cout << "Missing point " << i << std::endl;
        }
    }
#endif

    std::cout << "Set number of threads to: " << threadsBefore << std::endl;
    omp_set_num_threads(threadsBefore);

// combine edges
#ifdef CROSS_REFERENCE
    std::vector<Edge> myEdges;
    {
        for (unsigned int i = 0; i < noWorker; ++i)
            myEdges.insert(myEdges.end(), edges[i].cbegin(), edges[i].cend());
        __gnu_parallel::sort(myEdges.begin(), myEdges.end());
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
            rads.push_back(std::acosh(pt.r.cosh));
        }
        
        NetworKit::HyperbolicGenerator nkGen(confNoPoints, confAvgDeg, confAlpha, 0.0);
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
    {
        auto mi = myEdges.cbegin();
        auto ni = nkEdges.cbegin();
        EdgeId matches = 0;
        EdgeId missing = 0;
        EdgeId wrong  = 0;
        EdgeId repeats = 0;

        while(mi != myEdges.cend() || ni != nkEdges.cend()) {
            const bool mdone = (mi == myEdges.cend());
            const bool ndone = (ni == nkEdges.cend());

            auto print = [&] (bool my, bool nk) {
                if (my != nk) {
                    Edge e = my ? *mi : *ni;
                    std::cout << "["
                              << std::setw(10) << std::right << e.first << ", "
                              << std::setw(10) << std::right << e.second << "] "
                              << (my ? "My " : "   ")
                              << (nk ? "Nk " : "   ")
                              << "\n";
                }

                if (nk) ni++;
                if (my) mi++;
           };

            if (mdone || (!ndone && *ni < *mi)){
                print(false, true);
                missing++;
            } else if (ndone || *mi < *ni) {
                print(true, false);
                wrong++;
                if (mi != myEdges.cbegin()) {
                    const Edge& e = *(mi - 1);
                    std::cout << " Before " << e.first << ", " << e.second << "\n "
                              << points[e.first] << "\n "
                              << points[e.second] << ", "
                              << "distance: " << std::setw(20) << points[e.first].distanceTo(points[e.second]) << ", "
                              << "distance(H): " << std::setw(20) << points[e.first].distanceToHyper(points[e.second]) << ", "

                              << "R: " << std::setw(20) <<  gen.getGeometry().R
                              << std::endl;

                    repeats += (e == *mi);
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
            " Matches: " << matches << "\n"
            " Missing: " << missing << "\n"
            " Wrong:   " << wrong << "\n"
            " Repeat:  " << repeats
        << std::endl;
    }
#endif

    return 0;
}
