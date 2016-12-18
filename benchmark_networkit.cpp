#include <iostream>

#include <NetworKit/generators/HyperbolicGenerator.h>
#include <NetworKit/auxiliary/Timer.h>
#include <NetworKit/graph/GraphBuilder.h>
#include <NetworKit/auxiliary/Parallel.h>
#include <NetworKit/auxiliary/Random.h>

#include "include/ScopedTimer.hpp"


int main(int argc, char* argv[]) {
    unsigned int confNoPoints = 1000;
    double confAvgDeg = 5;
    double confAlpha = 2.1;
    unsigned int generatorAlgo = 0;
    unsigned int seed = 1234;

    for(unsigned int i=1; i+1 < argc; i+=2) {
        std::string key = argv[i];
        std::string value = argv[i+1];

        if      (key == "-n") {confNoPoints = stoll(value);}
        else if (key == "-d") {confAvgDeg = stod(value);}
        else if (key == "-a") {confAlpha = stod(value);}
        else if (key == "-g") {generatorAlgo = stoi(value);}
        else if (key == "-s") {seed = stoi(value);}
        else {
            std::cerr << "Unknown argument: " << key << std::endl;
            abort();
        }
    }

    Aux::Random::setSeed(seed, false );

    {
        ScopedTimer timer;
        NetworKit::HyperbolicGenerator gen(confNoPoints, confAvgDeg, confAlpha, 0.0);

        const double R = gen.getR();
        const auto pointsData = gen.generatePoints(confNoPoints, R, confAlpha, 0.0);

        std::string key;
        switch (generatorAlgo) {
            case 0:
                gen.generateColdOrig(pointsData.first, pointsData.second, R);
                break;

            case 1:
                gen.generateColdOpt(pointsData.first, pointsData.second, R);
                break;

            default:
                std::cerr << "Unknown Generator Algorithm." << std::endl;
                abort();
        }
    }

    return 0;

}

