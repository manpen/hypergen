#include <iostream>
#include <NetworKit/generators/HyperbolicGenerator.h>
#include "NetworKit/auxiliary/Timer.h"


int main() {
    std::cout << "simple demonstration of NetworKit as a library\n";

    const uint64_t points = 10000000;
    const double avgdeg = 50.0;
    const double alpha = 3.0;
    const double T = 0.0;
    
    NetworKit::HyperbolicGenerator gen(points, avgdeg, alpha, T);
    
    const double R = gen.getR();
    const auto pointsData = gen.generatePoints(points, R, alpha, T);
    
    for(unsigned int m=0; m<2; m++) {
        for(unsigned int r=0; r<3; r++) {
            Aux::Timer timer;
            timer.start();
            if (m) {
                gen.generateColdOpt(pointsData.first, pointsData.second, R);
            } else {
                gen.generateColdOrig(pointsData.first, pointsData.second, R);
            }
            timer.stop();
            std::cout << "Method: " << (m ? "Opt" : "Org") << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
        }
    }
    

    return 0;

}

