#include <iostream>
#include <NetworKit/generators/HyperbolicGenerator.h>
#include <NetworKit/auxiliary/Timer.h>
#include <NetworKit/graph/GraphBuilder.h>
#include <NetworKit/auxiliary/Parallel.h>



int main() {
    std::cout << "simple demonstration of NetworKit as a library\n";

    const uint64_t points = 10000000;
    const double avgdeg = 50.0;
    const double alpha = 3.0;
    const double T = 0.0;
    
    NetworKit::HyperbolicGenerator gen(points, avgdeg, alpha, T);
    
    const double R = gen.getR();
    const auto pointsData = gen.generatePoints(points, R, alpha, T);
 
    for(unsigned int m=2; m; --m) {
        for(unsigned int r=0; r<1; r++) {
            Aux::Timer timer;
            timer.start();

            std::string key;
            switch(m=1) {
                case 1:
                    key = "Org";
                    gen.generateColdOrig(pointsData.first, pointsData.second, R);
                    break;

                case 2:
                    key = "Opt1";
                    gen.generateColdOpt(pointsData.first, pointsData.second, R);
                    break;

                case 3:
                    key = "Opt2";
                    gen.generateColdOrig(pointsData.first, pointsData.second, R);
                    break;
                    
                default:
                    abort();
            }
            timer.stop();
            std::cout << "Method: " << key << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
        }
    }
    

    return 0;

}

