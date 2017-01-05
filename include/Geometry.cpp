#include "Geometry.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

// borrowed from NetworKit
Coord Geometry::_computeTargetRadius(const Count n, const Coord avgDeg, const Coord alpha) const {
    const double gamma = 2*alpha+1;
    const double xi = (gamma-1)/(gamma-2);
    const double xiInv = ((gamma-2)/(gamma-1));
    const double v = avgDeg * (M_PI/2)*xiInv*xiInv;
    
    constexpr double epsilon = 1e-10;
    
    Coord currentR = 2.0*log(n / v);
    Coord lowerBound = currentR/2;
    Coord upperBound = currentR*2;
    
    assert(getExpectedDegree(lowerBound, alpha, n) > avgDeg);
    assert(getExpectedDegree(upperBound, alpha, n) < avgDeg);

    unsigned int iteration = 0;
    double currentDev = 2*epsilon;

    do {
        ++iteration;

        currentR = (lowerBound + upperBound)/2;
        const double currentK = getExpectedDegree(currentR, alpha, n);
        currentDev = std::abs(getExpectedDegree(currentR, alpha, n) / avgDeg - 1.0);

        if (0) {
            std::cout << "ComputeTargetRadius " << iteration
                      << " CurrentR: " << currentR
                      << " CurrentK: " << currentK
                      << " Eps: " << currentDev << std::endl;
        }

        if (iteration > 100) {
            std::cerr << "[WARNING] ComputeTargetRadius seems not to converge; use current value" << std::endl;
            return currentR;
        }

        if (currentK < avgDeg) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while(currentDev > epsilon);
    
    return currentR;    
}
