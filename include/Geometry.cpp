#include "Geometry.hpp"
#include <cmath>
#include <cassert>

// borrowed from NetworKit
Coord Geometry::_computeTargetRadius(const Count n, const Coord avgDeg, const Coord alpha) const {
    const Coord gamma = 2*alpha+1;
    const Coord xi = (gamma-1)/(gamma-2);
    const Coord xiInv = ((gamma-2)/(gamma-1));
    const Coord v = avgDeg * (M_PI/2)*xiInv*xiInv;
    
    constexpr Coord epsilon = 1e-10;
    
    Coord currentR = 2*log(n / v);
    Coord lowerBound = currentR/2;
    Coord upperBound = currentR*2;
    
    assert(getExpectedDegree(lowerBound, alpha, n) > avgDeg);
    assert(getExpectedDegree(upperBound, alpha, n) < avgDeg);
    
    do {
        currentR = (lowerBound + upperBound)/2;
        const Coord currentK = getExpectedDegree(currentR, alpha, n);
        if (currentK < avgDeg) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while(std::abs(getExpectedDegree(currentR, alpha, n) - avgDeg) > epsilon);
    
    return currentR;    
}
