#include "Geometry.hpp"
#include <cmath>
#include <cassert>

// borrowed from NetworKit
double Geometry::_computeTargetRadius(const Count n, const double avgDeg, const double alpha) const {
    const double gamma = 2*alpha+1;
    const double xi = (gamma-1)/(gamma-2);
    const double xiInv = ((gamma-2)/(gamma-1));
    const double v = avgDeg * (M_PI/2)*xiInv*xiInv;
    
    constexpr double epsilon = 1e-10;
    
    double currentR = 2*log(n / v);
    double lowerBound = currentR/2;
    double upperBound = currentR*2;
    
    auto getExpectedDegree = [&] (double R) {
        double firstSumTerm = exp(-R/2);
        double secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((M_PI/4)*pow((1/alpha),2)-(M_PI-1)*(1/alpha)+(M_PI-2))-1);
        return (2.0 / M_PI) * xi * xi * n *(firstSumTerm + secondSumTerm);
    };
        
    
    assert(getExpectedDegree(lowerBound) > avgDeg);
    assert(getExpectedDegree(upperBound) < avgDeg);
    
    do {
        currentR = (lowerBound + upperBound)/2;
        const double currentK = getExpectedDegree(currentR);
        if (currentK < avgDeg) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while(std::abs(getExpectedDegree(currentR) - avgDeg) > epsilon);
    
    return currentR;    
}
