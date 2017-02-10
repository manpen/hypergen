/**
 * @file
 * @brief Implementation of Geometry
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
#include "Geometry.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

Coord Geometry::_computeTargetRadius(const Count n, const Coord avgDeg, const Coord alpha) const {
    // algorithm adopted from NetworKIT
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
