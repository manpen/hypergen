#pragma once
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "Definitions.hpp"

#include <cmath>
#include <tuple>

struct Geometry {
    const Coord alpha;
    const Coord avgDeg;
    const Coord R;

    const Coord coshR;
    const Coord coshAlphaR;
    const Coord poincareR;
    
    Geometry() = delete;
    
    Geometry(Coord alpha, Coord avgDeg, Coord R) 
        : alpha(alpha)
        , avgDeg(avgDeg)
        , R(R)
        , coshR(std::cosh(R))
        , coshAlphaR(std::cosh(alpha*R))
        , poincareR(coshR * 0.5 - 0.5)
    {}    
    
    Geometry(Count nodes, Coord avgDeg, Coord alpha)
        : Geometry(alpha, avgDeg, _computeTargetRadius(nodes, avgDeg, alpha))
    {}
    
    Coord radDensity(Coord r) const {
        return alpha * std::sinh(r*alpha) / (coshAlphaR - 1.0);
    }
    
    Coord radCdf(Coord r) const {
        return (std::cosh(alpha * r) - 1.0) / (coshAlphaR - 1.0);
    }
    
    Coord radInvCdf(Coord x) const {
        return std::acosh(x * coshAlphaR - x + 1) / alpha;
    }

    Coord deltaPhi(const SinhCosh pt, const SinhCosh band) const {
        Coord deltaPhiCos = (pt.cosh * band.cosh - coshR) * pt.invsinh * band.invsinh;
        if (deltaPhiCos >  1.0) deltaPhiCos =  1.0;
        if (deltaPhiCos < -1.0) deltaPhiCos = -1.0;
    
        return std::acos(deltaPhiCos);
    }
    
    
private: 
    Coord _computeTargetRadius(Count nodes, Coord avgDeg, Coord alpha) const;
};

#endif
