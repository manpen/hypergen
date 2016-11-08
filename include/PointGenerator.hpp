#pragma once
#ifndef POINT_GENERATOR_HPP
#define POINT_GENERATOR_HPP

#include "Geometry.hpp"
#include "Summation.hpp"
#include <random>

class PointGenerator {
public:    
    PointGenerator(Node id0, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed) 
        : _geometry(geometry)
        , _random(seed)
        
        , _nodeId(id0)
        , _nodesLeft(nodes)
        
        // radial
        , _distrRad(geometry.radCdf(radRange.first), geometry.radCdf(radRange.second))
        , _paramsRad(std::cosh(geometry.alpha * geometry.R) - 1.0, 1.0 / geometry.R)
        
        // angular
        , _paramAngular(phiRange.first - phiRange.second, phiRange.second)
        
        , _nextData(_computeNextData())
    {}

    std::pair<std::vector<Point>, std::vector<Request>>
    generate(Count noPoints, const Coord threshold);
    
private:
    const Geometry& _geometry;
    
    std::mt19937_64 _random;

    Node  _nodeId;
    Count _nodesLeft;

    // Radial component
    std::uniform_real_distribution<Coord> _distrRad;
    // r = acosh(random * (cosh (alpha * R) - 1) + 1) / alpha
    //   = acosh(random * p.first) * p.second
    const std::pair<Coord, Coord> _paramsRad;
    
    // Angular component
    std::uniform_real_distribution<Coord> _distrAngular{0.0, 1.0};
    const std::pair<Coord, Coord> _paramAngular;
    KahnSummation<Coord> _sumAngular{0.0};

    // Tmp
    std::pair<Point, Request> _nextData;
    
    std::pair<Point, Request> _computeNextData();
};

#endif
