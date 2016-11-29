#pragma once
#ifndef POINT_GENERATOR_HPP
#define POINT_GENERATOR_HPP

#include <random>

#include "Geometry.hpp"
#include "Summation.hpp"
#include "Assert.hpp"

class PointGenerator {
public:    
    PointGenerator(Node id0, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed) 
        : _geometry(geometry)
        , _phiRange(phiRange)
        , _radRange(radRange)
        , _random(seed)
        
        , _nodeId(id0)
        , _nodesLeft(nodes+1)
        
        // radial
        , _distrRad(geometry.radCdf(radRange.first), geometry.radCdf(radRange.second))
        , _paramsRad(_geometry.coshAlphaR - 1.0, 1.0 / geometry.alpha)
        
        // angular
        , _paramAngular(phiRange.first - phiRange.second, phiRange.second)
        , _minDeltaPhi(_geometry.deltaPhi(radRange.second, radRange.second))
        
        , _nextData(_computeNextData())
    {}

    
    /**
     * @brief
     * Generates the requests number of points and requesets.
     * @param noPoints p_noPoints: Number of points to be generated: 0 means all remaining
     */
    std::pair<std::vector<Point>, std::vector<Request>> generate(Count noPoints);

    void generate(Count noPoints, std::vector<Point>& points, std::vector<Request>& requests);



    /**
     * @brief 
     * Number of points not yet output
     */
    Count nodesLeft() const {
        return _nodesLeft;
    }

    Node nextNodeId() const {
        return _nodeId;
    }

    Coord nextRequestLB() const {
        return  _nodesLeft
          ? (_paramAngular.second + _paramAngular.first * std::exp(_sumAngular.sum()))
          : _phiRange.second;
    }

    Coord nextPointLB() const {
        return  _minDeltaPhi + nextRequestLB();
    }


private:
    const Geometry& _geometry;
    const CoordInter _radRange;
    const CoordInter _phiRange;

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
    const Coord _minDeltaPhi;

    // Tmp
    std::pair<Point, Request> _nextData;
    
    std::pair<Point, Request> _computeNextData();
};

#endif
