/**
 * @file
 * @brief PointGenerator
 *
 * Provides a random point/request stream with increasing requests
 * beginnings. Radial and polar interval may be constrained.
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
#pragma once
#ifndef POINT_GENERATOR_HPP
#define POINT_GENERATOR_HPP

#include <random>

#include "Definitions.hpp"
#include "Geometry.hpp"
#include "Point.hpp"
#include "Summation.hpp"
#include "Assert.hpp"
#include "RandomHelper.hpp"

class PointGenerator {
public:    
    PointGenerator(Node id0, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, DefaultPrng rg)
        : _geometry(geometry)
        , _phiRange(phiRange)
        , _radRange(radRange)
        , _random(rg)
        
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
    void generate(Count noPoints, std::vector<Point>& points);

    template<typename Callback>
    void generate(Count noPoints, Callback cb) {
        // if the threshold is negative, we will ignore it and produce all points available
        if (!noPoints || noPoints > _nodesLeft)
            noPoints = _nodesLeft;

        if (!noPoints)
            return;


        for(Count i=0; i < noPoints; ++i) {
            cb(_nextData.first, _nextData.second);
            _nextData = _computeNextData();
        }
    }

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

    DefaultPrng _random;

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
    

    std::pair<Point, Request> _computeNextData() {
        ASSERT_GT(_nodesLeft, 0);

        if (!--_nodesLeft)
            return {Point{}, Request{}};

        _sumAngular.push(std::log(_distrAngular(_random)) / _nodesLeft);

        const Coord phi = _paramAngular.second + _paramAngular.first * std::exp(_sumAngular.sum());
        ASSERT_GE(phi, _phiRange.first);
        ASSERT_LS(phi, _phiRange.second);

        const Coord rad = std::acosh(_distrRad(_random) * _paramsRad.first + 1.0) * _paramsRad.second;
        ASSERT_GE(rad, _radRange.first);
        ASSERT_LS(rad, _radRange.second);

        const SinhCosh radSC(rad);
        ASSERT_LS(radSC.cosh, _geometry.coshR);

        const Coord deltaPhi = _geometry.deltaPhi(radSC, radSC);

        Point   pts(_nodeId, phi+deltaPhi, rad);
        Request req(_nodeId, phi+deltaPhi, {phi, phi+deltaPhi+deltaPhi}, pts.r);

        _nodeId++;

        return {pts, req};
    }

    Point _computeNextPoint() {
        ASSERT_GT(_nodesLeft, 0);

        if (!--_nodesLeft)
            return {};

        _sumAngular.push(std::log(_distrAngular(_random)) / _nodesLeft);

        const Coord phi = _paramAngular.second + _paramAngular.first * std::exp(_sumAngular.sum());
        ASSERT_GE(phi, _phiRange.first);
        ASSERT_LS(phi, _phiRange.second);

        const Coord rad = std::acosh(_distrRad(_random) * _paramsRad.first + 1.0) * _paramsRad.second;
        ASSERT_GE(rad, _radRange.first);
        ASSERT_LS(rad, _radRange.second);

        const SinhCosh radSC(rad);
        ASSERT_LS(radSC.cosh, _geometry.coshR);

        const Coord deltaPhi = _geometry.deltaPhi(radSC, radSC);

        return {_nodeId++, phi+deltaPhi, rad};
    }
};

#endif
