/**
 * @file
 * @brief Implementation of PointGenerator
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
#include "PointGenerator.hpp"
#include <iostream>


void PointGenerator::generate(Count noPoints, std::vector<Point>& points, std::vector<Request>& requests) {
    // if the threshold is negative, we will ignore it and produce all points available
    if (!noPoints || noPoints > _nodesLeft)
        noPoints = _nodesLeft;

    if (!noPoints)
        return;

    // produce points
    points.reserve(points.size() + noPoints);
    requests.reserve(requests.size() + noPoints);

    for(Count i=0; i < noPoints; ++i) {
        points.push_back(_nextData.first);
        requests.push_back(_nextData.second);
        _nextData = _computeNextData();
    }
}


void PointGenerator::generate(Count noPoints, std::vector<Point>& points) {
    // if the threshold is negative, we will ignore it and produce all points available
    if (!noPoints || noPoints > _nodesLeft)
        noPoints = _nodesLeft;

    if (!noPoints)
        return;

    // produce points
    points.reserve(points.size() + noPoints);

    Point nextPoint = _nextData.first;
    for(Count i=0; i < noPoints; ++i) {
        points.push_back(nextPoint);
        nextPoint = _computeNextPoint();
    }

    _nextData.second = Request(nextPoint, _geometry, nextPoint.r);
}

std::pair<std::vector<Point>, std::vector<Request> > PointGenerator::generate(Count noPoints) {
    std::vector<Point> points;
    std::vector<Request> requests;

    generate(noPoints, points, requests);

    return {points, requests};
}
