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

std::pair<std::vector<Point>, std::vector<Request> > PointGenerator::generate(Count noPoints) {
    std::vector<Point> points;
    std::vector<Request> requests;

    generate(noPoints, points, requests);

    return {points, requests};
}

std::pair<Point, Request> PointGenerator::_computeNextData() {
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
