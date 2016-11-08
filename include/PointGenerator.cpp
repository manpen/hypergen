#include "PointGenerator.hpp"

std::pair<std::vector<Point>, std::vector<Request> > PointGenerator::generate(Count noPoints, const Coord threshold)
{
    std::vector<Point> points;
    points.reserve(noPoints);

    std::vector<Request> requests;
    requests.reserve(noPoints);
    
    if (noPoints > _nodesLeft)
        noPoints = _nodesLeft + 1;
    
    for(Count i=0; i < noPoints && _nextData.second.range.first < threshold; ++i) {
        points.push_back(_nextData.first);
        requests.push_back(_nextData.second);
        _nextData = _computeNextData();
    }
    
    return {points, requests};
}

std::pair<Point, Request> PointGenerator::_computeNextData()
{
    _sumAngular.push(std::log(_distrAngular(_random)) / _nodesLeft--);
    
    const Coord phi = _paramAngular.second + _paramAngular.first * std::exp(_sumAngular.sum()); 
    const Coord rad = std::acosh(_distrRad(_random) * _paramsRad.first + 1.0) * _paramAngular.second;
    
    const SinhCosh radSC(rad);
    const Coord deltaPhi = _geometry.deltaPhi(radSC, radSC);
    
    Point   pts(_nodeId, phi+deltaPhi, rad);
    Request req(_nodeId, phi+deltaPhi, {phi, phi+deltaPhi+deltaPhi}, pts.r);
    
    return {pts, req};
}

