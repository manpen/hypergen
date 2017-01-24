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
