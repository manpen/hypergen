#include "BandSegment.hpp"
#include <iostream>
#include <utility>
#include <tuple>
#include <algorithm>
#include "Assert.hpp"

BandSegment::BandSegment(Node firstNode, Count nodes,
                         const CoordInter phiRange, CoordInter radRange,
                         const Geometry& geometry,
                         const Configuration& config,
                         Coord_b streamingBound,
                         Seed seed,
                         unsigned int bandIdx,
                         unsigned int firstStreamingBand
)
    : _firstNode(firstNode)
    , _numberOfNodes(nodes)
    , _geometry(geometry)
    , _config(config)
    , _phiRange(phiRange)
    , _radRange(radRange)
    , _upperLimit(radRange.second)
    , _batchSize( (bandIdx >= firstStreamingBand) * geometry.avgDeg )
    , _active(geometry.avgDeg)
    , _maxDeltaPhi(_geometry.deltaPhi(_radRange.first, _radRange.first))
    , _bandIdx(bandIdx)
    , _propagatedUntil(_phiRange.first)
    , _coverFullCircle(std::abs(_phiRange.second - _phiRange.first) > 2*M_PI - 1e-6)
    , _streamingBound(streamingBound)
    , _seed(seed)
{}

void BandSegment::clear() {
    _points.clear(); _points.shrink_to_fit();
    _requests.clear(); _requests.shrink_to_fit();
    _generator.reset(nullptr);
    _active.clear();
}


// FIX-ME: Check whether it is more efficient to NOT de/allocate vector but rather use
//         member vectors
// TODO: We should check whether there at most (1+epsilon)_batchSize left
void BandSegment::generatePoints() {
    // generate points and merge new requests into _active
    _generator->generate(_batchSize, [&] (const Point& pt, const Request& req) {
        _points.push_back(pt);
        _active.addRequest(req);
    });
    sortPoints();


    // points are in general not sorted; we do it here as a preparation for SIMD accesses
    _no_valid_requests = _points.empty();
}

void BandSegment::generateGlobalPoints() {
    // generate points and merge new requests into _active
    _generator->generate(_numberOfNodes, _points, _requests);

    // points are in general not sorted; we do it here as a preparation for SIMD accesses
    _no_valid_requests = _points.empty();
}

void BandSegment::sortPoints() {
    std::sort(_points.begin(), _points.end(), [] (const auto& a, const auto &b) {
        return std::tie(a.phi, a.id) < std::tie(b.phi, b.id);}
    );
}

void BandSegment::copyGlobalState(const BandSegment& band, Coord deltaPhi) {
    const Coord threshold = _phiRange.first + deltaPhi;
    _active.copyFrom(band._active, false, threshold);
}

Coord_b BandSegment::prepareEndgame(BandSegment& band) {
    Coord_b maxPhi = _phiRange.first;
    _checkInvariants();

// copy points
    //assert(_points.empty());
    for(auto pt : band.getPoints()) {
        if (pt.phi > 2*M_PI) {
            pt.phi -= 2 * M_PI;
            ASSERT_EQ(_phiRange.first, 0.0);
        }

        _points.push_back(pt);
        _points.back().setOld();

        ASSERT_GE(pt.phi, _phiRange.first);
        ASSERT_LE(pt.phi, _phiRange.second);

        if (pt.phi > maxPhi)
            maxPhi = pt.phi;
    }

    band.getPoints().clear();
    band.getPoints().shrink_to_fit();

    //band.getActive().clear(true);
    _active.copyFrom(band._active, true);
    if (_phiRange.first < 1e-10)
        _active.fixRange(2*M_PI);

    _checkInvariants();

    return std::max<Coord_b>(maxPhi, _active.maxReplayRange());
}


#ifndef NDEBUG
void BandSegment::_checkInvariants() const {
    ASSERT(std::is_sorted(_points.cbegin(), _points.cend(),
        [] (const auto& a, const auto &b) {return a.phi < b.phi;}));
};
#endif
