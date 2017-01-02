#include "BandSegment.hpp"
#include <iostream>
#include <utility>
#include <tuple>
#include <algorithm>
#include "Assert.hpp"

template <typename T>
static void merge_vectors(std::vector<T>& srcdest, const std::vector<T>& add, std::vector<T>& helper) {
    helper.resize(srcdest.size() + add.size());
    std::merge(srcdest.cbegin(), srcdest.cend(), add.cbegin(), add.cend(), helper.begin());
    srcdest.swap(helper);
}


BandSegment::BandSegment(Node firstNode, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed)
    : _firstNode(firstNode)
    , _numberOfNodes(nodes)
    , _geometry(geometry)
    , _phiRange(phiRange)
    , _radRange(radRange)
    , _upperLimit(radRange.second)
    , _batchSize(1*geometry.avgDeg)
    , _generator(_firstNode, _numberOfNodes, _phiRange, _radRange, _geometry, seed)
{}


// FIX-ME: Check whether it is more efficient to NOT de/allocate vector but rather use
//         member vectors
// TODO: We should check whether there at most (1+epsilon)_batchSize left
void BandSegment::generatePoints(BandSegment& bandAbove) {
    // generate points and merge new requests into _active
    std::vector<Request> newRequests;
    _generator.generate(_batchSize, _points, newRequests);

    // points are in general not sorted; we do it here as a preparation for SIMD accesses
    std::sort(_points.begin(), _points.end(), [] (const auto& a, const auto &b) {
        return std::tie(a.phi, a.id) < std::tie(b.phi, b.id);}
    );
    
    _no_valid_requests = _points.empty();
    _propagateRequests(newRequests, bandAbove);
    
    // merge all new requests into active
    merge_vectors(_active, newRequests, _merge_helper);
}

void BandSegment::generateGlobalPoints() {
    _generator.generate(0, _points, _active);

    std::sort(_points.begin(), _points.end(), [] (const auto& a, const auto &b) {
        return std::tie(a.phi, a.id) < std::tie(b.phi, b.id);}
    );
}

void BandSegment::_propagateRequests(const std::vector<Request>& reqs, BandSegment& bandAbove) {
    // in the up-most band, we don't have to propagate
    if (&bandAbove == this)
        return;
    
    // first compute new ranges ...
    std::vector<Request> mapped(reqs.size());
    std::transform(reqs.cbegin(), reqs.cend(), mapped.begin(), 
                   [&] (const Request r) {return Request(r, _geometry, _upperLimit);});
    
    std::sort(mapped.begin(), mapped.end(), // TODO: can we do without sorting?
              [] (const Request& a, const Request &b) {return a.range.first < b.range.first;});
    
    // ... and then merge them into the insertion buffer above
    if (bandAbove._insertion_buffer.empty()) {
        bandAbove._insertion_buffer.swap(mapped);
    } else {
        merge_vectors(bandAbove.getInsertionBuffer(), mapped, _merge_helper);
    }
}

void BandSegment::_mergeInsertionBuffer(BandSegment& bandAbove) {
    if (_insertion_buffer.empty())
        return;

    if (&bandAbove != this)
        _propagateRequests(_insertion_buffer, bandAbove);
    
    merge_vectors(_active, _insertion_buffer, _merge_helper);

    // save to _active and free _insertion_buffer
    _insertion_buffer.clear();
}

void BandSegment::copyGlobalState(const BandSegment& band, Coord deltaPhi) {
    const Coord threshold = _phiRange.first + deltaPhi;

// copy requests
    {
        ASSERT(_active.empty());
        ASSERT(_insertion_buffer.empty());
        ASSERT(band._insertion_buffer.empty());
        const auto end = std::upper_bound(band.getRequests().cbegin(), band.getRequests().cend(), threshold,
                                          [] (const Coord thresh, const Request& req) {return thresh < req.range.first;}
        );
        _active.insert(_active.begin(), band.getRequests().cbegin(), end);
   }
}

Coord BandSegment::prepareEndgame(const BandSegment& band) {
    Coord maxPhi = _phiRange.first;
    ASSERT(_insertion_buffer.empty());

// copy points
    //assert(_points.empty());
    for(auto pt : band.getPoints()) {
        if (pt.phi > 2*M_PI) {
            pt.phi -= 2 * M_PI;
            ASSERT_EQ(_phiRange.first, 0.0);
        }

        _points.push_back(pt);
        _points.back().setOld();

        if (pt.phi > maxPhi)
            maxPhi = pt.phi;
    }


// copy requests
    auto extractRequests = [&] (const std::vector<Request>& reqs) {
        std::vector<Request> newReqs;
        for(auto req : reqs) {
            if (req.range.second >= 2*M_PI) {
                req.range.first = 0;
                req.range.second -= 2*M_PI;

                ASSERT_EQ(_phiRange.first, 0.0);
            }

            if (req.phi >= 2*M_PI)
                req.phi -= 2*M_PI;

            if (req.range.second > _phiRange.first) {
                newReqs.push_back(req);
                newReqs.back().setOld();
            }

            if (req.range.second > maxPhi)
                maxPhi = req.range.second;
        }

        return newReqs;
    };

    merge_vectors(_active, extractRequests(band.getRequests()), _merge_helper);
    if (!band.getInsertionBuffer().empty())
        merge_vectors(_insertion_buffer, extractRequests(band.getInsertionBuffer()), _merge_helper);

    return maxPhi;
}


#ifndef NDEBUG
void BandSegment::_checkInvariants() const {
    ASSERT(std::is_sorted(_points.cbegin(), _points.cend(),
        [] (const auto& a, const auto &b) {return a.phi < b.phi;}));
    ASSERT(std::is_sorted(_active.cbegin(), _active.cend(),
        [] (const auto& a, const auto &b) {return a.range.first < b.range.first;}));
    ASSERT(std::is_sorted(_insertion_buffer.cbegin(), _insertion_buffer.cend(),
        [] (const auto& a, const auto &b) {return a.range.first < b.range.first;}));
};
#endif
