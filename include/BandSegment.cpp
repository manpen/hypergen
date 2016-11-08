#include "BandSegment.hpp"
#include <utility>
#include <tuple>
#include <algorithm>
#include <cassert>

// FIX-ME: Check whether it is more efficient to NOT de/allocate vector but rather use
//         member vectors
void BandSegment::_prepare_advance(const Coord threshold, BandSegment& bandAbove) {
    // generate points and merge new requests into _active
    {
        std::vector<Request> newRequests;
        std::tie(_points, newRequests) = _generator.generate(_batchSize, threshold);
        
        // merge insertion buffer into new requests
        if (!_insertion_buffer.empty()) {
            std::vector<Request> merged(newRequests.size() + _insertion_buffer.size());
            std::merge(_active.cbegin(), _active.cend(), newRequests.cbegin(), newRequests.cend(), merged.begin());
            std::swap(newRequests, merged);
        }
        
        // propagate requests to band above (if not top-most band)
        if (_upperLimit.cosh + std::numeric_limits<Coord>::epsilon() < _geometry.coshR) {
            auto insBuf = bandAbove._insertion_buffer;
            assert(insBuf.empty());
            insBuf.resize(newRequests.size());
            
            for(unsigned int i=0; i < newRequests.size(); ++i) {
                insBuf[i] = Request(newRequests[i], _geometry, _upperLimit);
            }
        }
        
        // merge all new requests into active
        std::vector<Request> reqs(_active.size() + newRequests.size());
        std::merge(_active.cbegin(), _active.cend(), newRequests.cbegin(), newRequests.cend(), reqs.begin());
        std::swap(_active, reqs);
    }
    
    // transform AoS to SoA
    _req_ids      .resize(_active.size());
    _req_phi      .resize(_active.size());
    _req_phi_start.resize(_active.size());
    _req_phi_stop .resize(_active.size());
    _req_sinh     .resize(_active.size());
    _req_cosh     .resize(_active.size());
    
    for(unsigned int i=0; i < _active.size(); ++i) {
        const Request& req = _active[i];
        
        _req_ids      [i] = req.id;
        _req_phi      [i] = req.phi;
        _req_phi_start[i] = req.range.first;
        _req_phi_stop [i] = req.range.second;
        _req_sinh     [i] = req.r.sinh;
        _req_cosh     [i] = req.r.cosh;
    }
}
 
