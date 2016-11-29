#pragma once
#ifndef BAND_SEGMENT_HPP
#define BAND_SEGMENT_HPP

#include <iostream>
#include <vector>
#include <algorithm>

#include "Geometry.hpp"
#include "PointGenerator.hpp"
#include "Assert.hpp"


class BandSegment {
public:
    BandSegment() = delete;
    BandSegment(Node firstNode, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed, bool endGame);

    // produce new points and request; delete exhausted requests from _active
    void generatePoints(BandSegment& bandAbove);
    void generateGlobalPoints();
    
    template <typename EdgeCallback>
    void generateEdges(EdgeCallback edgeCallback, BandSegment& bandAbove, double threshold = 2*M_PI) {
        if (_stats) {
            _stats_data.batches++;
            _stats_data.activeSize += _active.size();
        }

        if (_points.empty() || _points.front().phi > threshold) { //} _phiRange.second) {
            _no_valid_requests = true;
            return;
        }

        if (_verbose) {
            for (const auto &pt : getPoints())
                std::cout << " Pt:  " << pt << "\n";
            for (const auto &r: getRequests())
                std::cout << " Req: " << r << "\n";
            for (const auto &r: getInsertionBuffer())
                std::cout << " InB: " << r << "\n";
        }

        //if (threshold < 0)
        //
        const Coord bandThreshold = _generator.nodesLeft() ? (_generator.nextRequestLB() - _geometry.deltaPhi(_radRange.first, _radRange.first)) : _phiRange.second;
        threshold = std::min(threshold, bandThreshold);

        if (_verbose)
            std::cout << "  Threshold: " << threshold << std::endl;

        _checkInvariants();
        _mergeInsertionBuffer(bandAbove);
        const auto upper = _AosToSoa(_points.back().phi);
        const auto maxNode = _firstNode + _numberOfNodes;
        
        if (!upper) {
            _no_valid_requests = true;
            return;
        }
        
        _no_valid_requests = false;
        const Coord bandLowerCosh = std::cosh(_radRange.first);
        
        const auto startEnd = std::upper_bound(_req_phi_start.cbegin(), _req_phi_start.cbegin() + upper, _points.back().phi,
            [] (const Coord thresh, const Coord& start) {return thresh < start;}
        );
        
        if (_verbose) {
            std::cout << " generate edges with " 
                        << _points.size() << " points and " 
                        << _req_ids.size() << " requests of which " 
                        << std::distance(_req_phi_start.cbegin(), startEnd) << " will be considered" 
            << std::endl;
        }

        auto p = _points.cbegin();
        for(; p != _points.cend(); ++p) {
            const Point& pt = *p;
            if (pt.phi > threshold)
                break;

            const uint32_t this_upper = static_cast<uint32_t>(std::distance(_req_phi_start.cbegin(), 
                                       std::lower_bound(_req_phi_start.cbegin(), startEnd, pt.phi)));
            if (_verbose) {
                std::cout << "  compare point ("
                    "id: " << pt.id << ", "
                    "phi: " << pt.phi << ", "
                    "r: " << std::acosh(pt.r.cosh) << ") "
                    "with " << (this_upper) << " requests" 
                << std::endl;
            }

            const bool pointFromLastSegment = !_endgame || pt.old();

            for(uint32_t i = 0; i < this_upper; ++i) {
                const bool reqFromLastSegment = !_endgame || _req_old[i];
                //const bool ptIsSmaller = (_req_phi_start[i] <= pt.phi) && (_req_phi_stop[i] > pt.phi) && ( (_req_cosh[i] < bandLowerCosh) || (_req_ids[i] < pt.id) || pointFromLastSegment != reqFromLastSegment);
                const bool ptIsSmaller = (std::tie(_req_cosh[i], _req_phi_start[i]) < std::tie(pt.r.cosh, pt.phi)) && pt.id != _req_ids[i] && _req_phi_stop[i] > pt.phi;

                if (_stats) {
                    _stats_data.compares++;
                }

                if (_verbose) {
                    const Coord cosDist = (pt.r.cosh * _req_cosh[i] - _geometry.coshR) / (pt.r.sinh * _req_sinh[i]);
                    std::cout << "   request ("
                        "id: " << _req_ids[i] << ", "
                        "phi: " <<_req_phi[i] << ", "
                        "range: [" << _req_phi_start[i] << ", " << _req_phi_stop[i] << "], "
                        "r: " << std::acosh(_req_cosh[i]) << ") "
                        << " ptIsSmaller: " << ptIsSmaller
                        << " deltaPhi: " << _req_phi[i] - pt.phi
                        << " cosDist: " << cosDist
                        << " cos(deltaPhi): " << cos(_req_phi[i] - pt.phi)
                        << " oldPt: " << pointFromLastSegment
                        << " oldReq: " << reqFromLastSegment
                        << (ptIsSmaller && (cosDist < cos(_req_phi[i] - pt.phi)) ? " <------------" : "")
                        << (i == this_upper ? " !!!!!!!!" :"")
                        << "  ." << pt.id << "-" << _req_ids[i] << "." << _req_ids[i] << "-" << pt.id << "."

                    << std::endl;
                }

                //if (unlikely(deltaPhi > 2*M_PI)) deltaPhi -= 2*M_PI;
                //deltaPhi -= (deltaPhi > 2*M_PI) * 2.0 * M_PI;

                Coord deltaPhi = fabs(_req_phi[i] - pt.phi);
                const Coord dist = (pt.r.cosh * _req_cosh[i] - _geometry.coshR) / (pt.r.sinh * _req_sinh[i]);


                if (ptIsSmaller && (reqFromLastSegment || pointFromLastSegment)) {
                    if (dist < cos(deltaPhi)) {
                        edgeCallback({pt.id, _req_ids[i]});
                        if (_stats)
                            _stats_data.edges++;
                    }
                }
            }
        }

        // remove completed requests
        {
            const auto thresh = std::min<Coord>(_phiRange.second,
                        (p == _points.cend()) ? _generator.nextPointLB() : p->phi);
            if (_verbose) {
                std::cout << "Clean up with thresh = " << thresh << std::endl;
            }
            const auto end =  std::remove_if(_active.begin(), _active.end(), [&] (const Request& a) {return a.range.second < thresh;});
            _active.erase(end, _active.end());
        }

        // remove complete points
        if (p == _points.cend()) {
            _points.clear();
        } else {
            _points.erase(_points.begin(), p);
        }
    }

    void copyGlobalState(const BandSegment& band, Coord deltaPhi);
    Coord prepareEndgame(const BandSegment& band);

    template <typename EdgeCallback>
    void endgame(EdgeCallback edgeCallback) {

        ASSERT(!_active.empty());

        _checkInvariants();
        
        // compute the active request range, i.e. the request whos range ends the latest
        const Coord maxRange = std::max_element(_active.cbegin(), _active.cend(), [](const Request& a, const Request b) {return a.range.second < b.range.second;})->range.second;
        ASSERT_GT(maxRange, _phiRange.first);
        ASSERT_LE(maxRange, _phiRange.second);
        
        // estimate the current progress (lastPhi) and the number of points we still need to generate
        Coord lastPhi = _phiRange.first;
        while(!_active.empty()) {
            const Coord expDistance =  (_phiRange.second - _phiRange.first) / _numberOfNodes;
            ASSERT_LE(lastPhi, maxRange);
            
            const auto batchSize = std::max(100u, std::min<unsigned int>(10 * _batchSize, 
                1.1 * (maxRange - lastPhi) / expDistance
            ));
            
            // generate points and discard requests
            {
                std::vector<Request> newRequests;
                std::tie(_points, newRequests) = _generator.generate(batchSize);
            }
            
            if (_points.empty())
                return;
            
            lastPhi = _points.back().phi;

            // do the job ;)
            generateEdges(edgeCallback, *this);
        }
    }
    
    
    Node getFirstNode() const {return _firstNode;}
    Node getNumberOfNodes() const {return _numberOfNodes;}
    
    const CoordInter& getPhiRange() const {return _phiRange;}
    const CoordInter& getRadRange() const {return _radRange;}
    
    std::vector<Point>& getPoints() {return _points;}
    const std::vector<Point>& getPoints() const {return _points;}
    std::vector<Request>& getRequests() {return _active;}
    const std::vector<Request>& getRequests() const {return _active;}
    std::vector<Request>& getInsertionBuffer() {return _insertion_buffer;}
    const std::vector<Request>& getInsertionBuffer() const {return _insertion_buffer;}
    
    bool done() const {
        return _no_valid_requests && !_generator.nodesLeft();
    }

    bool done(const Coord thresh) const {
        return (_points.empty() || _points.front().phi > thresh) && (!_generator.nodesLeft() || _generator.nextRequestLB() > thresh);
    }

    struct Statistics {
        Count batches{0};
        EdgeId activeSize{0};
        EdgeId edges{0};
        EdgeId compares{0};
        
        Statistics operator+(const Statistics& o) const {
            Statistics s;
            
            s.batches    = batches    + o.batches;
            s.activeSize = activeSize + o.activeSize;
            s.edges      = edges      + o.edges;
            s.compares   = compares   + o.compares;
            
            return s;
        }
    };
    
    const Statistics& getStatistics() const {
        return _stats_data;
    }

    const Coord nextRequestLB() const {
        return _generator.nextRequestLB();
    }

private:
    // Geometry and Parameter
    const bool _endgame;
    const Node _firstNode;
    const Count _numberOfNodes;
    const Geometry& _geometry;
    const CoordInter _phiRange;
    const CoordInter _radRange;
    const SinhCosh   _upperLimit;
    
    const Count _batchSize;
    
    // Statistics
    static constexpr bool _verbose{false};
    static constexpr bool _stats{true};
    Statistics _stats_data;

    // Generator and AoS data
    PointGenerator _generator;

    std::vector<Request> _insertion_buffer;
    std::vector<Request> _active;
    bool _no_valid_requests{false};
    
    // Internal data structures (SoA)
    std::vector<Point> _points;
    std::vector<Node > _req_ids;
    std::vector<Coord> _req_phi;     
    std::vector<Coord> _req_phi_start;
    std::vector<Coord> _req_phi_stop; 
    std::vector<Coord> _req_sinh;     
    std::vector<Coord> _req_cosh;
    std::vector<bool>  _req_old;


    unsigned int _AosToSoa(const Coord thresh);
    void _mergeInsertionBuffer(BandSegment& bandAbove);
    void _propagateRequests(const std::vector<Request>& reqs, BandSegment& bandAbove);

#ifdef NDEBUG
    void _checkInvariants() const {}
#else    
    void _checkInvariants() const;
#endif
};

#endif
