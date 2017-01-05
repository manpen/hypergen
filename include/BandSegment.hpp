#pragma once
#ifndef BAND_SEGMENT_HPP
#define BAND_SEGMENT_HPP

#include <iostream>
#include <vector>
#include <algorithm>

#include "Definitions.hpp"
#include "Point.hpp"
#include "Geometry.hpp"

#include "PointGenerator.hpp"
#include "Assert.hpp"
#include "Histogram.hpp"

#include <Vc/Vc>
#include <Vc/Allocator>

//#define SKIP_DIST_COMP

class BandSegment {
public:
    BandSegment() = delete;
    BandSegment(Node firstNode, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed);

    // produce new points and request; delete exhausted requests from _active
    void generatePoints(BandSegment& bandAbove);
    void generateGlobalPoints();
    
    template <bool Endgame, typename EdgeCallback>
    void generateEdges(EdgeCallback edgeCallback, BandSegment& bandAbove, Coord threshold = 2*M_PI) {
        _checkInvariants();
        _mergeInsertionBuffer(bandAbove);

        if (_stats) {
            _stats_data.activeSizes.addPoint( _active.size() );
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

        const Coord bandThreshold = _generator.nodesLeft()
            ? (_generator.nextRequestLB() - _geometry.deltaPhi(_radRange.first, _radRange.first))
            : _phiRange.second;
        threshold = std::min(threshold, bandThreshold);

        const auto upper = _AosToSoa<Endgame>(_points.front().phi, _points.back().phi);

        if (!upper) {
            _no_valid_requests = true;
            return;
        }

        const auto maxNode = _firstNode + _numberOfNodes;
        _no_valid_requests = false;
        const auto startEnd = std::upper_bound(_req_phi_start.cbegin(), _req_phi_start.cbegin() + upper, _points.back().phi,
            [] (const Coord thresh, const Coord_v& start) {return thresh < start[0];}
        );
        
        if (_verbose) {
            std::cout << " generate edges with " 
                        << _points.size() << " points and " 
                        << _req_ids.size() << " requests of which " 
                        << std::distance(_req_phi_start.cbegin(), startEnd) << " will be considered"
            << std::endl;
        }

        auto p = _points.cbegin();

        if (_stats) {
            _stats_data.pointSizes.addPoint( std::distance(p, _points.cend() ) );
        }

#ifdef POINCARE
    #ifndef LOG_TRANSFORM
        const Coord_v threshR = static_cast<Coord_b>(_geometry.poincareR);
    #endif
#else
        const Coord_v threshR = static_cast<Coord_b>(_geometry.coshR);
#endif

        auto binsearch_begin = _req_phi_start.cbegin();

        for(; p != _points.cend(); ++p) {
            const Point& pt = *p;
            if (pt.phi > threshold)
                break;

            const Coord_v ptrcosh(static_cast<Coord_b>(pt.r.cosh));
            const Coord_v ptphi(static_cast<Coord_b>(pt.phi));
#ifdef POINCARE
            const Coord_v poinX(static_cast<Coord_b>(pt.poinX));
            const Coord_v poinY(static_cast<Coord_b>(pt.poinY));
    #ifdef LOG_TRANSFORM
            const Coord_v threshR = static_cast<Coord_b>(_geometry.logPoincareR - pt.poinLogInvLen);
    #else
            const Coord_v poin_invr = static_cast<Coord_b>(pt.poinInvLen);
    #endif
#else
            const Coord_b ptrsinh = pt.r.invsinh;
#endif

#ifndef SKIP_DIST_COMP
            while(binsearch_begin != startEnd && (*binsearch_begin)[0] < pt.phi)
                ++binsearch_begin;

            const auto this_upper = static_cast<uint32_t>(
                    std::distance(_req_phi_start.cbegin(), binsearch_begin));

            if (_verbose) {
                std::cout << "  compare point ("
                    "id: " << pt.id << ", "
                    "phi: " << pt.phi << ", "
                    "r: " << std::acosh(pt.r.cosh) << ") "
                    "with " << (this_upper) << " requests"
                << std::endl;
            }

            const Coord_m pointFromLastSegment(!Endgame || pt.old());

            unsigned int noCandidates = 0;
            unsigned int noNeighbors = 0;

            noCandidates = this_upper * _req_phi[0].size();

#ifndef NDEBUG
            auto active_iter = _active.cbegin();
#endif

            for(unsigned int i = 0; i < this_upper; ++i) {
                auto prelimChecks =
                           (_req_phi_stop[i] > ptphi)
                        && (_req_phi_start[i] <= ptphi)
                        && ((_req_cosh[i] < ptrcosh) || (_req_cosh[i] == ptrcosh && _req_phi[i] < ptphi));
                
                if (Endgame)
                        prelimChecks = prelimChecks && (pointFromLastSegment || _req_old[i]);

                if (Statistics::enablePrelimCheck)
                    _stats_data.prelimCheck.addPoint( prelimChecks.count() );

                if (prelimChecks.isEmpty()) {
                    continue;
                }


                // computations
#ifdef POINCARE
                const auto deltaX = _req_poin_x[i] - poinX;
                const auto deltaY = _req_poin_y[i] - poinY;

                #ifdef LOG_TRANSFORM
                const auto dist = deltaX*deltaX + deltaY*deltaY;
                const auto isEdge = prelimChecks && (dist < Vc::exp(threshR - _req_poin_loginvr[i]));
                #else
                const auto dist = (deltaX*deltaX + deltaY*deltaY) * _req_poin_invr[i] * poin_invr;
                const auto isEdge = prelimChecks && (dist < threshR);
                #endif

#else
                const auto dist = (_req_cosh[i] * ptrcosh - threshR) * (_req_invsinh[i] * ptrsinh);
                const auto cosDist = Vc::cos(_req_phi[i] - ptphi);
                const auto isEdge = prelimChecks && (dist < cosDist);
#endif

                _stats_data.compares += CoordPacking;
                for(unsigned int j=0; j < isEdge.size(); ++j) {
                    const auto & neighbor = _req_ids[i*CoordPacking + j];

                    if (isEdge[j] && pt.id != neighbor) {
#ifndef NDEBUG
                        while(active_iter->id != neighbor)
                            ++active_iter;

                        if (_debugSuperParanoid && (pt.coshDistanceToPoincare(*active_iter) - _geometry.coshR) / _geometry.coshR > 1e-4) {
                            std::cerr << "pt:     " << pt << "\n"
                                         "req:    " << *active_iter << "\n"
                                         "req_ps: " << _req_phi_start[i][j] << "\n"
                                         "poi:    " << pt.coshDistanceToPoincare(*active_iter) << "\n"
                                         "hyp:    " << pt.coshDistanceToHyper(*active_iter) << "\n"
                                         "coshR:  " << _geometry.coshR << "\n"
                                         "j:      " << j << "\n"
                                         "dist:   " << dist << "\n"
                                         "thresh: " << threshR << "\n" <<
                            std::endl;
                            abort();
                        }
#endif

                        edgeCallback({pt.id, neighbor});
                        _stats_data.edges++;
                        noNeighbors += Statistics::enableNeighbors;
                    }

                    //if (Statistics::enableCandidates)
                    //    noCandidates += prelimChecks[j] && pt.id != neighbor;
                }
            }

            _stats_data.candidates.addPoint(noCandidates);
            _stats_data.neighbors.addPoint(noNeighbors);
#endif
        }

        // remove completed requests
        {
            auto thresh = std::min<Coord_b>(_phiRange.second, _generator.nextPointLB());

            if (p == _points.end() && thresh > _points.back().phi)
                thresh = _points.back().phi;

            if (p != _points.end() && thresh > p->phi)
                thresh = p->phi;

            if (_verbose) {
                std::cout << "Clean up with thresh = " << thresh << std::endl;
            }
            const auto end = std::remove_if(_active.begin(), _active.end(),
                                            [&] (const Request& a) {return a.range.second < thresh;});
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
        constexpr static bool enableHistograms { false };
        constexpr static bool enableActiveSizes{ enableHistograms && true };
        constexpr static bool enablePointSizes { enableHistograms && true };
        constexpr static bool enableCandidates { enableHistograms && true };
        constexpr static bool enableNeighbors  { enableHistograms && true };
        constexpr static bool enablePrelimCheck{ enableHistograms && true };

        Count batches{0};
        EdgeId edges{0};
        EdgeId compares{0};

        Histogram<enableActiveSizes> activeSizes;
        Histogram<enablePointSizes>  pointSizes;
        Histogram<enableCandidates>  candidates;
        Histogram<enableNeighbors>   neighbors;
        Histogram<enablePrelimCheck> prelimCheck;


        Statistics operator+(const Statistics& o) const {
            Statistics s;

            s.batches    = batches    + o.batches;
            s.activeSizes= activeSizes+ o.activeSizes;
            s.candidates = candidates + o.candidates;
            s.pointSizes = pointSizes + o.pointSizes;
            s.neighbors  = neighbors  + o.neighbors;
            s.edges      = edges      + o.edges;
            s.compares   = compares   + o.compares;
            s.prelimCheck= prelimCheck+ o.prelimCheck;

            return s;
        }
    };

    const Statistics& getStatistics() const {
        return _stats_data;
    }

    Coord nextRequestLB() const {
        return _generator.nextRequestLB();
    }

private:
    // Geometry and Parameter
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
    static constexpr bool _debugSuperParanoid {
#ifndef NDEBUG
        true
#else
        false
#endif
    };
    Statistics _stats_data;

    // Generator and AoS data
    PointGenerator _generator;

    std::vector<Request> _insertion_buffer;
    std::vector<Request> _active;
    bool _no_valid_requests{false};

    // Internal data structures (SoA)
    std::vector<Point> _points;

    std::vector<Node> _req_ids;
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_phi;
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_phi_start;
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_phi_stop;
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_cosh;

#ifdef POINCARE
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_poin_x;
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_poin_y;
    #ifdef LOG_TRANSFORM
        std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_poin_loginvr;
    #else
        std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_poin_invr;
    #endif

    #else
    std::vector<Coord_v, Vc::Allocator<Coord_v> > _req_invsinh;
#endif
    std::vector<Coord_m, Vc::Allocator<Coord_m> > _req_old;

    template <bool Endgame>
    unsigned int _AosToSoa(const Coord minThresh, const Coord maxThresh) {
        if (_active.empty())
            return 0;

        // transform AoS to SoA
        if (_req_ids.size() < _active.size()) {
            const size_t vectorSize = 2 * ((_active.size() + CoordPacking - 1) / CoordPacking);


            _req_ids.resize(vectorSize * CoordPacking);
            _req_phi.resize(vectorSize);
            _req_phi_start.resize(vectorSize);
            _req_phi_stop.resize(vectorSize);
            _req_cosh.resize(vectorSize);

            #ifdef POINCARE
            _req_poin_x.resize(vectorSize);
            _req_poin_y.resize(vectorSize);
            #ifdef LOG_TRANSFORM
                _req_poin_loginvr.resize(vectorSize);
            #else
                _req_poin_invr.resize(vectorSize);
            #endif
            #else
            _req_invsinh  .resize(vectorSize);
            #endif

            if (Endgame)
                _req_old.resize(vectorSize);
        }

        auto id_it = _req_ids.begin();

        unsigned int i=0;
        unsigned int j=0;

        for(const Request& req : _active) {
            if (req.range.second < minThresh || req.range.first > maxThresh)
                continue;

            *(id_it++) = req.id;
            _req_phi[i][j] = req.phi;
            _req_phi_start[i][j] = req.range.first;
            _req_phi_stop[i][j] = req.range.second;
#ifdef POINCARE
            _req_poin_x[i][j]    = req.poinX;
            _req_poin_y[i][j]    = req.poinY;
    #ifdef LOG_TRANSFORM
            _req_poin_loginvr[i][j] = req.poinLogInvLen;
    #else
            _req_poin_invr[i][j] = req.poinInvLen;
    #endif
#else
            _req_invsinh[i][j] = req.r.invsinh;
#endif
            _req_cosh[i][j] = req.r.cosh;

            if (Endgame) {
                _req_old[i][j] = req.old();
            }

            if (++j == CoordPacking) {
                j = 0;
                i++;
            }
        }

        if (j) {
            for (; j < CoordPacking; ++j) {
                _req_cosh[i][j] = std::numeric_limits<Coord>::max();
                *(id_it++) = std::numeric_limits<Node>::max();
            }
            ++i;
        }

        return i;
    }

    void _mergeInsertionBuffer(BandSegment& bandAbove);
    void _propagateRequests(const std::vector<Request>& reqs, BandSegment& bandAbove);


    std::vector<Request> _merge_helper;


#ifdef NDEBUG
    void _checkInvariants() const {}
#else
    void _checkInvariants() const;
#endif
};

#endif
