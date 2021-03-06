/**
 * @file
 * @brief BandSegment
 *
 * A band segment generates random points and requests, performs
 * distance calculation and propagates requests.
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

#include "ActiveManager.hpp"

#include "Configuration.hpp"
#include "RandomHelper.hpp"

#include <Vc/Vc>

//#define SKIP_DIST_COMP

class BandSegment {
public:
    struct Statistics {
        constexpr static bool enableHistograms { false };
        constexpr static bool enableActiveSizes{ true && enableHistograms};
        constexpr static bool enablePointSizes { true && enableHistograms };
        constexpr static bool enableCandidates { true && enableHistograms};
        constexpr static bool enableNeighbors  { !true && enableHistograms };
        constexpr static bool enablePrelimCheck{ !true && enableHistograms };

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



    BandSegment() = delete;
    BandSegment(Node firstNode, Count nodes,
                const CoordInter phiRange, CoordInter radRange,
                const Geometry& geometry,
                const Configuration& config,
                Coord_b streamingBound,
                Seed seed,
                unsigned int bandIdx,
                unsigned int firstStreamingBand,
                Statistics& stats
    );

    void enable() {
        _generator.reset(new PointGenerator(_firstNode, _numberOfNodes, _phiRange, _radRange, _geometry, DefaultPrng{_seed}));
    }

    // produce new points and request; delete exhausted requests from _active
    void generatePoints();

    void addRequest(const Request& req) {
        if (_verbose > 1)
            std::cout << "Add to band " << _bandIdx << " request " << req << std::endl;

        ASSERT_LE(req.range.first, req.range.second);

#ifndef NDEBUG
        if (req.r.cosh < _streamingBound) {
            //ASSERT_GE(req.range.first, _phiRange.first);
            ASSERT_LE(req.range.second, _phiRange.second);
        }

        if (_active.requestAlreadyPending(req)) {
            std::cout << "Double insertion of " << req << std::endl;
        }
#endif

        if (req.range.first < _processedUntil) {
            if (_verbose)
                std::cout << "Ignore " << req << " as before processing limit of " << _processedUntil << std::endl;
            return;
        }

        _active.addRequest(req);
    }

    template <bool Endgame, bool GlobalPhase>
    void propagate(Coord_b del, Coord_b ins, BandSegment& bandAbove) {
        if (!GlobalPhase && &bandAbove != this) {
            _active.update<Endgame>(del, ins, [&](const Request &req) {
                auto nreq = Request(req, _geometry, _upperLimit);

                if (_verbose > 2)
                    std::cout << "Propagate " << nreq << std::endl;

                if (likely(req.r.cosh >= _streamingBound)) {
                    // streaming request which cannot be split
                    bandAbove.addRequest(nreq);
                } else {
                    const auto& old_range = nreq.range;
                    auto range = nreq.range;

                    if (range.first < 0) {
                        range.first += 2*M_PI;
                        range.second += 2*M_PI;
                    }

                    ASSERT_LE(range.first, 2*M_PI);

                    if (range.second > range.first + 2*M_PI)
                        range.second = range.first + 2*M_PI;

                    const auto width = _phiRange.second - _phiRange.first;

                    if (range.second - range.first > 2*M_PI - width / 4) {
                        bandAbove.addRequest(Request(nreq, _phiRange));
                    } else {
                        if ((range.first <= _phiRange.first && range.second >= _phiRange.second) ||
                            (range.first-2*M_PI <= _phiRange.first && range.second-2*M_PI >= _phiRange.second)){
                            // segment completely enclosed
                            bandAbove.addRequest(Request(nreq, _phiRange));
                        } else if (range.first >= _phiRange.first && range.second <= _phiRange.second) {
                            // request completely enclosed
                            bandAbove.addRequest(Request(nreq, range));
                        } else {
                            if (range.first >= _phiRange.first && range.first <= _phiRange.second) {
                                bandAbove.addRequest(Request(nreq, {range.first, _phiRange.second}));
                            }

                            if (_phiRange.first <= range.second && range.second <= _phiRange.second) {
                                bandAbove.addRequest(Request(nreq, {_phiRange.first, range.second}));
                            } else if (_phiRange.first <= range.second - 2*M_PI && range.second - 2*M_PI <= _phiRange.second) {
                                bandAbove.addRequest(Request(nreq, {_phiRange.first, range.second-2*M_PI}));
                            }
                        }
                    }
                }
            });
        } else {
            _active.update<Endgame>(del, ins, [&](const Request &req) {});
        }
    }

    template <bool Endgame, bool GlobalPhase = false, typename EdgeCallback>
    void generateEdges(EdgeCallback edgeCallback, BandSegment& bandAbove, Coord threshold = 2*M_PI) {
        // Handle active updates
        const bool bandAboveExists = (&bandAbove != this);

        auto performUpdates = [&] (Coord_b del, Coord_b ins) {
            if (_verbose > 2)
                std::cout << "Band " << _bandIdx << " issue active update:" << std::endl;

            propagate<Endgame, GlobalPhase>(del, ins, bandAbove);
            _propagatedUntil = ins;
        };

        if (_points.empty() || _points.front().phi > threshold) {
            _no_valid_requests = true;
            return;
        }



        const Coord bandThreshold = _generator->nodesLeft()
                                    ? (_generator->nextRequestLB() - _maxDeltaPhi)
                                    : _phiRange.second;
        threshold = std::min(threshold, bandThreshold);

        if (_active.empty()) {
            std::cout << "Active empty" << std::endl;
            _no_valid_requests = true;
            return;
        }

        _no_valid_requests = false;

        _stats_data.pointSizes.addPoint( _points.size() );


        unsigned int pointIdx = 0;

        auto p = _points.cbegin();
        auto nextUpdate = p;
        for(; p != _points.cend(); ++p,++pointIdx) {
            const Point& pt = *p;

            if (unlikely(pt.phi > threshold))
                break;

            _processedUntil = pt.phi;

            _stats_data.pointSizes.addPoint( _points.size() );

#ifndef SKIP_DIST_COMP
            if (nextUpdate == p) {
                nextUpdate = _points.cbegin() +
                             std::min<size_t>(_points.size(), pointIdx+_config.activeUpdateInterval);

                ASSERT_GE((nextUpdate-1)->phi, pt.phi);

                performUpdates(pt.phi, (nextUpdate-1)->phi);
            }


            _stats_data.activeSizes.addPoint( _active.requestsPending() );
            _stats_data.candidates.addPoint( _active.size() );



            const Coord_v pt_phi(static_cast<Coord_b>(pt.phi));

            const Coord_v pt_poin_x(static_cast<Coord_b>(pt.poinX));
            const Coord_v pt_poin_y(static_cast<Coord_b>(pt.poinY));

        #ifdef LOG_TRANSFORM
            const Coord_v threshR = static_cast<Coord_b>(_geometry.logPoincareR - pt.poinLogInvLen);
            const Coord_v pt_poin_r(static_cast<Coord_b>(pt.poinLogInvLen));
        #else
            const Coord_v pt_poin_r(static_cast<Coord_b>(pt.poinInvLen));
            const Coord_v threshR = static_cast<Coord_b>(_geometry.poincareR / pt.poinInvLen);
        #endif

            bool pointFromLastSegment;
            if (Endgame)
                pointFromLastSegment = pt.old();

            unsigned int noCandidates = 0;
            unsigned int noNeighbors = 0;

            for(unsigned int i = 0; i < _active.end(); ++i) {
                auto prelimChecks = ((_active.req_poin_r(i) < pt_poin_r)
                                     || (_active.req_poin_r(i) == pt_poin_r && _active.req_phi(i) < pt_phi));
                
                if (Statistics::enablePrelimCheck)
                    _stats_data.prelimCheck.addPoint( prelimChecks.count() );

                if (prelimChecks.isEmpty())
                    continue;

                // computations
                const auto deltaX = _active.req_poin_x(i) - pt_poin_x;
                const auto deltaY = _active.req_poin_y(i) - pt_poin_y;

                #ifdef LOG_TRANSFORM
                const auto dist = deltaX*deltaX + deltaY*deltaY;
                const auto isEdge = prelimChecks && (dist < Vc::exp(threshR - _active.req_poin_r(i)));
                #else
                const auto dist = (deltaX*deltaX + deltaY*deltaY) * _active.req_poin_r(i);
                const auto isEdge = prelimChecks && (dist < threshR);
                #endif

                _stats_data.compares += CoordPacking;
                for(unsigned int j=0; j < isEdge.size(); ++j) {
                    const auto & neighbor = _active.req_id(i*CoordPacking + j);

                    if (_verbose > 2 && _active.req_poin_r(i)[j] <= _geometry.R) {
                        std::cout << (Endgame ? "End-" : "")
                                  << "Compare in band " << _bandIdx << " " << pt << " against " << (neighbor & Point::NODE_MASK)
                                  << " Prelim: " << prelimChecks[j] << " isEdge: " << isEdge[j]
                                  << " Point-Old: " << pt.old() << " Req-Old: " << bool(neighbor & Point::OLD_MASK)
                                  << std::endl;
                    }

                    bool cond = isEdge[j];

                    if (Endgame) {
                        cond &= (pt.id & Point::NODE_MASK) != (neighbor & Point::NODE_MASK);
                        cond &= pointFromLastSegment || (neighbor & Point::OLD_MASK);
                    } else {
                        cond &= pt.id != neighbor;
                    }


                    if (cond) {
                        ASSERT_NE(neighbor, std::numeric_limits<Node>::max());
                        if (Endgame)
                            edgeCallback({pt.id & Point::NODE_MASK, neighbor & Point::NODE_MASK});
                        else
                            edgeCallback({pt.id, neighbor});

                        _stats_data.edges++;
                        noNeighbors += Statistics::enableNeighbors;
                    }
                }
            }

            _stats_data.neighbors.addPoint(noNeighbors);
#endif
        }

#ifdef SKIP_DIST_COMP
        performUpdates( (p-1)->phi, (p-1)->phi );
#endif

        // remove complete points
        if (p == _points.cend()) {
            _points.clear();
            if (!_generator->nodesLeft())
                performUpdates(_phiRange.second, _phiRange.second);

        } else {
            _points.erase(_points.begin(), p);
        }
    }

    void copyGlobalState(const BandSegment& band, Coord deltaPhi);
    Coord_b prepareEndgame(BandSegment& band);

    template <typename EdgeCallback>
    void endgame(EdgeCallback edgeCallback) {

        ASSERT(!_active.empty());

        _checkInvariants();

        // compute the active request range, i.e. the request whos range ends the latest
        const Coord maxRange = _active.maxRange();
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
            _generator->generate(batchSize, _points);

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
    void sortPoints();

    ActiveManager& getActive() {return _active;}


    bool done() const {
        return _no_valid_requests && !_generator->nodesLeft();
    }

    bool done(const Coord thresh) const {
        return (_points.empty() || _points.front().phi > thresh) && (!_generator->nodesLeft() || _generator->nextRequestLB() > thresh);
    }


    Coord nextRequestLB() const {
        return _generator->nextRequestLB();
    }

    Coord_b propagatedUntil() const {
        return _propagatedUntil;
    }

    void clear();

private:
    // Geometry and Parameter
    const Node _firstNode;
    const Count _numberOfNodes;
    const Geometry& _geometry;
    const Configuration& _config;
    const CoordInter _phiRange;
    const CoordInter _radRange;
    const SinhCosh   _upperLimit;

    const Count _batchSize;
    const Coord_b _maxDeltaPhi;
    const unsigned int _bandIdx;

    // Statistics
    static constexpr unsigned int _verbose{VERBOSITY(3)};
    static constexpr bool _stats{true};

    Statistics& _stats_data;

    // Generator and AoS data
    std::unique_ptr<PointGenerator> _generator;

    bool _no_valid_requests{false};

    std::vector<Point> _points;

    ActiveManager _active;

    Coord_b _propagatedUntil  {-1.0};
    Coord_b _processedUntil {-1.0};

    const bool _coverFullCircle;
    const Coord_b _streamingBound;

    const Seed _seed;

#ifdef NDEBUG
    void _checkInvariants() const {}
#else
    void _checkInvariants() const;
#endif
};

#endif
