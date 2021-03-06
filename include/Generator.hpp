/**
 * @file
 * @brief Generator
 *
 * Setup geometry and perform global phase
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
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <vector>
#include <memory>

#include "RandomHelper.hpp"
#include "Geometry.hpp"
#include "Configuration.hpp"
#include "Segment.hpp"

#include "ScopedTimer.hpp"

#include <omp.h>

class Generator {
public:
    Generator(const Configuration& config);
    Generator(const Generator& g) = default;
    
    template<typename EdgeCallback, typename PointCallback>
    void generate(EdgeCallback edgeCB, 
                  PointCallback pointCB
    ) {
        const unsigned int noBands = _bandLimits.size()-1;
        std::vector<Coord> maxPhis(_segments.size(), 0.0);

        std::vector<double> timer_global(_segments.size());
        std::vector<double> timer_main(_segments.size());
        std::vector<double> timer_endgame(_segments.size());



        // GLOBAL PHASE
        {
            ScopedTimer timer("Preparation");
            _prepareGlobalPoints();
        }

        #pragma omp parallel
        {
            #pragma omp for
            for(unsigned int s=0; s < _segments.size(); ++s) {
                for (unsigned int b = 0; b < _firstStreamingBand; ++b) {
                    for (const auto &pt : _segments[s]->getBand(b).getPoints()) {
                        pointCB(pt, s);
                    }
                }

                {
                    ScopedTimer timer(timer_global[s]);

                    unsigned int firstRequestBand = 0;
                    ASSERT_GE(_bandLimits[1], _geometry.R / 2);

                    const auto &points = _segments[s]->getBand(0).getPoints();
                    for (auto it = points.cbegin(); it != points.cend(); ++it) {
                        for (auto nbr = it + 1; nbr != points.cend(); ++nbr) {
                            edgeCB({it->id, nbr->id}, s);
                        }
                    }

                    for (auto ns = s + 1; ns < _segments.size(); ++ns) {
                        for (const auto &nbr : _segments[ns]->getBand(0).getPoints()) {
                            for (const auto &pt : points) {
                                edgeCB({pt.id, nbr.id}, s);
                            }
                        }
                    }

                    firstRequestBand = 1;

                    for (unsigned int b = firstRequestBand; b < _firstStreamingBand; ++b) {
                        auto &band = _segments[s]->getBand(b);

                        {
                            auto &above = _segments[s]->getBand(b + 1);
                            above.getActive().copyFromBelow(band.getActive(), _geometry, above.getRadRange().first, band.getPhiRange());
                        }

                        if (!band.getPoints().empty()) {
                            band.enable();
                            band.sortPoints();

                            band.generateEdges<false, true>(
                                    [&](const Edge &e) { edgeCB(e, s); },
                                    _segments[s]->getBand(b) // the highest global band is not supposed to propagate requests
                            );
                        }

                        _segments.at(s)->releaseBand(b);
                    }
                }
            }

            #pragma omp for
            for(unsigned int s=0; s < _segments.size(); ++s) {
                _endgame_segments[s]->getBand(_firstStreamingBand).copyGlobalState(
                        _segments[s]->getBand(_firstStreamingBand), _maxRepeatRange);
            }

            #pragma omp for
            for(unsigned int s=0; s < _segments.size(); ++s) {
                {
                    ScopedTimer timer(timer_main[s]);

                    auto &segment = *_segments[s];
                    auto &firstBand = segment.getBand(_firstStreamingBand);

                    for(unsigned int b=_firstStreamingBand; b<noBands; ++b)
                        segment.getBand(b).enable();

                    // the main job: recursively merge all bands
                    {
                        bool finalize = true;
                        do {
                            finalize = !finalize;
                            segment.advance<false>(
                                    _firstStreamingBand,
                                    firstBand.getPhiRange().second,
                                    finalize,
                                    [&](const Edge &e) { edgeCB(e, s); },
                                    [&](const Point &p) { pointCB(p, s); }
                            );
                        } while (!finalize);
                    }

                    // clean up requests
                    auto phiEnd = segment.getPhiRange().second;
                    for (unsigned int b = _firstStreamingBand; b < noBands; ++b) {
                        segment.getBand(b).getActive().update<false>(phiEnd, phiEnd, [](const Request &) {});
                    }
                }

                {
                    ScopedTimer timer(timer_endgame[s]);
                    const auto endgameSeg = (s + 1) % _segments.size();

                    // prepare endgame

                    for (unsigned int b = _firstStreamingBand; b < noBands; ++b) {
                        auto &oldBand = _segments.at(s)->getBand(b);
                        auto &endgameBand = _endgame_segments.at(endgameSeg)->getBand(b);

                        const Coord maxPhi = endgameBand.prepareEndgame(oldBand);
                        if (maxPhi > maxPhis[endgameSeg])
                            maxPhis[endgameSeg] = maxPhi;

                        //oldBand.clear();
                        _segments.at(s)->releaseBand(b);
                        endgameBand.enable();
                    }

                    // report and assert replay width
                    {
                        const auto width = (maxPhis[endgameSeg] - _endgame_segments[endgameSeg]->getPhiRange().first);
                        if (_verbose) {
                            std::cout << "Replay-width of endgameSeg id " << endgameSeg << ": " << width << ". "
                                    "Conservative: " << _maxRepeatRange << std::endl;
                        }
                        #ifndef NDEBUG
                        for (unsigned i = _firstStreamingBand; i < noBands; ++i) {
                            std::cout << "Band " << i << "\n" << _endgame_segments[endgameSeg]->getBand(i).getActive() << std::endl;
                            for (const auto &pt: _endgame_segments[endgameSeg]->getBand(i).getPoints())
                                std::cout << pt << std::endl;
                        }
                        #endif
                        ASSERT_LE(width, _maxRepeatRange);
                    }

                    // execute endgame
                    {
                        bool finalize = true;

                        do {
                            finalize = !finalize;
                            _endgame_segments[endgameSeg]->advance<true>(
                                    _firstStreamingBand,
                                    maxPhis[endgameSeg],
                                    finalize,
                                    [&](const Edge &e) { edgeCB(e, s); },
                                    [&](const Point &p) { abort(); }
                            );
                        } while (!finalize);
                    }

                    for (unsigned int b = _firstStreamingBand; b < noBands; ++b) {
                        _endgame_segments[endgameSeg]->releaseBand(b);
                    }
                }
            }
        }

        if (_stats) {
            auto reportTimerStats = [] (const std::string& key, std::vector<double> timers) {
                if (timers.size() == 1) {
                    std::cout << timers.front() << "ms";
                    return;
                }

                double avg = std::accumulate(timers.cbegin(), timers.cend(), 0.0) / timers.size();

                std::sort(timers.begin(), timers.end());
                std::cout << key << " min: " << timers.front() << "ms\n"
                          << key << " max: " << timers.back() << "ms\n"
                          << key << " avg: " << avg << "ms\n";
            };

            reportTimerStats("Global  Phase Timer", timer_global);
            reportTimerStats("Main    Phase Timer", timer_main);
            reportTimerStats("Endgame Phase Timer", timer_endgame);
        }

        _reportEndStats();
    }
    
    const Geometry getGeometry() const { return _geometry; }
    
private:
    // statistics
    static constexpr bool _verbose {VERBOSITY(false)};
    static constexpr bool _stats {true};

    const Configuration& _config;
    const Geometry _geometry;
    const Count _noNodes;
    Node _globalNodes;
    DefaultPrng _randgen;

    // geometry
    const Coord _maxRepeatRange;
    const std::vector<Coord> _bandLimits;
    unsigned int _firstStreamingBand;
    
    std::vector<std::unique_ptr<Segment>> _segments;
    std::vector<std::unique_ptr<Segment>> _endgame_segments;

    std::vector<std::unique_ptr<DefaultPrng>> _randgens;

    
    // helper functions
    std::vector<Coord> _computeBandLimits() const;
    unsigned int _computeFirstStreamingBand(double thresholdSize) const;

    // generation
    void _prepareGlobalPoints();
    void _reportEndStats() const;
};

#endif
