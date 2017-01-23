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
        _prepareGlobalPoints();

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
                    if (_config.bandLimits == Configuration::BandLimitType::BandLin) {
                        ASSERT_GE(_bandLimits[1], _geometry.R/2);

                        const auto& points = _segments[s]->getBand(0).getPoints();
                        for(auto it = points.cbegin(); it != points.cend(); ++it) {
                            for(auto nbr = it+1; nbr != points.cend(); ++nbr) {
                                edgeCB({it->id, nbr->id}, s);
                            }
                        }

                        for(auto ns=s+1; ns<_segments.size(); ++ns) {
                            for(const auto& nbr : _segments[ns]->getBand(0).getPoints()) {
                                for(const auto& pt : points) {
                                    edgeCB({pt.id, nbr.id}, s);
                                }
                            }
                        }

                        firstRequestBand = 1;
                    }

                    for (unsigned int b = firstRequestBand; b < _firstStreamingBand; ++b) {
                        auto &band = _segments[s]->getBand(b);

                        if (!band.getPoints().empty())
                            band.generateEdges<false, true>(
                                    [&](const Edge &e) { edgeCB(e, s); },
                                    (b < _firstStreamingBand)
                                    ? _segments[s]->getBandAbove(b)
                                    : _segments[s]->getBand(b) // the highest global band is not supposed to propagate requests
                            );

                        if (b < _firstStreamingBand)
                            band.propagate<false, true>(band.getPhiRange().second,
                                                        band.getPhiRange().second,
                                                        _segments[s]->getBand(b)
                            );
                    }
                }


                {
                    ScopedTimer timer(timer_main[s]);

                    auto &segment = *_segments[s];
                    auto &firstBand = segment.getBand(_firstStreamingBand);

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

                    // prepare endgame
                    const auto endgameSeg = (s + 1) % _segments.size();

                    for (unsigned int b = _firstStreamingBand; b < noBands; ++b) {
                        const auto &oldBand = _segments.at(s)->getBand(b);
                        auto &endgameBand = _endgame_segments.at(endgameSeg)->getBand(b);

                        const Coord maxPhi = endgameBand.prepareEndgame(oldBand);
                        if (maxPhi > maxPhis[endgameSeg])
                            maxPhis[endgameSeg] = maxPhi;

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
                                    [&](const Edge &e) { edgeCB(e, endgameSeg); },
                                    [&](const Point &p) { pointCB(p, endgameSeg); }
                            );
                        } while (!finalize);
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

    
    // helper functions
    std::vector<Coord> _computeBandLimits() const;
    unsigned int _computeFirstStreamingBand(double thresholdSize) const;

    // generation
    void _prepareGlobalPoints();
    void _reportEndStats() const;
    
    void _dumpAllPointsAndRequests(const std::vector<std::unique_ptr<Segment>>& segments, std::string key = "") const;
};

#endif
