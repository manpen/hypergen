#pragma once
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <vector>
#include <memory>

#include "Geometry.hpp"
#include "Segment.hpp"

#include "ScopedTimer.hpp"

#include <omp.h>

class Generator {
public:
    Generator(Count n, Coord avgDeg, Coord alpha, Seed seed, uint32_t worker);
    Generator(const Generator& g) = default;

    
    
    template<typename EdgeCallback, typename PointCallback>
    void generate(EdgeCallback edgeCB, 
                  PointCallback pointCB
    ) {
        const unsigned int noBands = _bandLimits.size()-1;

        // GLOBAL PHASE
        {
            ScopedTimer timer("Global Generation Phase");
            // compute all global points and request
            _prepareGlobalPoints();
            for(unsigned int s=0; s < _segments.size(); ++s) {
                for(unsigned int b=0; b < _firstStreamingBand; ++b) {
                    for(const auto & pt : _segments[s]->getBand(b).getPoints()) {
                        pointCB(pt, s);
                    }
                }
            }
            
            // handle their edges
            for(unsigned int s=0; s < _segments.size(); ++s) { // TODO: in parallel?
                for(unsigned int b=0; b < _firstStreamingBand; ++b) {
                    auto& band =_segments[s]->getBand(b);

                    if (!band.getPoints().empty())
                        band.generateEdges<false>(
                            [&] (const Edge& e) {edgeCB(e, s);},
                            _segments[s]->getBand(b)    // Pass band itself in order to prevent propagation (which happened earlier)
                        );
                    band.getRequests().clear();
                }
            }
        }

        _dumpAllPointsAndRequests(_segments, "s");
        _dumpAllPointsAndRequests(_endgame_segments, "e");


        // Streaming Phase
        std::vector<Coord> maxPhis(_segments.size(), 0.0);

        {
            ScopedTimer timer("Main task");
//            omp_set_num_threads(_segments.size());

            #pragma omp parallel
            {
                #pragma omp for
                for (unsigned int s = 0; s < _segments.size(); ++s) {
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
                }

                _dumpAllPointsAndRequests(_segments, "e");

                #pragma omp for
                for (unsigned int oldSeg = 0; oldSeg < _segments.size(); ++oldSeg) {
                    // prepare endgame
                    const auto endgameSeg = (oldSeg + 1) % _segments.size();

                    for (unsigned int b = _firstStreamingBand; b < noBands; ++b) {
                        const auto &oldBand = _segments.at(oldSeg)->getBand(b);
                        auto &endgameBand = _endgame_segments.at(endgameSeg)->getBand(b);

                        const Coord maxPhi = endgameBand.prepareEndgame(oldBand);
                        if (maxPhi > maxPhis[endgameSeg])
                            maxPhis[endgameSeg] = maxPhi;

                    }

                    if (0 && _stats) {
                        std::cout << "maxPhi: " << (maxPhis[endgameSeg] - _endgame_segments[endgameSeg]->getPhiRange().first)
                                  << ", conservative: " << +_maxRepeatRange << std::endl;
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

        _reportEndStats();
    }
    
    const Geometry getGeometry() const { return _geometry; }
    
private:
    // statistics
    static constexpr bool _verbose {false};
    static constexpr bool _stats {true};
    
    const Geometry _geometry;
    const Count _noNodes;

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
