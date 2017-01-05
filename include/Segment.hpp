#pragma once
#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <vector>
#include <memory>
#include <cassert>

#include "Geometry.hpp"
#include "BandSegment.hpp"
#include "Configuration.hpp"

class Segment {
public:
    Segment(Node firstNode, Count nodes,
        CoordInter phiRange, const Geometry& geometry,
        const Configuration& config,
        const std::vector<Coord>& limits,
        const unsigned int firstStreamingBand,
        Seed seed
    );

    BandSegment& getBand(unsigned int i) {
        assert(i < _bands.size());
        return *_bands[i];
    }

    const BandSegment& getBand(unsigned int i) const {
        assert(i < _bands.size());
        return *_bands[i];
    }

    BandSegment& getBandAbove(unsigned int i) {
        assert(i < _bands.size());
        if (++i == _bands.size())
            i = _bands.size() - 1;
        return *_bands[i];
    }

    const BandSegment& getBandAbove(unsigned int i) const {
        assert(i < _bands.size());
        if (++i == _bands.size())
            i = _bands.size() - 1;
        return *_bands[i];
    }


    template<bool Endgame, typename EdgeCallback, typename PointCallback>
    void advance(
        const unsigned int bandIdx,
        Coord threshold,
        const bool finalize,
        EdgeCallback edgeCB,
        PointCallback pointCB
    ) {
        auto& band = *_bands.at(bandIdx);

        // In the finalization step we need to recurse into all upper bands.
        // So even if this band is done, we cannot stop just yet
        if (band.done() && !finalize)
            return;

        unsigned int i=0;
        do {
            if (_verbose) {
                std::cout << std::string(bandIdx, ' ')
                          << "\x1B[31m"
                                "Advance (" << i++ << ") band " << bandIdx << " to threshold " << threshold << " in with final=" << finalize
                          << "\x1B[39m"
                          << std::endl;
            }

            if (i == 1000) {
                std::cerr << "Loop ?" << std::endl;
                abort();
            }

            // Generate points. This is necessary, if the last time
            // we went outside of our allowance
            if (1) {
                band.generatePoints();

                // invoke point callback
                if (!Endgame) {
                    for (const auto &pt : band.getPoints()) {
                        pointCB(pt);
                    }
                }
            }

            if (band.getPoints().empty() || band.getPoints().front().phi > threshold) {
                if (!finalize)
                    break;
            }

            // call even if there are no points, because there could be pending requests in the inbuf
            // that need to be propagated
            band.generateEdges<Endgame, false, EdgeCallback>(edgeCB, getBandAbove(bandIdx), threshold);

            if (bandIdx+1 < _bands.size()) {
                auto th = std::min<Coord_b>(std::min(band.nextRequestLB(), threshold), band.propagatedUntil());
                if (finalize && band.done(threshold)) {
                    th = threshold;
                    band.propagate<Endgame, false>(0.0, th, getBandAbove(bandIdx));
                }

                ASSERT_GE(band.getActive().nextInsertionPending(), th);
                advance<Endgame>(bandIdx+1, th, finalize, edgeCB, pointCB);
            }

        } while(finalize ? !band.done(threshold) : !band.done());

        //if (finalize)
        //    band.propagate<Endgame, false>(threshold, threshold, getBandAbove(bandIdx));
    }

    const CoordInter& getPhiRange() const {return _phiRange;}

private:
    static constexpr bool _verbose {VERBOSITY(true)};

    const Geometry _geometry;
    const CoordInter _phiRange;
    const unsigned int _firstStreamingBand;
    const Configuration& _config;

    std::vector<std::unique_ptr<BandSegment>> _bands;

};

#endif
