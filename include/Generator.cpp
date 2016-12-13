#include "Generator.hpp"
#include "RandomHelper.hpp"

#include <iostream>
#include <iomanip>

#include <algorithm>
#include <random>

Generator::Generator(Count n, Coord avgDeg, Coord alpha, Seed seed, uint32_t worker)
    : _geometry(n, avgDeg, alpha)
    , _noNodes(n)
    , _maxRepeatRange(2.0 * M_PI / 8 / worker)
    , _bandLimits(_computeBandLimits())
    , _firstStreamingBand(_computeFirstStreamingBand(_maxRepeatRange))
{
    DefaultPrng randgen(seed);
    auto nodesPerSegment = RandomHelper::sampleMultinomial(n, worker, randgen);
    
    // print out global stats
    if (_stats) { 
        std::cout << "Number of nodes requested: " << n << "\n"
                     "Average Degree: " << avgDeg << "\n"
                     "Number of Worker/Segments: " << worker << "\n"
                     "Number of bands: " << (_bandLimits.size()-1) << "\n"
                     "TargetRadius: " << _geometry.R
        << std::endl;
        
        std::cout << "Nodes per segment:";
        for(const auto& n : nodesPerSegment)
            std::cout << " " << n;
        std::cout << std::endl;
    }

    // instantiate workers
    Node n0 = 0;
    const Coord segWidth = 2.0 * M_PI / worker;

    assert(segWidth > _maxRepeatRange);

    for(unsigned int s=0; s<worker; ++s) {
        Seed seed = static_cast<Seed>(randgen());

        _segments.emplace_back( new Segment {
            n0, nodesPerSegment[s],
            {s*segWidth, (1+s)*segWidth},
            _geometry, _bandLimits,
            seed,
            false
        });

        _endgame_segments.emplace_back( new Segment {
            n0, nodesPerSegment[s],
            {s*segWidth, (1+s)*segWidth},
            _geometry, _bandLimits,
            seed,
            true
        });

        n0 += nodesPerSegment[s];
    }

    // print out statistics
    if (_stats) {
        std::cout
                << "MaxAllowRepeat: " << _maxRepeatRange << "\n"
                << "Band statistics [LowLim, UpLim, MaxDeltaPhi, PtsInBand, CumPts]"
        << std::endl;

        Count cumPoints = 0;
        for(unsigned int b=0; b < _bandLimits.size() - 1; ++b) {
            const auto pointsInBand =
                    std::accumulate(_segments.cbegin(), _segments.cend(), Count(0),
                                    [&] (Count s, const auto& seg) {return s + seg->getBand(b).getNumberOfNodes();});
            cumPoints += pointsInBand;

            std::cout << std::setw(4) << std::right << b
                      << " " << std::setw(10) << std::right << _bandLimits[b]
                      << " " << std::setw(10) << std::right << _bandLimits[b+1]
                      << " " << std::setw(10) << std::right << _geometry.deltaPhi(_bandLimits[b], _bandLimits[b])
                      << " " << std::setw(10) << std::right << pointsInBand
                      << " " << std::setw(10) << std::right << cumPoints
                      << " " << (b < _firstStreamingBand ? "global" : "stream")
                      << std::endl;
        }
    }
}

void Generator::_prepareGlobalPoints() {
    // compute all points in global section
    for(unsigned int sIdx = 0; sIdx < _segments.size(); ++sIdx) {
        for(unsigned int actBandIdx = 0; actBandIdx < _firstStreamingBand; ++actBandIdx) {
            _segments[sIdx]->getBand(actBandIdx).generateGlobalPoints();
        }
    }

    _dumpAllPointsAndRequests(_segments);

    const unsigned int noBands = _bandLimits.size() - 1;

    // compute sinh/cosh band limits
    std::vector<SinhCosh> bandLimits;
    for(unsigned int s=0; s <= noBands; ++s)
        bandLimits.emplace_back(SinhCosh(_bandLimits[s]));

    // distribute requests
    const Coord segWidth = (2.0*M_PI) / _segments.size();
    const Coord invSegWidth = 1.0 / segWidth; 
    
    for(unsigned int sIdx = 0; sIdx < _segments.size(); ++sIdx) {
        for(unsigned int actBandIdx = 0; actBandIdx <= _firstStreamingBand; ++actBandIdx) {
            auto& actBand = _segments[sIdx]->getBand(actBandIdx);
            
            // insert the request into all segments (partially) covered by range
            auto addRequest = [&] (const CoordInter range,
                                   const Request& req,
                                   const unsigned int bandIdx,
                                   bool maySkip = true)
            {
                const unsigned int startSeg = range.first * invSegWidth;
                const unsigned int stopSeg = std::min<unsigned int>(ceil(range.second * invSegWidth + 0.1), _segments.size());
                assert(startSeg <= stopSeg);
                
                //std::cout
                //    << " add " << req << " to band " << actBandIdx << " in segs " << startSeg << " to " << stopSeg << " originating from " << sIdx
                //    << " restricted to " << range.first << " to " << range.second
                //<< std::endl;
                
                for(unsigned int s = startSeg; s != stopSeg; ++s) {
                    // skip issuing segment, since the request is already in active requests
                    if (maySkip && (s == sIdx))
                        continue;
                    
                    auto & band = _segments.at(s)->getBand(bandIdx);

                    // insert point and reduce query range to the limits of this segment
                    band.getInsertionBuffer().push_back(req);
                    band.getInsertionBuffer().back().range = {
                        std::max<Coord>(range.first, band.getPhiRange().first),
                        std::min<Coord>(range.second, band.getPhiRange().second)
                    };
                }
            };
            
            // iterate over all requests and distribute them to all segments involved
            for(Request & req : actBand.getRequests()) {
                if (req.phi > 2*M_PI)
                    req.phi -= 2*M_PI;

                // distribute to other segments
                if (req.range.first >= 0.0 && req.range.second <= 2.0 * M_PI) {
                    addRequest(req.range, req, actBandIdx, true);
                } else {
                    const Coord minPolar = fmod(req.range.first + 2.0 * M_PI, 2.0 * M_PI);
                    const Coord maxPolar = fmod(req.range.second, 2.0 * M_PI);
                    addRequest({0, maxPolar}, req, actBandIdx, false);
                    addRequest({minPolar, 2 * M_PI}, req, actBandIdx, true);
                }

                // distribute to higher bands
                for(unsigned int insBand = actBandIdx+1; insBand < noBands; ++insBand) {
                    Request nReq(req, _geometry, bandLimits[insBand]);

                    // distribute to other bands
                    if (nReq.range.first >= 0.0 && nReq.range.second <= 2.0 * M_PI) {
                        addRequest(nReq.range, nReq, insBand, false);
                    } else {
                        const Coord minPolar = fmod(nReq.range.first + 2.0 * M_PI, 2.0 * M_PI);
                        const Coord maxPolar = fmod(nReq.range.second, 2.0 * M_PI);
                        addRequest({0, maxPolar}, nReq, insBand, false);
                        addRequest({minPolar, 2 * M_PI}, nReq, insBand, false);
                    }
                }

                // reduce range to lay within the boundaries of the current segment
                req.range = {
                    std::max<Coord>(req.range.first,  actBand.getPhiRange().first),
                    std::min<Coord>(req.range.second, actBand.getPhiRange().second)
                };
            }
            
            // iterate over all points and distribute them to the right segment
            for(Point & pt : actBand.getPoints()) {
                pt.phi = fmod(pt.phi, 2.0*M_PI);
                const unsigned int ptSeg = pt.phi * invSegWidth;
                if (ptSeg != sIdx) {
                    _segments.at(ptSeg)->getBand(actBandIdx).getPoints().push_back(pt);
                }
            }
        }
    }


    _dumpAllPointsAndRequests(_segments);

    // sort insertion buffers, clean points and sort them
    for(unsigned int s=0; s < _segments.size(); s++) {
        auto & seg = _segments[s];
        for(unsigned int actBandIdx = 0; actBandIdx < noBands; ++actBandIdx) {
            auto& band = seg->getBand(actBandIdx);

            // merge from insertions buffers
            auto& insBuf = band.getInsertionBuffer();
            std::sort(insBuf.begin(), insBuf.end());
            if (band.getRequests().empty()) {
                band.getRequests().swap(insBuf);
            } else {
                std::vector<Request> merged(band.getRequests().size() + insBuf.size());
                std::merge(insBuf.cbegin(), insBuf.cend(), band.getRequests().cbegin(), band.getRequests().cend(), merged.begin());
                band.getRequests().swap(merged);
                insBuf.clear();
            }

            // correct manipulated points
            auto& pts = band.getPoints();
            const auto end = std::remove_if(pts.begin(), pts.end(), [&seg] (const auto& pt) {
                return pt.phi < seg->getPhiRange().first || pt.phi > seg->getPhiRange().second;
            });
            pts.erase(end, pts.end());

            std::sort(pts.begin(), pts.end(), [] (const auto&a, const auto&b) {
                return a.phi < b.phi;
            });

            _endgame_segments[s]->getBand(actBandIdx).copyGlobalState(band, _maxRepeatRange);
        }
    }
}


std::vector<Coord> Generator::_computeBandLimits() const {
    // based on NetworKIT
    std::vector<Coord> bandRadius;
    bandRadius.push_back(0.0);

    constexpr Coord seriesRatio = 0.91;
    const uint32_t logn = std::ceil(std::log(_noNodes));
    const Coord a = _geometry.R * (1.0-seriesRatio) / (1.0 - std::pow(seriesRatio, logn));

    for (uint32_t i = 1; i < logn; i++) {
        Coord c_i = a * (1.0-std::pow(seriesRatio, i)) / (1.0-seriesRatio);
        bandRadius.push_back(c_i);
    }
    
    bandRadius.push_back(_geometry.R);
    return bandRadius;
}

unsigned int Generator::_computeFirstStreamingBand(const double thresholdSize) const {
    unsigned int firstStreamingBand = 0;
    Count globalPoints = 0;
    while(++firstStreamingBand < _bandLimits.size()) {
        if (2*_geometry.deltaPhi(_bandLimits[firstStreamingBand], _bandLimits[firstStreamingBand]) < thresholdSize)
            break;
    }

    assert(firstStreamingBand+1 < _bandLimits.size());
    return firstStreamingBand;
}

void Generator::_reportEndStats() const {
    if (!_stats)
        return;
  
    std::cout << "Band-wise statistics: " << std::endl;
    std::cout << std::setw(4) << std::right << "Band"
              << " " << "      #Comps"
              << " " << "      #Edges"
              << " " << "  #Comp/Edge" << " |"
              << " " << "      #Comps"
              << " " << "      #Edges"
              << " " << "  #Comp/Edge"
    << std::endl;

    BandSegment::Statistics tot_stats;
    std::vector<Histogram> activeSizes;
    std::vector<Histogram> pointSizes;

    for(unsigned int b=0; b < _bandLimits.size() - 1; ++b) {
        const BandSegment::Statistics stats =
            std::accumulate(_segments.cbegin(), _segments.cend(), BandSegment::Statistics{},
                            [&] (const auto& s, const auto& seg) {return s + seg->getBand(b).getStatistics();});

        const BandSegment::Statistics eg_stats =
                std::accumulate(_endgame_segments.cbegin(), _endgame_segments.cend(), BandSegment::Statistics{},
                                [&] (const auto& s, const auto& seg) {return s + seg->getBand(b).getStatistics();});

        tot_stats = tot_stats + stats + eg_stats;

        std::cout << std::setw(4) << std::right << b  
                << " " << std::setw(12) << std::right << stats.compares
                << " " << std::setw(12) << std::right << stats.edges
                << " " << std::setw(12) << std::right << (static_cast<double>(stats.compares) / stats.edges) << " |"
                << " " << std::setw(12) << std::right << eg_stats.compares
                << " " << std::setw(12) << std::right << eg_stats.edges
                << " " << std::setw(12) << std::right << (static_cast<double>(eg_stats.compares) / eg_stats.edges)
        << std::endl;

        activeSizes.push_back(stats.activeSizes + eg_stats.activeSizes);
        pointSizes.push_back(stats.pointSizes + eg_stats.pointSizes);
    }

    std::cout <<
        "Number of edges:    " << std::setw(12) << std::right << tot_stats.edges << "\n"
        "Number of compares: " << std::setw(12) << std::right << tot_stats.compares << "\n"
        "Compares per edge:  " << std::setw(12) << std::right << (static_cast<double>(tot_stats.compares) / tot_stats.edges) << "\n"
        "Average Degree:     " << std::setw(12) << std::right << (2.0 * tot_stats.edges / _noNodes)
    << std::endl;

    for(unsigned int b=0; b < activeSizes.size(); ++b) {
        std::cout << "Active size in band " << b << "\n";
        activeSizes[b].toStream(std::cout, "ACT-SZE-" + std::to_string(b));
        std::cout << std::endl;
    }

    std::cout << "Total points size:\n";
    tot_stats.pointSizes.toStream(std::cout, "TOT-PTS-SZE");
    std::cout << std::endl;

    std::cout << "Total active size:\n";
    tot_stats.activeSizes.toStream(std::cout, "TOT-ACT-SZE");
    std::cout << std::endl;

}

void Generator::_dumpAllPointsAndRequests(const std::vector<std::unique_ptr<Segment>>& segments, const std::string key) const {
    if (!_verbose)
        return;

    for(unsigned int s=0; s<segments.size(); ++s) {
        std::cout << key << " Segment " << s << std::endl;
        for(unsigned int b=0; b <  _bandLimits.size()-1; ++b) {
            std::cout << key << " | Band " << b << std::endl;
            for(const auto& p : segments[s]->getBand(b).getPoints())
                std::cout << key << " | | " << p << std::endl;
            
            
            const auto& r1 = segments[s]->getBand(b).getRequests();
            const auto& r2 = segments[s]->getBand(b).getInsertionBuffer();

            auto i1 = r1.cbegin();
            auto i2 = r2.cbegin();

            while(i1 != r1.cend() || i2 != r2.cend()) {
                if (i1 == r1.cend() || (i2 != r2.cend() && *i2 < *i1)) {
                    std::cout << key << " | | InBu: " << *i2 << std::endl;
                    i2++;
                } else {
                    std::cout << key << " | | Reqs: " << *i1 << std::endl;
                    i1++;
                }
            }
        }
    }
}


