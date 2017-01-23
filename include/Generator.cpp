#include "Generator.hpp"
#include "RandomHelper.hpp"

#include <iostream>
#include <iomanip>

#include <algorithm>
#include <random>
#include "PointGenerator.hpp"

Generator::Generator(const Configuration& config)
    : _config(config)
    , _geometry(config.nodes, config.avgDegree, config.alpha, config.R)
    , _noNodes(config.nodes)
    , _randgen(config.seed)
    , _maxRepeatRange(2.0 * M_PI / 4 / config.noSegments)
    , _bandLimits(_computeBandLimits())
    , _firstStreamingBand(_computeFirstStreamingBand(_maxRepeatRange))
{
    // compute number of nodes in global section
    {
        const auto prob = _geometry.radCdf(_bandLimits[_firstStreamingBand]);
        std::binomial_distribution<Node> distr(_noNodes, prob);
        _globalNodes = distr(_randgen);
    }

    auto nodesPerSegment = RandomHelper::sampleMultinomial(config.nodes - _globalNodes,
                                                           config.noSegments,
                                                           _randgen);

    // print out global stats
    if (_stats) { 
        std::cout << "Number of bands: " << (_bandLimits.size()-1) << "\n"
                     "TargetRadius: " << _geometry.R
        << std::endl;
        
        std::cout << "Streaming nodes per segment:";
        for(const auto& n : nodesPerSegment)
            std::cout << " " << n;
        std::cout << std::endl;
    }

    // instantiate workers
    Node n0 = _globalNodes;
    const Coord segWidth = 2.0 * M_PI / config.noSegments;

    assert(segWidth > _maxRepeatRange);

    for(unsigned int s=0; s<config.noSegments; ++s) {
        Seed seed = static_cast<Seed>(_randgen());

        _segments.push_back( std::make_unique<Segment>(
            n0, nodesPerSegment[s],
            CoordInter{s*segWidth, (1+s)*segWidth},
            _geometry, _config, _bandLimits, _firstStreamingBand,
            seed
        ));

        _endgame_segments.push_back( std::make_unique<Segment>(
            n0, nodesPerSegment[s],
            CoordInter{s*segWidth, (1+s)*segWidth},
            _geometry, _config, _bandLimits, _firstStreamingBand,
            seed
        ));

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
                      << " " << std::setw(10) << std::right <<  pointsInBand
                      << " " << std::setw(10) << std::right << cumPoints
                      << " " << (b < _firstStreamingBand ? "global" : "stream")
                      << std::endl;
        }
    }
}

void Generator::_prepareGlobalPoints() {
    // compute all points in global section
    std::vector<Request> requests;
    std::vector<Point> points;
    {
        PointGenerator ptgen(0, _globalNodes,
                             {0, 2 * M_PI}, {0, _bandLimits[_firstStreamingBand]},
                             _geometry, _randgen());

        ptgen.generate(_globalNodes, points, requests);
    }



    // distribute points
    const Coord segWidth = (2.0*M_PI) / _segments.size();
    const Coord invSegWidth = 1.0 / segWidth;

    // compute sinh/cosh band limits
    std::vector<SinhCosh> bandLimits;
    const auto noBands = _bandLimits.size() - 1;
    for(unsigned int s=0; s <= noBands; ++s)
        bandLimits.emplace_back(SinhCosh(_bandLimits[s]));
    const auto limitEnd = bandLimits.cbegin() + _firstStreamingBand + 1;

    //TODO: If we keep linear bands we can replace this by simple division
    auto getBandIdx = [&] (const Point& pt) {
        unsigned int band = std::distance(bandLimits.cbegin(),
            std::lower_bound(bandLimits.cbegin(), limitEnd, pt.r.cosh,
                             [] (const SinhCosh& l, const Coord& c) {return l.cosh < c;}));

        ASSERT_GT(band, 0);
        ASSERT_LE(band, _firstStreamingBand);

        return band-1;
    };


    // Insert Points
    {
        for (auto &pt : points) {
            if (pt.phi > 2 * M_PI)
                pt.phi -= 2 * M_PI;

            const unsigned int seg = pt.phi * invSegWidth;
            const unsigned int band = getBandIdx(pt);

            _segments[seg]->getBand(band).getPoints().push_back(pt);
        }

        for (unsigned int sIdx = 0; sIdx < _segments.size(); ++sIdx) {
            for (unsigned int band = 0; band < _firstStreamingBand; ++band)
                _segments[sIdx]->getBand(band).sortPoints();
        }
    }

    // Insert Requests
    {
        // insert the request into all segments (partially) covered by range
        auto addRequest = [&](const CoordInter range,
                              const Request &req,
                              const unsigned int bandIdx) {
            const unsigned int startSeg = range.first * invSegWidth;
            const unsigned int stopSeg = std::min<unsigned int>(ceil(range.second * invSegWidth), _segments.size());
            assert(startSeg <= stopSeg);

            for (unsigned int s = startSeg; s != stopSeg; ++s) {
                auto &band = _segments.at(s)->getBand(bandIdx);
                band.addRequest(
                        Request(req, {
                                std::max<Coord>(range.first, band.getPhiRange().first),
                                std::min<Coord>(range.second, band.getPhiRange().second)
                        })
                );
            }
        };

        auto repairAndInsertRequest = [&](Request req, const unsigned bandIdx) {
            // distribute to other segments
            if (req.range.first >= 0.0 && req.range.second <= 2.0 * M_PI) {
                addRequest(req.range, req, bandIdx);
            } else if (req.range.second - req.range.first >= 2 * M_PI - 1e-3) {
                addRequest({0, 2 * M_PI}, req, bandIdx);
            } else if (req.range.first <= 0) {
                addRequest({0, req.range.second}, req, bandIdx);
                addRequest({2 * M_PI + req.range.first, 2 * M_PI}, req, bandIdx);
            } else {
                ASSERT_GT(req.range.second, 2 * M_PI);
                addRequest({0, req.range.second - 2 * M_PI}, req, bandIdx);
                addRequest({req.range.first, 2 * M_PI}, req, bandIdx);
            }
        };

        for (auto &req : requests) {
            // iterate over all requests and distribute them to all segments involved
            if (req.phi > 2 * M_PI)
                req.phi -= 2 * M_PI;

            const unsigned int band = getBandIdx(req);
            repairAndInsertRequest(req, band);

            for (unsigned int b = band + 1; b <= _firstStreamingBand; ++b)
                repairAndInsertRequest(Request(req, _geometry, bandLimits[b]), b);
        }
    }

    // copy state of first streaming band into endgame
    for(unsigned int s=0; s < _segments.size(); s++) {
        auto & seg = _segments[s];
        auto& band = seg->getBand(_firstStreamingBand);
        _endgame_segments[s]->getBand(_firstStreamingBand).copyGlobalState(band, _maxRepeatRange);
    }

    if (_stats) {
        Node cumPoints = 0;
        for(unsigned int b=0; b<_firstStreamingBand; b++) {
            const auto points = std::accumulate(_segments.cbegin(), _segments.cend(), 0,
                [b] (const Node s, const std::unique_ptr<Segment>& seg) {return s+ seg->getBand(b).getPoints().size();}
            );
            cumPoints += points;
            std::cout << std::setw(4) << std::right << b
                      << " " << std::setw(10) << std::right << points
                      << " " << std::setw(10) << std::right << cumPoints
                      << std::endl;
        }

    }
}


std::vector<Coord> Generator::_computeBandLimits() const {
    std::vector<Coord> bandRadius;

    if (_config.bandLimits == Configuration::BandLimitType::BandLin) {
        Coord seriesSpacing = _config.bandLinFactor;
        const unsigned int len = std::ceil(_geometry.R / 2 / seriesSpacing);
        seriesSpacing = _geometry.R / 2 / len;

        std::cout << "Use " << len << " bands with a spacing of " << seriesSpacing << std::endl;

        constexpr auto baseBands = 1;

        bandRadius.reserve(baseBands + len + 1);

        for(unsigned int k=0; k <= baseBands; k++)
            bandRadius.push_back(_geometry.R / 2 / baseBands * k);

        for(unsigned int k=1; k <= len; k++)
            bandRadius.push_back(_geometry.R / 2 + seriesSpacing * k);

        bandRadius.back() = _geometry.R;

    } else {
        // based on NetworKIT
        bandRadius.push_back(0.0);

        constexpr Coord seriesRatio = 0.91;
        const uint32_t logn = ceil(_geometry.R) * _config.bandExpFactor;
        const Coord a = _geometry.R * (1.0 - seriesRatio) / (1.0 - std::pow(seriesRatio, logn));

        for (uint32_t i = 1; i < logn; i++) {
            Coord c_i = a * (1.0 - std::pow(seriesRatio, i)) / (1.0 - seriesRatio);
            bandRadius.push_back(c_i);
        }

        bandRadius.push_back(_geometry.R);
    }

    return bandRadius;
}

unsigned int Generator::_computeFirstStreamingBand(const double thresholdSize) const {
    unsigned int firstStreamingBand = 0;
    Count globalPoints = 0;
    while(++firstStreamingBand < _bandLimits.size()) {
        if (2*_geometry.deltaPhi(_bandLimits[firstStreamingBand], _bandLimits[firstStreamingBand]) < thresholdSize)
            break;
    }

    if (firstStreamingBand+1 >= _bandLimits.size()) {
        std::cerr << "No streaming band; increase graph size" << std::endl;
        abort();
    }
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
    std::vector<BandSegment::Statistics> bandwiseStats;

    for(unsigned int b=0; b < _bandLimits.size() - 1; ++b) {
        const BandSegment::Statistics stats =
            std::accumulate(_segments.cbegin(), _segments.cend(), BandSegment::Statistics{},
                            [&] (const BandSegment::Statistics& s, const std::unique_ptr<Segment>& seg) {return s + seg->getBand(b).getStatistics();});

        const BandSegment::Statistics eg_stats =
                std::accumulate(_endgame_segments.cbegin(), _endgame_segments.cend(), BandSegment::Statistics{},
                                [&] (const BandSegment::Statistics& s, const std::unique_ptr<Segment>& seg) {return s + seg->getBand(b).getStatistics();});

        tot_stats = tot_stats + stats + eg_stats;

        std::cout << std::setw(4) << std::right << b  
                << " " << std::setw(12) << std::right << stats.compares
                << " " << std::setw(12) << std::right << stats.edges
                << " " << std::setw(12) << std::right << (static_cast<double>(stats.compares) / stats.edges) << " |"
                << " " << std::setw(12) << std::right << eg_stats.compares
                << " " << std::setw(12) << std::right << eg_stats.edges
                << " " << std::setw(12) << std::right << (static_cast<double>(eg_stats.compares) / eg_stats.edges)
        << std::endl;

        bandwiseStats.push_back(stats);
    }

    std::cout <<
        "Number of edges:    " << std::setw(12) << std::right << tot_stats.edges << "\n"
        "Number of compares: " << std::setw(12) << std::right << tot_stats.compares << "\n"
        "Compares per edge:  " << std::setw(12) << std::right << (static_cast<double>(tot_stats.compares) / tot_stats.edges) << "\n"
        "Average Degree:     " << std::setw(12) << std::right << (2.0 * tot_stats.edges / _noNodes)
    << std::endl;

// Dump Histograms
    for(unsigned int b=0; b < bandwiseStats.size(); ++b) {
        bandwiseStats[b].activeSizes.toStream(std::cout, "HIST-ACT-" + std::to_string(b));
        bandwiseStats[b].candidates.toStream(std::cout, "HIST-CND-" + std::to_string(b));
        bandwiseStats[b].pointSizes.toStream(std::cout, "HIST-PTS-" + std::to_string(b));
        bandwiseStats[b].neighbors.toStream(std::cout, "HIST-NBR-" + std::to_string(b));
    }

    tot_stats.activeSizes.toStream(std::cout, "HIST-TOT-ACT");
    tot_stats.candidates.toStream(std::cout, "HIST-TOT-CND");
    tot_stats.pointSizes.toStream(std::cout, "HIST-TOT-PTS");
    tot_stats.neighbors.toStream(std::cout, "HIST-TOT-NBR");
    tot_stats.prelimCheck.toStream(std::cout, "HIST-PRE-CHK");
}

void Generator::_dumpAllPointsAndRequests(const std::vector<std::unique_ptr<Segment>>& segments, const std::string key) const {
    if (!_verbose)
        return;
/*
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
*/
}


