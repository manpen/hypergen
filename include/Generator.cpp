/**
 * @file
 * @brief Implementation of Generator
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
    , _maxRepeatRange(2.0 * M_PI / (config.noSegments < 64 ? 8.0 : 1.5) / config.noSegments)
    , _bandLimits(_computeBandLimits())
    , _firstStreamingBand(_computeFirstStreamingBand(_maxRepeatRange))
{
    std::cout << "Max Repeat Range: " << _maxRepeatRange << "\n"
                 "Max Repeat Range / Seg: " << (_maxRepeatRange * 2.0 * M_PI / config.noSegments) 
              << std::endl;   

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
    const Coord segWidth = 2.0 * M_PI / config.noSegments;

    assert(segWidth > _maxRepeatRange);

    _randgens.resize(config.noSegments);
    _segments.resize(config.noSegments);
    _endgame_segments.resize(config.noSegments);

    std::vector<Node> firstNode(config.noSegments);
    firstNode[0] = _globalNodes;
    for(unsigned int i=0; i<config.noSegments-1; i++)
        firstNode[i+1] = firstNode[i] + nodesPerSegment[i];

    #pragma omp parallel for
    for(unsigned int s=0; s<config.noSegments; ++s) {
        Seed seed = static_cast<Seed>(firstNode[s] * 12345678 ^ 97712327 ^ config.seed);

        _segments[s].reset(new Segment(
            firstNode[s], nodesPerSegment[s],
            CoordInter{s*segWidth, (1+s)*segWidth},
            _geometry, _config, _bandLimits, _firstStreamingBand,
            seed, false
        ));

        _endgame_segments[s].reset(new Segment(
            firstNode[s], nodesPerSegment[s],
            CoordInter{s*segWidth, (1+s)*segWidth},
            _geometry, _config, _bandLimits, _firstStreamingBand,
            seed, true
        ));
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

    // Insert Points
    PointGenerator ptgen(0, _globalNodes,
                         {0, 2 * M_PI}, {0, _bandLimits[_firstStreamingBand]},
                         _geometry, _randgen);

    ptgen.generate(_globalNodes, [&] (Point& pt, const Request&) {
        if (pt.phi > 2 * M_PI)
            pt.phi -= 2 * M_PI;

        const unsigned int seg = pt.phi * invSegWidth;
        const unsigned int band = getBandIdx(pt);

        if (likely(band)) {
            repairAndInsertRequest(Request(pt, _geometry, pt.r), band);
        } else {
            repairAndInsertRequest(Request(pt, _geometry, bandLimits[1]), 1);
        }

        _segments[seg]->getBand(band).getPoints().push_back(pt);
    });

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
        Coord seriesSpacing = std::log(_config.bandLinFactor) /  (_geometry.alpha + 0.5);
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
                            [&] (const BandSegment::Statistics& s, const std::unique_ptr<Segment>& seg) {
                                return s + seg->getStatistics(b);});

        BandSegment::Statistics eg_stats;
        if (b >= _firstStreamingBand) {
            eg_stats =
                std::accumulate(_endgame_segments.cbegin(), _endgame_segments.cend(), BandSegment::Statistics{},
                                [&] (const BandSegment::Statistics& s, const std::unique_ptr<Segment>& seg) {
                                    return s + seg->getStatistics(b);});
        }

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
