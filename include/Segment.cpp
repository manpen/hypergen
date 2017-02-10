/**
 * @file
 * @brief Implementation of Segment
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
#include "Segment.hpp"

#include <cassert>
 
Segment::Segment(Node firstNode, Count nodes, 
        CoordInter phiRange, const Geometry& geometry,
        const Configuration& config,
        const std::vector<Coord>& limits,
        const unsigned int firstStreamingBand,
        Seed seed,
        bool streamingOnly
) 
    : _geometry(geometry)
    , _phiRange(phiRange)
    , _firstStreamingBand(firstStreamingBand)
    , _config(config)
{
    assert(!limits.empty());
    DefaultPrng randgen(seed);

    // distribute the requested number of points to bands
    // according to the radial distribution function
    std::vector<Count> pointsInBand;
    {
        std::vector<Coord> probs(firstStreamingBand, 0.0);

        Coord last_cdf = _geometry.radCdf(limits[firstStreamingBand-1]);
        Coord norm = 1.0 / (1.0 - last_cdf);

        for(unsigned int i=firstStreamingBand+1; i < limits.size(); ++i) {
            assert(limits[i] <= _geometry.R);
            const Coord cdf = _geometry.radCdf(limits[i]);
            probs.push_back( (cdf - last_cdf) * norm ) ;
            last_cdf = cdf;
        }

        pointsInBand = RandomHelper::sampleMultinomial(nodes, probs, randgen);

        ASSERT_EQ(pointsInBand.size(), limits.size()-1);
        ASSERT_EQ(std::accumulate(pointsInBand.cbegin(), pointsInBand.cend(), 0), nodes);
    }

    // instantiate an instance for each band
    Node n0 = firstNode;

    _stats.resize(pointsInBand.size());

    for(unsigned int i=0; i < pointsInBand.size(); ++i) {
        const Seed mSeed = randgen();

        if (streamingOnly && i < _firstStreamingBand) {
            _bands.push_back( std::unique_ptr<BandSegment>() );
        } else {
            _bands.push_back(std::make_unique<BandSegment>(
                    n0, pointsInBand[i],
                    phiRange, CoordInter{limits[i], limits[i + 1]},
                    _geometry, _config,
                    std::cosh(limits[firstStreamingBand]),
                    mSeed,
                    i, firstStreamingBand,
                    _stats[i]
            ));
        }
        n0 += pointsInBand[i];
    }
}
