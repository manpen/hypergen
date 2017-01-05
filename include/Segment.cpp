#include "Segment.hpp"

#include "RandomHelper.hpp"
#include <cassert>
 
Segment::Segment(Node firstNode, Count nodes, 
        CoordInter phiRange, const Geometry& geometry,
        const Configuration& config,
        const std::vector<Coord>& limits,
        const unsigned int firstStreamingBand,
        Seed seed
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
    for(unsigned int i=0; i < pointsInBand.size(); ++i) {
        _bands.push_back( std::make_unique<BandSegment>(
            n0, pointsInBand[i],
            phiRange, CoordInter{limits[i], limits[i+1]},
            _geometry, _config,
            std::cosh(limits[firstStreamingBand]),
            static_cast<uint32_t>(randgen()),
            i
        ));
        n0 += pointsInBand[i];
    }
}
