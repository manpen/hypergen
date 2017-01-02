#include "Segment.hpp"

#include "RandomHelper.hpp"
#include <cassert>
 
Segment::Segment(Node firstNode, Count nodes, 
        CoordInter phiRange, const Geometry& geometry,
        const std::vector<Coord>& limits,
        Seed seed
) 
    : _geometry(geometry)
    , _phiRange(phiRange)
{
    assert(!limits.empty());

    DefaultPrng randgen(seed);
    
    // distribute the requested number of points to bands
    // according to the radial distribution function
    std::vector<Count> pointsInBand;
    {
        std::vector<Coord> probs;
        Coord last_cdf = 0.0;
        for(unsigned int i=1; i < limits.size(); ++i) {
            assert(limits[i] <= _geometry.R);
            const Coord cdf = _geometry.radCdf(limits[i]);
            probs.push_back(cdf - last_cdf);
            last_cdf = cdf;
        }
        
        pointsInBand = RandomHelper::sampleMultinomial(nodes, probs, randgen);
    }
    
    // instantiate an instance for each band
    Node n0 = firstNode;
    for(unsigned int i=0; i < pointsInBand.size(); ++i) {
        _bands.push_back( std::make_unique<BandSegment>(
            n0, pointsInBand[i],
            phiRange, CoordInter{limits[i], limits[i+1]},
            _geometry,
            static_cast<uint32_t>(randgen())
        ));
        n0 += pointsInBand[i];
    }
}
