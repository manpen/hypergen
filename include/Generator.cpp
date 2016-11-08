#include "Generator.hpp"
#include "RandomHelper.hpp"
#include <random>

Generator::Generator(Count n, Coord avgDeg, Coord alpha, Seed seed, uint32_t worker)
    : _geometry(n, avgDeg, alpha)
    , _noNodes(n)
    , _bandLimits(_computeBandLimits())
{
    DefaultPrng randgen(seed);
    auto nodesPerSegment = RandomHelper::sampleMultinomial(worker, worker, randgen);

    Node n0 = 0;
    const Coord segWidth = 2.0 * M_PI / worker;
    for(uint32_t i=0; i<worker; ++i) {
        _segments.emplace_back( new Segment {
            n0, nodesPerSegment[i],
            {i*segWidth, (1+i)*segWidth},
            _geometry, _bandLimits,
            static_cast<uint32_t>(randgen())
        });
        n0 += nodesPerSegment[i];
    }
}


std::vector<Coord> Generator::_computeBandLimits() const {
    // based on NetworKIT
    std::vector<Coord> bandRadius;
    bandRadius.push_back(0.0);

    constexpr Coord seriesRatio = 0.9;
    const uint32_t logn = std::ceil(std::log(_noNodes));
    const Coord a = _geometry.R * (1.0-seriesRatio) / (1.0 - std::pow(seriesRatio, logn));

    for (uint32_t i = 1; i < logn; i++) {
        Coord c_i = a * (1.0-pow(seriesRatio, i)) / (1.0-seriesRatio);
        bandRadius.push_back(c_i);
    }
    bandRadius.push_back(_geometry.R);
    return bandRadius;
}

