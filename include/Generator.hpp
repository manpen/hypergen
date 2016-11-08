#pragma once
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <vector>
#include <memory>

#include <omp.h>

#include "Geometry.hpp"
#include "Segment.hpp"
#include "RandomHelper.hpp"

#include <random>


class Generator {
public:
    Generator(Count n, Coord avgDeg, Coord alpha, Seed seed, uint32_t worker);
    
    
private:
    const Geometry _geometry;
    const Count _noNodes;
    
    const std::vector<Coord> _bandLimits;
    
    std::vector<std::unique_ptr<Segment>> _segments;
    
    std::vector<Coord> _computeBandLimits() const;
};


#endif
