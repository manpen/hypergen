#pragma once
#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <vector>
#include <memory>
#include <cassert>

#include "Geometry.hpp"
#include "BandSegment.hpp"

class Segment {
public:    
    Segment(Node firstNode, Count nodes, 
            CoordInter phiRange, const Geometry& geometry,
            const std::vector<Coord>& limits,
            Seed seed
    ); 

    BandSegment& getBand(unsigned int i) {
        assert(i < _bands.size());
        return *_bands[i];
    }
    
    
private:
    const Geometry _geometry;
    const CoordInter _phiRange;
    
    std::vector<std::unique_ptr<BandSegment>> _bands;
    
};

#endif
