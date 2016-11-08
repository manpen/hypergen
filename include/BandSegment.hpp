#pragma once
#ifndef BAND_SEGMENT_HPP
#define BAND_SEGMENT_HPP

#include "Geometry.hpp"
#include "PointGenerator.hpp"
#include <vector>
#include <algorithm>


class BandSegment {
public:
    BandSegment() = delete;
    
    BandSegment(Node firstNode, Count nodes, CoordInter phiRange, CoordInter radRange, const Geometry& geometry, Seed seed) 
        : _seed(seed)
        , _firstNode(firstNode)
        , _numberOfNodes(nodes)
        , _geometry(geometry)
        , _phiRange(phiRange)
        , _radRange(radRange)
        , _upperLimit(radRange.second) 
        , _batchSize(geometry.avgDeg * 2)
        , _generator(_newGenerator())
    {}
    
    template <typename EdgeCallback, typename PointCallback>
    void advance(const Coord threshold,
                 BandSegment& bandAbove,
                 EdgeCallback& edgeCallback,
                 PointCallback& pointCallback = [] (Point&) {}
    ) {
        _prepare_advance(threshold, bandAbove);
        
        if (_req_ids.empty())
            return;
        
        for(const Point& pt : _points) {
            const uint32_t lower = static_cast<uint32_t>(std::distance(_req_phi_start.cbegin(), 
                                       std::lower_bound(_req_phi_start.cbegin(), _req_phi_start.cend(), pt.phi)));
            const uint32_t upper = static_cast<uint32_t>(std::distance(_req_phi_stop.cbegin(), 
                                       std::upper_bound(_req_phi_stop.cbegin(), _req_phi_stop.cend(), pt.phi)));
            
            for(uint32_t i = lower; i < upper; ++i) {
                const bool ptIsSmaller =  std::tie(pt.r.cosh, pt.phi, pt.id) < std::tie(_req_cosh[i], _req_phi, _req_ids[i]);
                const Coord deltaPhi = M_PI - std::abs(M_PI - std::abs(_req_phi[i] - pt.phi));
                const Coord cosDist = (pt.r.cosh * _req_cosh[i] - _geometry.coshR) / (pt.r.sinh * _req_sinh[i]);
                
                if (ptIsSmaller && (cosDist < cos(deltaPhi))) {
                    edgeCallback(pt.id, _req_ids[i]);
                }
            }
        }
        
        for(const Point& pt : _points)
            pointCallback(pt);
        
        _points.clear();
    }
    
    Seed getSeed() const {return _seed;}
    Node getFirstNode() const {return _firstNode;}
    Node getNumberOfNodes() const {return _numberOfNodes;}
    
private:
    // Geometry and Parameter
    const Seed _seed;
    const Node _firstNode;
    const Count _numberOfNodes;
    const Geometry& _geometry;
    const CoordInter _phiRange;
    const CoordInter _radRange;
    const SinhCosh   _upperLimit;
    
    const Count _batchSize;
    
    PointGenerator _generator;
    
    std::vector<Request> _insertion_buffer;
    std::vector<Request> _active;
    
    
    // internal data structures
    std::vector<Point> _points;
    std::vector<Node > _req_ids;
    std::vector<Coord> _req_phi;     
    std::vector<Coord> _req_phi_start;
    std::vector<Coord> _req_phi_stop; 
    std::vector<Coord> _req_sinh;     
    std::vector<Coord> _req_cosh;         
    
    PointGenerator _newGenerator() const {
        return {_firstNode, _numberOfNodes, _phiRange, _radRange, _geometry, _seed};
    };
    
    void _prepare_advance(const Coord threshold, BandSegment& bandAbove);
    
};

#endif
