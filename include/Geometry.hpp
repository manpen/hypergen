#pragma once
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <cmath>
#include <cstdint>
#include <utility>
#include <tuple>

#include <iostream>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)       __builtin_expect((x),0)


using Coord = double;
using Node = uint64_t;
using Count = Node;
using Seed = uint32_t;

using EdgeId = uint64_t;
using Edge = std::pair<Node, Node>;

template <typename T>
using Interval = std::pair<T, T>;

using CoordInter = Interval<Coord>;

struct SinhCosh {
    Coord sinh;
    Coord cosh;
    
    SinhCosh() {}
    
    SinhCosh(const Coord r)
    {
        const Coord e  = std::exp(r);
        const Coord ie = 1.0 / e;

        sinh = 0.5 * (e - ie);
        cosh = 0.5 * (e + ie);     
    }
    
    SinhCosh(const Coord sinh, const Coord cosh)
        : sinh(sinh), cosh(cosh)
    {}
        
};

struct Point {
    Node id;
    Coord phi;
    SinhCosh r;

    bool _old{false};

    
    Point() {}
    
    Point(const Node id, const Coord phi, const SinhCosh r)
        : id(id), phi(phi), r(r)
    {}

    void setOld() {
        _old = true;
    }

    bool old() const {
        return _old;
    }
    
    bool operator<(const Point& o) const {
        return std::tie(r.cosh, phi) < std::tie(o.r.cosh, o.phi);
    }
    
    friend std::ostream& operator <<(std::ostream& stream, const Point& o) {
        return stream << "Point(id:" << o.id << ", phi: " << o.phi << ", r: " << std::acosh(o.r.cosh) << ", old:" << o.old() << ")";
    }

    Coord distanceTo(const Point& o) const {
        return std::acosh(r.cosh*o.r.cosh - r.sinh*o.r.sinh*cos(o.phi - phi));
    }
};

struct Geometry {
    const Coord alpha;
    const Coord avgDeg;
    const Coord R;
    const Coord coshR;
    const Coord coshAlphaR;
    
    Geometry() = delete;
    
    Geometry(Coord alpha, Coord avgDeg, Coord R) 
        : alpha(alpha)
        , avgDeg(avgDeg)
        , R(R)
        , coshR(std::cosh(R))
        , coshAlphaR(std::cosh(alpha*R))
    {}    
    
    Geometry(Count nodes, Coord avgDeg, Coord alpha)
        : Geometry(alpha, avgDeg, _computeTargetRadius(nodes, avgDeg, alpha))
    {}
    
    Coord radDensity(Coord r) const {
        return alpha * std::sinh(r*alpha) / (coshAlphaR - 1.0);
    }
    
    Coord radCdf(Coord r) const {
        return (std::cosh(alpha * r) - 1.0) / (coshAlphaR - 1.0);
    }
    
    Coord radInvCdf(Coord x) const {
        return std::acosh(x * coshAlphaR - x + 1) / alpha;
    }

    Coord deltaPhi(const SinhCosh pt, const SinhCosh band) const {
        Coord deltaPhiCos = (pt.cosh * band.cosh - coshR) / (pt.sinh * band.sinh);
        if (deltaPhiCos >  1.0) deltaPhiCos =  1.0;
        if (deltaPhiCos < -1.0) deltaPhiCos = -1.0;
    
        return std::acos(deltaPhiCos);
    }
    
    
private: 
    Coord _computeTargetRadius(Count nodes, Coord avgDeg, Coord alpha) const;
};

struct Request : public Point {
    CoordInter range;
    
    Request() : Point() {}
    
    Request(const Node id, const Coord phi, const CoordInter range, const SinhCosh r)
        : Point(id, phi, r), range(range)
    {}
    
    Request(const Point& pt, const Geometry& geo,  const SinhCosh band)
        : Point(pt)
    {
        const Coord deltaPhi = geo.deltaPhi(pt.r, band);

        Coord phi = pt.phi;
        //if (phi > 2*M_PI)
        //    phi -= 2*M_PI; // TODO: benchmark against fmod

        range = {phi - deltaPhi, phi + deltaPhi};
    }

    bool operator< (const Request& o) const {
        return std::tie(range.first, id) < std::tie(o.range.first, o.id);
    }

    friend std::ostream& operator <<(std::ostream& stream, const Request& o) {
        return stream << "Request(id:" << o.id << ", phi: " << o.phi << ", range: [" << o.range.first << ", " << o.range.second << "], r: " << std::acosh(o.r.cosh) << ", old:" << o.old() << ")";
    }
};




#endif
