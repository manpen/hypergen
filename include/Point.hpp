#pragma once
#ifndef POINT_HPP
#define POINT_HPP

#include "Definitions.hpp"
#include "Geometry.hpp"
#include <iostream>

struct Point {
    Node id;
    Coord phi;
    SinhCosh r;

#ifdef POINCARE
    Coord poinR;
    Coord poinX;
    Coord poinY;
    Coord poinInvLen;
#endif

    bool _old{false};


    Point() {}

    Point(const Node id, const Coord phi, const SinhCosh r)
            : id(id), phi(phi), r(r)
#ifdef POINCARE
            , poinR(std::sqrt( (r.cosh - 1.0) / (r.cosh + 1.0) ))
            , poinX(poinR * std::sin(phi))
            , poinY(poinR * std::cos(phi))
            , poinInvLen(1.0 / (1.0 - poinR * poinR))
#endif
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

    Coord distanceToHyper(const Point& o) const {
        return std::acosh(r.cosh * o.r.cosh - cos(o.phi - phi) / r.invsinh / o.r.invsinh);
    }

#ifdef POINCARE
    Coord distanceToPoincare(const Point& o) const {
        return std::acosh(1.0 + 2.0*(poinX*o.poinX + poinY*o.poinY) * poinInvLen * o.poinInvLen);
    }

    Coord distanceTo(const Point& o) const {
        return distanceToPoincare(o);
    }
#else
    Coord distanceTo(const Point& o) const {
        return distanceToHyper(o);
    }
#endif
};

struct Request : public Point {
    CoordInter range;

    Request() : Point() {}

    Request(const Node id, const Coord phi, const CoordInter range, const SinhCosh r)
            : Point(id, phi, r), range(range)
    {}

    Request(const Point& pt, const Geometry& geo, const SinhCosh band)
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

    static Request faraway() {
        return {0, 0.0, {0.0, 0.0}, 1e10};
    }
};

#endif
