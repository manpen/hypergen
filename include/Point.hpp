#pragma once
#ifndef POINT_HPP
#define POINT_HPP

#include "Definitions.hpp"
#include "Geometry.hpp"
#include "Assert.hpp"
#include <iostream>

struct Point {
    Node id;
    Coord phi;
    SinhCosh r;

#ifdef POINCARE
//    Coord poinR;
    Coord poinX;
    Coord poinY;
    Coord poinInvLen;
#endif

    bool _old{false};


    Point() {}

    Point(const Node id, const Coord phi, const SinhCosh r)
            : id(id), phi(phi), r(r)
    {
#ifdef POINCARE
        const Coord poinRR = (r.cosh - 1.0) / (r.cosh + 1.0);
        const Coord poinR = (std::sqrt( poinRR ));

        ASSERT_LS(poinR, 1.0);
        ASSERT_LS(poinRR, 1.0);

        poinX = (poinR * std::sin(phi));
        poinY = (poinR * std::cos(phi));

        poinInvLen = (1.0 / (1.0 - poinRR));
#endif
    }

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
        return stream << "Point(id:" << o.id << ", "
                "phi: " << o.phi << ", "
                "r: " << std::acosh(o.r.cosh) << ", "
#ifdef CROSS_REFERENCE
                "r.r: " << o.r.r << ", "
#endif
                "old:" << o.old() << ")";
    }

    Coord distanceToHyper(const Point& o) const {
        return std::acosh(coshDistanceToHyper(o));
    }

    Coord coshDistanceToHyper(const Point& o) const {
        return r.cosh * o.r.cosh - std::cos(o.phi - phi) / r.invsinh / o.r.invsinh;
    }


#ifdef POINCARE
    Coord distanceToPoincare(const Point& o) const {
        return std::acosh(coshDistanceToPoincare(o));
    }

    Coord coshDistanceToPoincare(const Point& o) const {
        const auto diffX = poinX - o.poinX;
        const auto diffY = poinY - o.poinY;
        return 1.0 + 2.0*(diffX*diffX + diffY*diffY) * poinInvLen * o.poinInvLen;
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
        return stream << "Request(id:" << o.id << ", "
                         "phi: " << o.phi << ", "
                         "range: [" << o.range.first << ", " << o.range.second << "], "
                         "r: " << std::acosh(o.r.cosh) << ", "
#ifdef CROSS_REFERENCE
                         "r.r: " << o.r.r << ", "
#endif
                         "old:" << o.old() << ")";
    }

    static Request faraway() {
        return {0, 0.0, {0.0, 0.0}, std::numeric_limits<Coord>::max()};
    }
};

#endif
