/**
 * @file
 * @brief Point
 *
 * A point in hyperbolic geometry. Depending on the definitions
 * in Definitions.hpp, paramters such as Poincare mappings are
 * computed.
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
#pragma once
#ifndef POINT_HPP
#define POINT_HPP

#include "Definitions.hpp"
#include "Geometry.hpp"
#include "Assert.hpp"
#include <iostream>

struct Point {
    Node id;
    float phi;
    SinhCosh r;

#ifdef POINCARE
    Coord_b poinX;
    Coord_b poinY;

#ifdef LOG_TRANSFORM
    Coord_b poinLogInvLen;
#else
    Coord poinInvLen;
#endif
#endif

    static constexpr Node OLD_MASK  = Node(1) << (sizeof(Node)*8-1);
    static constexpr Node NODE_MASK = ~OLD_MASK;

    Point() {}

    Point(const Node id, const Coord phi, const SinhCosh r)
            : id(id), phi(phi), r(r)
    {
        /*
         * poinLen = sqrt[(cosh-1)/(cosh+1)] = sqrt[1 - 2/(cosh+1)]
         *         = sqrt[1 - 1.0 / X] with X = (cosh + 1)/2
         * dist = (...) / (1 - poinLen**2)
         *      = (...) / (1 - (1-1/X))
         *      = (...) / (1 - (1-2/(cosh+1))
         *      = (...) / (2 / (cosh + 1)
         *      = (...) * 0.5 * (cosh + 1)
         *      = (...) * X
         */
#ifndef SKIP_DIST_COMP
#ifdef POINCARE
        const Coord pInvLen = (r.cosh + 1.0) * 0.5;
        const Coord poinR = (std::sqrt( (1.0 - 1.0 / pInvLen)));

#ifdef LOG_TRANSFORM
        poinLogInvLen = std::log(pInvLen);
#else
        poinInvLen = pInvLen;
#endif

        ASSERT_LS(poinR, 1.0);

        poinX = (poinR * std::sin(phi));
        poinY = (poinR * std::cos(phi));

#endif // POINCARE
#endif //SKIP_DIST_COMP
    }

    void setOld() {
        id |= OLD_MASK;
    }

    bool old() const {
        return id & OLD_MASK;
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
#ifdef LOG_TRANSFORM
        return 1.0 + std::exp(2.0*(diffX*diffX + diffY*diffY) - poinLogInvLen - o.poinLogInvLen);
#else
        return 1.0 + 2.0*(diffX*diffX + diffY*diffY) * poinInvLen * o.poinInvLen;
#endif
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

    Request(const Point& pt, const CoordInter& range)
            : Point(pt), range(range)
    {}

    bool operator< (const Request& o) const {
        return std::tie(range.first, id) < std::tie(o.range.first, o.id);
    }

    bool operator> (const Request& o) const {
        return std::tie(range.first, id) > std::tie(o.range.first, o.id);
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
