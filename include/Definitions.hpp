#pragma once
#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <utility>
#include <Vc/Vc>

#include <cstdint>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)       __builtin_expect((x),0)

#define POINCARE

using Coord = double;
using Node = unsigned long int;

using Coord_b = double;
using Coord_v = Vc::Vector<Coord_b>;
using Coord_m = Vc::Mask<Coord_b>;
using Node_b = uint64_t;
using Node_v = Vc::Vector<Node_b>;


constexpr unsigned int CoordPacking = sizeof(Coord_v) / sizeof(Coord_b);
constexpr unsigned int NodePacking = sizeof(Node_v) / sizeof(Node_b);

constexpr unsigned int MinPacking = CoordPacking > NodePacking ? NodePacking : CoordPacking;
constexpr unsigned int MaxPacking = CoordPacking < NodePacking ? NodePacking : CoordPacking;

using Count = Node;
using Seed = uint32_t;

using EdgeId = uint64_t;
using Edge = std::pair<Node, Node>;

template <typename T>
using Interval = std::pair<T, T>;

using CoordInter = Interval<Coord>;

struct SinhCosh {
    Coord invsinh;
    Coord cosh;

    SinhCosh() {}

    SinhCosh(const Coord r)
    {
        const Coord e  = std::exp(r);
        const Coord ie = 1.0 / e;

        invsinh = 2.0 / (e - ie);
        cosh = 0.5 * (e + ie);
    }

    SinhCosh(const Coord invsinh, const Coord cosh)
            : invsinh(invsinh), cosh(cosh)
    {}

    bool operator< (const SinhCosh& o) const {
        return cosh < o.cosh;
    }
};

#endif