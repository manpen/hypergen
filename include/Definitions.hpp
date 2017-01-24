#pragma once
#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <utility>
#include <Vc/Vc>

#include <cstdint>

#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)

#define POINCARE
//#define LOG_TRANSFORM

#if !defined(NDEBUG)
#define CROSS_REFERENCE
#endif

#ifdef NDEBUG
#define VERBOSITY(X) 0
#else
#define VERBOSITY(X) X
#endif

using Node = uint32_t; //unsigned long int;
using Coord = double;

#ifdef LOG_TRANSFORM
    using Coord_b = float;
#else
    using Coord_b = Coord;
#endif

using Coord_v = Vc::Vector<Coord_b>;
using Coord_m = Vc::Mask<Coord_b>;
using Node_b = uint32_t;
using Node_v = Vc::Vector<Node_b>;


constexpr unsigned int CoordPacking = sizeof(Coord_v) / sizeof(Coord_b);
constexpr unsigned int NodePacking = sizeof(Node_v) / sizeof(Node_b);

using Count = Node;
using Seed = uint32_t;

using EdgeId = uint64_t;
using Edge = std::pair<Node, Node>;

template <typename T>
using Interval = std::pair<T,T>;

using CoordInter = Interval<Coord_b>;

struct SinhCosh {
    Coord invsinh;
    Coord cosh;

#ifdef CROSS_REFERENCE
    Coord r;
#endif

    SinhCosh() {}

    SinhCosh(const Coord r)
#ifdef CROSS_REFERENCE
        : r(r)
#endif
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