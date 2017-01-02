#pragma once
#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "Geometry.hpp"
#include <ostream>
#include <map>

template<bool Enabled=true>
class Histogram {
public:
    constexpr static bool enabled {Enabled};

    void addPoint(Count key) {
        if (Enabled)
            _map[key]++;
    }

    friend std::ostream& operator <<(std::ostream& stream, const Histogram& o) {
        if (Enabled)
            o._dump(stream, "");
        return stream;
    }

    void toStream(std::ostream& stream, std::string s) const {
        if (Enabled)
            _dump(stream, s);
    }

    Histogram operator+(const Histogram& o) const;

protected:
    std::map<Count, Count> _map;

    void _dump(std::ostream& stream, std::string s) const;
};


#endif