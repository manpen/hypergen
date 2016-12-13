#pragma once
#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "Geometry.hpp"
#include <ostream>
#include <map>


class Histogram {
public:

    Histogram() {}

    void addPoint(Count key) {
        _map[key]++;
    }

    friend std::ostream& operator <<(std::ostream& stream, const Histogram& o) {
        o._dump(stream, "");
        return stream;
    }

    void toStream(std::ostream& stream, std::string s) const {
        _dump(stream, s);
    }

    Histogram operator+(const Histogram& o) const;


protected:
    std::map<Count, Count> _map;

    void _dump(std::ostream& stream, std::string s) const;
};


#endif