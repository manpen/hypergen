/**
 * @file
 * @brief Histogram
 *
 * A integer histogram that can be disabled via a template parameter
 * with zero-overhead.
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

template<>
class Histogram<false> {
public:
    void addPoint(Count key) {}

    friend std::ostream& operator <<(std::ostream& stream, const Histogram<false>& o) {
        return stream;
    }

    void toStream(std::ostream& stream, std::string s) const {}

    Histogram operator+(const Histogram& o) const {}
};


#endif