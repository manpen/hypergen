/**
 * @file
 * @brief Implementation of Histogram
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

#include "Histogram.hpp"

#include <algorithm>
#include <numeric>

template<bool Enabled>
void Histogram<Enabled>::_dump(std::ostream &stream, std::string label) const {
    if (!Enabled) return;

    const Count totalCount = std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [] (const Count c, const auto & it) {return c + it.second;});

    const Count totalWeight = std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [] (const Count c, const auto & it) {return c + it.first * it.second;});

    const Count mean = static_cast<double>(totalWeight) / totalCount;

    const double var = sqrt(static_cast<double>(std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [mean] (const Count c, const auto & it) {return c + it.second * (it.first - mean)*(it.first - mean);})) / totalCount);

    for(const auto& v : _map) {
        Coord key, freq;
        std::tie(key, freq) = v;

        stream << key << " "
               << freq << " "
               //<< (static_cast<double>(freq) / totalCount)
               << "# " << label << "\n";
    }

    stream
            << "Total count: " << totalCount << " # SUM-" << label << "\n"
            << "Total weight:" << totalWeight << " # SUM-" << label << "\n"
            << "Mean:        " << mean << " # SUM-" << label << "\n"
            << "Variance:    " << var << " # SUM-" << label
    << std::endl;
}

template<bool Enabled>
Histogram<Enabled> Histogram<Enabled>::operator+(const Histogram<Enabled> &o) const {
    if (!Enabled) return {};

    Histogram result(*this);

    for(const auto& v : o._map) {
        result._map[v.first] += v.second;
    }

    return result;
}

// we have two options for the template parameter (no shit) ...
// so let's compile for both variants
template class Histogram<false>;
template class Histogram<true>;