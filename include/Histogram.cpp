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