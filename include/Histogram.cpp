#include "Histogram.hpp"
#include <algorithm>

void Histogram::_dump(std::ostream &stream, std::string label) const {
    const Count totalCount = std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [] (const Count c, const auto & it) {return c + it.second;});

    const Count mean = static_cast<double>(std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [] (const Count c, const auto & it) {return c + it.first * it.second;})) / totalCount;

    const double var = sqrt(static_cast<double>(std::accumulate(_map.cbegin(), _map.cend(), Count(0),
        [mean] (const Count c, const auto & it) {return c + it.second * (it.first - mean)*(it.first - mean);})) / totalCount);

    for(const auto& v : _map) {
        Coord key, freq;
        std::tie(key, freq) = v;

        stream << key << " "
               << freq << " "
               << (static_cast<double>(freq) / totalCount)
               << " # " << label << "-HIST \n";
    }

    stream
            << "Total count: " << totalCount << " # " << label << "\n"
            << "Mean:        " << mean << " # " << label << "\n"
            << "Variance:    " << var << " # " << label
    << std::endl;
}

Histogram Histogram::operator+(const Histogram &o) const {
    Histogram result(*this);

    for(const auto& v : o._map) {
        result._map[v.first] += v.second;
    }

    return result;
}