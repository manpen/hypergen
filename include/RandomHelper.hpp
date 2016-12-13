#pragma once
#ifndef RANDOM_HELPER_HPP
#define RANDOM_HELPER_HPP

#include <vector>
#include <cstdint>
#include <random>

#include "Geometry.hpp"

using DefaultPrng = std::mt19937_64;

class RandomHelper {
public:    
    static std::vector<Node> sampleMultinomial(const Node n, const std::vector<Coord>& probs, DefaultPrng& rg);
    static std::vector<Node> sampleMultinomial(const Node n, const Node groups, DefaultPrng& rg);
};


#endif
