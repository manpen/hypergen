#pragma once
#ifndef RANDOM_HELPER_HPP
#define RANDOM_HELPER_HPP

#include <vector>
#include <cstdint>
#include <random>

using DefaultPrng = std::mt19937_64;

class RandomHelper {
public:    
    static std::vector<uint64_t> sampleMultinomial(const uint64_t n, const std::vector<double>& probs, DefaultPrng& rg);
    static std::vector<uint64_t> sampleMultinomial(const uint64_t n, const uint64_t groups, DefaultPrng& rg);
};


#endif
