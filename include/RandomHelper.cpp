#include "RandomHelper.hpp"

#include <random>
#include <algorithm>
#include <cassert>
#include "Summation.hpp"

std::vector<uint64_t> RandomHelper::sampleMultinomial(const uint64_t n, const std::vector<double>& probs, DefaultPrng& rg) {
    assert(!probs.empty());
    if (probs.size() == 1) return {n};
    std::vector<uint64_t> result;
    
    assert(std::abs( std::accumulate(probs.cbegin(), probs.cend(), 0.0) - 1.0 ) < std::numeric_limits<double>::epsilon());
    
    auto elems_left = n;
    double prob_mass = 1.0;
    for(size_t i=0; i < probs.size()-1; i++) {
        std::binomial_distribution<uint64_t> distr(elems_left, probs[i] / prob_mass);
        prob_mass -= probs[i];
        
        const auto num = distr(rg);
        elems_left -= num;
        result.push_back(num);
    }
    result.push_back(elems_left);

    assert(n == std::accumulate(result.cbegin(), result.cend(), 0));
    
    return result;
}

std::vector<uint64_t> RandomHelper::sampleMultinomial(const uint64_t n, const uint64_t groups, DefaultPrng& rg) {
    // FIX-ME: This is not really performant
    return sampleMultinomial(n, std::vector<double>(groups, 1.0 / groups), rg);
}
