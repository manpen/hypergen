#include "RandomHelper.hpp"

#include <random>
#include <algorithm>
#include <cassert>
#include "Summation.hpp"

std::vector<uint64_t> RandomHelper::sampleMultinomial(const uint64_t n, const std::vector<double>& probs, DefaultPrng& rg) {
    assert(probs.size() > 1);
    std::vector<uint64_t> result(probs.size());
    
    double normalization;
    {
        KahnSummation<double> sum; 
        for(const auto p : probs) {
            assert(p >= 0.0);
            sum.push(p);
        }
        normalization = 1.0 / sum.sum();
    }
    
    auto elems_left = n;
    for(auto it = probs.cbegin(); std::distance(probs.cend(), it) > 1; ++it) {
        std::binomial_distribution<uint64_t> distr(elems_left, (*it * normalization));
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
