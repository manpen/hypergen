#include "RandomHelper.hpp"

#include <random>
#include <algorithm>
#include <cassert>
#include "Summation.hpp"

std::vector<Node> RandomHelper::sampleMultinomial(const Node n, const std::vector<Coord>& probs, DefaultPrng& rg) {
    assert(!probs.empty());
    if (probs.size() == 1) return {n};
    std::vector<Node> result;
    
    assert(std::abs( std::accumulate(probs.cbegin(), probs.cend(), 0.0) - 1.0 ) < std::numeric_limits<Coord>::epsilon());
    
    auto elems_left = n;
    double prob_mass = 1.0;
    for(size_t i=0; i < probs.size()-1; i++) {
        std::binomial_distribution<Node> distr(elems_left, probs[i] / prob_mass);
        prob_mass -= probs[i];
        
        const auto num = distr(rg);
        elems_left -= num;
        result.push_back(num);
    }
    result.push_back(elems_left);

    assert(n == std::accumulate(result.cbegin(), result.cend(), 0));
    
    return result;
}

std::vector<Node> RandomHelper::sampleMultinomial(const Node n, const Node groups, DefaultPrng& rg) {
    // FIX-ME: This is not really performant
    return sampleMultinomial(n, std::vector<Coord>(groups, 1.0 / groups), rg);
}
