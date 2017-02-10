/**
 * @file
 * @brief Implementation of RandomHelper
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
#include "RandomHelper.hpp"

#include <random>
#include <algorithm>
#include "Assert.hpp"

#include "Summation.hpp"

std::vector<Node> RandomHelper::sampleMultinomial(const Node n, const std::vector<Coord>& probs, DefaultPrng& rg) {
    ASSERT(!probs.empty());
    if (probs.size() == 1) return {n};
    std::vector<Node> result;
    
    ASSERT_LS(std::abs( std::accumulate(probs.cbegin(), probs.cend(), 0.0) - 1.0 ), 16 * std::numeric_limits<Coord>::epsilon());
    
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

    ASSERT_EQ(n, std::accumulate(result.cbegin(), result.cend(), 0));
    
    return result;
}

std::vector<Node> RandomHelper::sampleMultinomial(const Node n, const Node groups, DefaultPrng& rg) {
    // FIX-ME: This is not really performant
    return sampleMultinomial(n, std::vector<Coord>(groups, 1.0 / groups), rg);
}
