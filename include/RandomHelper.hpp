/**
 * @file
 * @brief RandomHelper
 *
 * Tools for some random distributions.
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
