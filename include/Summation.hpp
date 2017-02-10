/**
 * @file
 * @brief Summation
 *
 * A simple stream sum with push interface, supporting numerically stable
 * Kahn sums.
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
#ifndef SUMMATION_HPP
#define SUMMATION_HPP

template <typename T>
class Summation {
public:
    using value_type = T;

    Summation(const T& init = 0.0)
            : _sum(init)
    {}


    const value_type& sum() const {
        return _sum;
    }

    void push(const value_type& x) {
        _sum += x;
    }

protected:
    value_type _sum{0.0};

};


template <typename T>
class KahnSummation {
public:
    using value_type = T;

    KahnSummation(const T& init = 0.0)
            : _sum(init), _carry(0.0)
    {}

    const value_type& sum() const {
        return _sum;
    }

    void push(const value_type& in) {
        const auto y = in - _carry;
        const auto tmp = _sum + y;
        _carry = (tmp - _sum) - y;
        _sum = tmp;
    }
    
    template<typename Iterator>
    void push(const Iterator begin, const Iterator end) {
        for(auto it = begin; it != end; ++it)
            push(*it);
    }

protected:
    value_type _sum;
    value_type _carry;

};


#endif // SUMMATION_HPP
