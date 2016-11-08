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
