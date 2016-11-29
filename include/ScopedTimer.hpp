#pragma once
#ifndef SCOPED_TIMER_HPP
#define SCOPED_TIMER_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>

class ScopedTimer {
    using Clock = std::chrono::high_resolution_clock;

    std::string _prefix;
    Clock::time_point _begin;

    uint64_t _scale;
    double _offset;


public:
    ScopedTimer()
          : _begin(Clock::now()), _scale(1), _offset(0)
    {}

    ScopedTimer(const std::string& prefix, uint64_t scale=0, double offset=0.0)
          : _prefix(prefix), _begin(Clock::now()), _scale(scale), _offset(offset)
    {}

    ~ScopedTimer() {
        if (!_prefix.empty())
            report();
    }

    void start() {
        _begin = Clock::now();
    }


    double report() const {
        return report(_prefix);
    }

    double report(const std::string & prefix) const {
        const auto t2 = Clock::now();
        std::chrono::duration<double> time_span =
                std::chrono::duration_cast<std::chrono::duration<double>>(t2 - _begin);
                
        const double timeUs = (time_span.count()*1e3) - _offset;
        
        if (!_scale) {
            std::cout << prefix << " Time elapsed: " << timeUs << "ms" << std::endl;
        } else {
            std::cout << prefix << " Time elapsed: " << timeUs << "ms / " << _scale << " = " <<  (1e3*timeUs / _scale) << "us" << std::endl;
        }

        return timeUs;
    }
};

#endif
