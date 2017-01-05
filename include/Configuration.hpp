#pragma once
#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "Definitions.hpp"

struct Configuration {
    Node nodes {1000};
    Coord avgDegree {10.0};
    Coord degreeExp {3.0};
    Coord alpha {1.0};
    Seed seed {1234};

    unsigned int noWorker;
    unsigned int noSegments;

    unsigned int verbosity {0};

    unsigned int activeUpdateInterval {1};

    Coord R {-1};

    enum BandLimitType {BandExp, BandLin};
    BandLimitType bandLimits {BandLin};

    Coord bandExpFactor {1.0};
    Coord bandLinFactor {3.0};

    Configuration() = delete;
    Configuration(int argc, char* argv[]);

    void dump() const;
};


#endif // CONFIGURATION_HPP