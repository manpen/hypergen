/**
 * @file
 * @brief Configuration (Command-Line Parser)
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

    Coord bandLinFactor {3.0};

    Configuration() = delete;
    Configuration(int argc, char* argv[]);

    void dump() const;
};


#endif // CONFIGURATION_HPP