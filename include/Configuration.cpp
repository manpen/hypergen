/**
 * @file
 * @brief Implementation of Configuration (Command-Line Parser)
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

#include "Configuration.hpp"
#include "Geometry.hpp"
#include <iostream>
#include <omp.h>

Configuration::Configuration(int argc, char* argv[]) {
    noWorker = noSegments = omp_get_max_threads();

    for(unsigned int i=1; i+1 < argc; i+=2) {
        std::string key = argv[i];
        std::string value = argv[i + 1];
        if (key == "-n") { nodes = stoll(value); }
        else if (key == "-d") { avgDegree = stod(value); }
        else if (key == "-e") { degreeExp = stod(value); }
        else if (key == "-s") { seed = stoi(value); }
        else if (key == "-p") { noSegments = stoi(value); }
        else if (key == "-w") { noWorker = stoi(value); }
        else if (key == "-v") { verbosity = stoi(value); }
        else if (key == "-R") { R = stod(value); }
        else if (key == "-i") { activeUpdateInterval = stoi(value); }

        else if (key == "-m") {
            if (value == "l" || value == "L") {
                bandLimits = BandLimitType::BandLin;
            } else if (value == "e" || value == "E") {
                bandLimits = BandLimitType::BandExp;
            } else {
                std::cerr << "Only support Bandlimit types l/L or e/E" << std::endl;
                abort();
            }


        }
        else if (key == "-b") {
            bandExpFactor = bandLinFactor = stod(value);
        } else {
            std::cerr << "Unknown argument: " << key << std::endl;
            abort();
        }
    }

    alpha = 0.5 * (degreeExp - 1.0);

    if (R > 0) {
        avgDegree = Geometry::getExpectedDegree(R, alpha, nodes);
        std::cerr << "WARNING: Since R is provided, we update the avgDegree parameter to " << avgDegree << std::endl;
    }

    dump();
}

void Configuration::dump() const {
    std::cout <<
        "conf.nodes:               " << nodes                << "\n"
        "conf.avgDegree:           " << avgDegree            << "\n"
        "conf.degreeExp:           " << degreeExp            << "\n"
        "conf.alpha:               " << alpha                << "\n"
        "conf.R:                   " << R                    << "\n\n"

        "conf.seed:                " << seed                 << "\n"
        "conf.actveUpdateInterval: " << activeUpdateInterval << "\n"
        "conf.noWorker:            " << noWorker             << "\n"
        "conf.noSegments:          " << noSegments           << "\n"
        "conf.verbosity:           " << verbosity            << "\n"

        "conf.bandLimits:          " << (bandLimits == BandLimitType::BandLin ? "LIN" : "EXP") << "\n"
        "conf.bandExpFactor:       " << bandExpFactor        << "\n"
        "conf.bandLinFactor:       " << bandLinFactor        << "\n"
    << std::endl;
}