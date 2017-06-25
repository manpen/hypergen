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

        else if (key == "-b") {
            bandLinFactor = stod(value);
        }
        else if (key == "-h") {
            dump();
            abort();
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
    "Configuration\n"
        " conf.nodes:                        -n " << nodes                << "\n"
        " conf.avgDegree:                    -d " << avgDegree            << "\n"
        " conf.degreeExp:                    -e " << degreeExp            << "\n"
        " conf.alpha: (derived from -e)         " << alpha                << "\n"
        " conf.R: (default: derived from -d) -R " << R                    << "\n\n"

        " conf.seed:                         -s " << seed                 << "\n"
        " conf.activeUpdateInterval:         -i " << activeUpdateInterval << "\n"
        " conf.noWorker:                     -w " << noWorker             << "\n"
        " conf.noSegments:                   -p " << noSegments           << "\n"
        " conf.verbosity:                    -v " << verbosity            << "\n"

        " conf.bandLinFactor:                -b " << bandLinFactor        << "\n"
    << std::endl;
}