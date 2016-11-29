#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <algorithm>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

Graph HyperbolicGenerator::generateColdOpt(const vector<double> &angles, const vector<double> &radii, double R, double seriesRatio) const {
    INFO("Start Opt.");
    const count n = angles.size();
    assert(radii.size() == n);

    for (index i = 0; i < n; i++) {
        assert(radii[i] < R);
    }

    //Map points
    std::vector<double> pt_cosh(n);
    std::vector<double> pt_invsinh(n);

    #pragma omp parallel for simd
    for (index i = 0; i < n; i++) {
        const double e = std::exp(radii[i]);
        const double ie = 1.0 / e;
        pt_invsinh[i] = (2.0 / (e - ie));
        pt_cosh[i] = (0.5 * (e + ie));
    }

    //Put points to bands
    struct Band {
        std::vector<index>  pt_id;
        std::vector<double> pt_theta;
        std::vector<double> pt_cosh;
        std::vector<double> pt_invsinh;
        
        unsigned int band_id;
        double rad_cosh;
        double rad_invsinh;
        
        double avg_spacing;
    };
    
    auto bandRadii = getBandRadii(n, R, seriesRatio);
    std::vector<Band> bands(bandRadii.size());
    {
        Aux::Timer timer;
        timer.start();

        #pragma omp parallel for
        for (index j = 0; j < bands.size(); j++) {
            Band& band = bands[j];
            
            if (j+1 < bands.size()) {
                for (index i = 0; i < n; i++) {
                    if (radii[i] >= bandRadii[j] && radii[i] < bandRadii[j+1]) {
                        band.pt_id.push_back(i);
                        band.pt_theta.push_back(angles[i]);
                        band.pt_invsinh.push_back(pt_invsinh[i]);
                        band.pt_cosh.push_back(pt_cosh[i]);
                    }
                }
            }
            
            band.band_id = j;
            band.rad_cosh = std::cosh(bandRadii[j]);
            band.rad_invsinh = 1.0 / std::sinh(bandRadii[j]);
            band.avg_spacing = band.pt_id.size() / (2.0 * M_PI);
        }

        timer.stop();
        INFO("Distributing points took ", timer.elapsedMilliseconds(), "ms");
    }
    
    const count bandCount = bands.size()-1;
    const double coshR = cosh(R);
    assert(radii.size() == n);

    uint64_t numEdges = 0;
    uint64_t numCompares = 0;

    auto getMinMaxTheta = [&] (const double pt_cosh, const double pt_invsinh, const double pt_theta, const Band & band) {
        if (!band.band_id)
            return std::make_pair(0.0, 2.0*M_PI);

        double a;
        if (pt_cosh > band.rad_cosh)
            a = (pt_cosh * pt_cosh - coshR) * pt_invsinh * pt_invsinh;
        else
            a = (pt_cosh * band.rad_cosh - coshR) * pt_invsinh * band.rad_invsinh;

        if     (a < -1.0) a = -1.0;
        else if(a >  1.0) a = 1.0;

        a = acos(a);

        return std::make_pair(pt_theta - a, pt_theta + a);
    };

    Aux::Timer timer;
    timer.start();    
    #pragma omp parallel reduction(+:numEdges,numCompares)
    {
        index id = omp_get_thread_num();
        threadtimers[id].start();

        #pragma omp for schedule(guided) nowait
        for (index i = 0; i < n; i++) {
            for(index j = bandCount; j--;) {
                if(bands[j+1].rad_cosh < pt_cosh[i])
                    break;

                const Band& band = bands[j];

                double minTheta, maxTheta;
                std::tie (minTheta, maxTheta) = getMinMaxTheta(pt_cosh[i], pt_invsinh[i], angles[i], band);

                using cit = typename std::vector<double>::const_iterator;
                auto check_range = [&] (const cit& begin, const cit& end) {
                    const size_t ibegin = std::distance(band.pt_theta.begin(), begin);
                    const size_t iend   = std::distance(band.pt_theta.begin(), end);

                    // for each candidate
                    for(auto k = ibegin; k != iend; k++) {
                        ++numCompares;

                        if (band.pt_id[k] == i || std::tie(band.pt_cosh[k], band.pt_theta[k]) < std::tie(pt_cosh[i], angles[i]))
                            continue;
                        
                        

                        const double deltaPhi = M_PI - abs(M_PI - abs(angles[i] - band.pt_theta[k]));
                        if ( (pt_cosh[i]*band.pt_cosh[k] - coshR) * pt_invsinh[i]*band.pt_invsinh[k] <= cos(deltaPhi) ) {
                            ++numEdges;
                        }
                    }
                };
                
                //Case 1: We do not have overlap 2pi, simply put all the points between min and max to the list
                if(maxTheta <= 2.0*M_PI && minTheta >= 0) {
                    const auto low = std::lower_bound(band.pt_theta.cbegin(), band.pt_theta.cend(), minTheta);
                    const auto high = std::upper_bound(low, band.pt_theta.cend(), maxTheta);

                    // TODO: since we're touching each point anyways a linear search for the upper_bound in
                    // check_range may be faster;
                    check_range(low, high);
                }
                //Case 2: We have an overlap with 0 or 2pi
                else {
                    //1. Get points from minTheta to 2pi
                    minTheta = fmod(minTheta + 2.0*M_PI, 2.0*M_PI);
                    maxTheta = fmod(maxTheta           , 2.0*M_PI);

                    const auto low = std::lower_bound(band.pt_theta.cbegin(), band.pt_theta.cend(), minTheta);
                    const auto high = std::upper_bound(band.pt_theta.cbegin(), band.pt_theta.cend(), maxTheta);

                    check_range(low, band.pt_theta.cend());
                    check_range(band.pt_theta.cbegin(), high);
                }
            }
        }

        threadtimers[id].stop();
    }
    timer.stop();
    INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
    INFO("Needed ", numCompares, " compared (",  (100.0 * numEdges / numCompares), "% successful)");
    INFO("Produced ", numEdges, " edges");
    
    return {};
}

}
