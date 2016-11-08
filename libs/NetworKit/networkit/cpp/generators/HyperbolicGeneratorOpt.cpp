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

	vector<double> bandRadii = getBandRadii(n, R, seriesRatio);
        
        struct Point {
                using T = double;
                
                Point() {}
                
                Point(const T& angle, const T& r, const index& id)
                    : _angle(angle), _rad(r), _id(id)
                {
                    const T e = std::exp(r);
                    const T ie = 1.0 / e;
                    _sinh2 = e - ie;
                    _cosh2 = e + ie;
                }
                
                const index& id() const {return _id;}
                const T& angle() const {return _angle;}
                const T& radius() const {return _rad;}
                const T& sinh2() const {return _sinh2;}
                const T& cosh2() const {return _cosh2;}
                
                bool operator<(const Point& other) const {
                    return std::tie(radius(), angle(), id()) <= std::tie(other.radius(), other.angle(), id());
                }
                
                bool operator<(const T& o) const {
                    return angle() < o;
                }
                
        private:
                T _angle;
                T _rad;
                T _sinh2;
                T _cosh2;
                index _id;
        };
        
        
        std::vector<Point> mapped_points(n);
        {
            Aux::Timer timer;
            timer.start();

            #pragma omp parallel for
            for(index i=0; i < n; i++) {
                const index oid = i; //permutation[i];
                mapped_points[i] = Point(angles[oid], radii[oid], oid);
            }
        
            INFO("Mapping points took ", timer.elapsedMilliseconds(), " milliseconds.");
        }
        
	//Initialize empty bands
	vector<vector<Point>> bands(bandRadii.size() - 1);
	
        //Put points to bands
	#pragma omp parallel for
	for (index j = 0; j < bands.size(); j++) {
		for (const Point& pt : mapped_points){
                        if (pt.radius() >= bandRadii[j] && pt.radius() <= bandRadii[j+1]){
				bands[j].push_back( pt );
			}
		}
	}

	const count bandCount = bands.size();
	const double coshR4 = 4*cosh(R);
	assert(radii.size() == n);

	Aux::Timer bandTimer;
	bandTimer.start();

	//2.Insert edges
	Aux::Timer timer;
	timer.start();
	vector<double> empty;
	GraphBuilder result(n, false, false);

        uint64_t numEdges = 0;
        uint64_t numCompares = 0;

        auto tryAddEdge = [&] (const Point& pt, const Point& candidate) { 
                if (candidate < pt)
                    return false;
                
                double deltaPhi = M_PI - abs(M_PI - abs(pt.angle() - candidate.angle()));
                
                if (pt.cosh2() * candidate.cosh2() - pt.sinh2()*candidate.sinh2()*cos(deltaPhi) <= coshR4) {
                    return true;
                }
                
                return false;
        }; 
        
        
        const double coshR = std::cosh(R);
        std::vector<double> coshBandRadii(bandRadii);
        for(auto& x: coshBandRadii) x = 0.5 * std::cosh(x);
        std::vector<double> sinhBandRadii(bandRadii);
        for(auto& x: sinhBandRadii) x = 0.5 * std::sinh(x);
        
        auto getMinMaxTheta = [&] (const Point pt, unsigned int bandIdx ) -> std::pair<double, double> {
	  if (!bandIdx)
              return {0.0, 2.0*M_PI};

          double a;
          if (pt.radius() > bandRadii[bandIdx])
            a = (pt.cosh2() * pt.cosh2() * 0.25 - coshR) / (pt.sinh2() * pt.sinh2() * 0.25);
          else
            a = (pt.cosh2() * coshBandRadii[bandIdx] - coshR) / (pt.sinh2() * sinhBandRadii[bandIdx]);
          
	  if(a < -1.0) a = -1.0;
	  else if(a > 1.0) a = 1.0;
          
          a = acos(a);
          
	  return {pt.angle() - a, pt.angle() + a};
	};
        
        
        // FIX-ME: This shadows the threadTimes in order to ma
        //std::vector<Aux::Timer> threadtimers(omp_get_max_threads());
        #pragma omp parallel reduction(+:numEdges,numCompares)
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
                
		#pragma omp for schedule(guided) nowait
		for (index i = 0; i < n; i++) {
                        const auto& pt = mapped_points[i];
                        
			for(index j = bandCount; --j;) {
				if(bandRadii[j+1] > pt.radius()){
                                        const auto& bandPoints = bands[j];
                                    
					double minTheta, maxTheta;
					std::tie (minTheta, maxTheta) = getMinMaxTheta(pt, j);

                                        
                                        using cit = typename std::vector<Point>::const_iterator;
                                        auto check_range = [&] (const cit& begin, const cit& end) {
                                            for(auto it = begin; it != end; it++) {
                                                numEdges += tryAddEdge(pt, *it);
                                                ++numCompares;
                                            }
                                        };
                                        
                                        
                                        //Case 1: We do not have overlap 2pi, simply put all the points between min and max to the list
                                        if(maxTheta <= 2.0*M_PI && minTheta >= 0) {
                                                const auto low = std::lower_bound(bandPoints.begin(), bandPoints.end(), minTheta, 
                                                    [] (const Point& pt, const double& t) {return pt.angle() < t;}
                                                );
                                                
                                                const auto high = std::upper_bound(low, bandPoints.end(), maxTheta, 
                                                    [] (const double& t, const Point& pt) {return t < pt.angle();}
                                                );
                                                
                                                // TODO: since we're touching each point anyways a linear search for the upper_bound in
                                                // check_range may be faster;
                                                check_range(low, high);
                                            
                                        } 
                                        //Case 2: We have an overlap with 0 or 2pi
                                        else {
                                                //1. Get points from minTheta to 2pi
                                                const auto minTh = fmod(minTheta + 2.0*M_PI, 2.0*M_PI);
                                                const auto maxTh = fmod(maxTheta           , 2.0*M_PI);
                                                
                                                const auto low = std::lower_bound(bandPoints.begin(), bandPoints.end(), minTh, 
                                                    [] (const Point& pt, const double& t) {return pt.angle() < t;}
                                                );
                                                
                                                const auto high = std::upper_bound(bandPoints.begin(), low, maxTh, 
                                                    [] (const double& t, const Point& pt) {return t < pt.angle();}
                                                );


                                                check_range(low, bandPoints.end());
                                                check_range(bandPoints.begin(), high);
                                        }
                                } else {
                                    break;
                                }
			}
		}

		threadtimers[id].stop();
	}
	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
        INFO("Needed ", numCompares, " compared (",  (100.0 * numEdges / numCompares), " successful");
        INFO("Produced ", numEdges, " edges");
        
	return result.toGraph(!directSwap, true);
}

}
