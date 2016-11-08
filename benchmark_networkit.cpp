#include <iostream>
#include <NetworKit/generators/HyperbolicGenerator.h>
#include <NetworKit/auxiliary/Timer.h>
#include <NetworKit/graph/GraphBuilder.h>
#include <NetworKit/auxiliary/Parallel.h>

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <algorithm>
#include <utility>




namespace NetworKit {
Graph HyperbolicGenerator::generateColdOpt1(const vector<double>& angles, const vector<double>& radii, double R, double seriesRatio) const {
        INFO("Start Opt.");    
	const count n = pt_angles.size();
	const double coshR = std::cosh(R);


	for (index i = 0; i < n; i++) {
            assert(pt_radii[i] < R);
	}
	
        struct Request {
            double start;
            double stop;
            index id;
            double theta;
            double cosh;
            double sinh;
            
            Request(const double& start, const double& stop, const index& id, const double& theta, const double& cosh, const double& sinh)
                : start(start), stop(stop), id(id), theta(theta), cosh(cosh), sinh(sinh)
            {}
        };
        
	//Initialize empty bands
	struct Band {
            const double lower;
            const double lowerCosh;
            const double lowerSinh;
            
            Band(double l) 
                : lower(l), lowerCosh(0.5 * std::cosh(l)), lowerSinh(0.5 * std::sinh(l))
            {}
            
            // points in band
            std::vector<index>  pt_id;
            std::vector<double> pt_theta;
            std::vector<double> pt_sinh;
            std::vector<double> pt_cosh;

            void insertPoint(const index& id, const double& theta, const double& cosh, const double& sinh) {
                pt_id.push_back(id);
                pt_theta.push_back(theta);
                pt_cosh.push_back(cosh);
                pt_sinh.push_back(sinh);
            }
            
            std::vector<Request> requests;

            void complete() {
                std::sort(requets.begin(), requets.end(), [] (const auto& a, const auto&b) {return std::get<0>(a) < std::get<0>(b);});
            }

            void issueRequest(const index& id, const double& theta, const double& cosh, const double& sinh, const double& coshR) {
            }
            
            std::pair<double, double> getMinMaxTheta(const double& pt_angle, const double& pt_cosh, const double& pt_sinh, const double& coshR) const {
                if (!bandIdx)
                    return {0.0, 2.0*M_PI};

                double a;
                if (pt_cosh > lowerCosh)
                    a = (pt_cosh * pt_cosh - coshR) / (pt_sinh * pt_sinh);
                else
                    a = (pt_cosh * lowerCosh - coshR) / (pt_sinh * lowerSinh);
                
                if(a < -1.0) a = -1.0;
                else if(a > 1.0) a = 1.0;
                
                a = acos(a);
                
                return {pt_angle - a, pt_angle + a};
            };            
        };
        
        struct Segment {
            const double thetaStart;
            const double thetaStop;
            const double segmentWidth;
            const unsigned int worker_id;
            const unsigned int worker;
            
            const index numPoints;
            
            const unsigned int numBands;
            const std::vector<double>& bandRadii;
            std::vector<Band> bands;

            std::vector<std::vector<std::vector<Request>>> insBuffers;
            
            Segment(const double& start,
                    const double& stop,
                    const unsigned int& worker_id,
                    const unsigned int& worker,
                    std::vector<Segments>& segments,
                    const index numPoints,
                    const std::vector<double>& bandRadii)
            
                : thetaStart(start)
                , thetaStop(stop)
                , segmentWidth(stop - start)
                , worker_id(worker_id)
                , worker(worker)
                , numPoints(numPoints)
                , numBands(bandRadii.size() - 1)
                , bandRadii(bandRadii)
                , insBuffer(worker)
            {
                bands.reserve(numBands);
                for(unsigned int b=0; b<bands.size(); b++) {
                    bands.emplace_back( Band(bandRadii[b]) );
                }
                
                for(auto& w : insBuffers)
                    w.resize(worker);
            }
                
            void loadPoints(index id, const double* radii, const double* angles) {
                const index end = id + numPoints;
                
                for(; id != end; ++id, ++radii, ++angles) {
                    const unsigned int bandIdx = 
                        std::distance(bandRadii.cbegin(), std::lower_bound(bandRadii.cbegin(), bandRadii.cend(), *radii));

                    const double& theta = *angles;
                        
                    const double e = std::exp(*radii);
                    const double ie = 1.0 / e;
                    const double cosh = 0.5*(e - ie);
                    const double sinh = 0.5*(e + ie);
                    
                    bands[bandIdx].insertPoint(id, theta, cosh, sinh);
                    
                    // issue requests
                    for(unsigned int b=bandIdx; b < numBands; ++b) {
                        double minTheta, maxTheta;
                        std::tie(minTheta, maxTheta) = bands[b].getMinMaxTheta(theta, cosh, sinh, coshR);

                        if(thetaStart <= minTheta && thetaStop > maxTheta) {
                            // default, only care about our own segment
                            bands[b].requests.emplace_back(minTheta, maxTheta, id, theta, cosh, sinh);
                        } else {
                            auto insert = [&] (const double& start, const double& stop) {
                                const unsigned int s = start / segmentWidth;
                                double ss = s * segmentWidth;
                                double se = ss + segmentWidth;
                                
                                for(unsigned int seg = s; seg < worker; seg++, ss = se, se += segmentWidth) {
                                    if (ss > stop)
                                        break;
                                    
                                    insBuffers[seg][b].emplace_back(std::max(start, ss), std::min(stop, ss), id, theta, cosh, sinh);
                                }
                            };
                            
                        
                            if(maxTheta <= 2.0*M_PI && minTheta >= 0) {
                                insert(minTheta, maxTheta);
                            } else {
                                minTheta = fmod(minTheta + 2.0*M_PI, 2.0*M_PI);
                                maxTheta = fmod(maxTheta           , 2.0*M_PI);

                                insert(minTheta, 2.0 * M_2_PI);
                                insert(0.0     , maxTheta);
                            }
                        }
                    }
                }
            }
            
            void consolidate(std::vector<Segments> segments) {
                for(unsigned int b=0; b < bands.size(); b++) {
                    size_t pendingInserts = std::accumulate(segments.cbegin(), segments.cend(), 0, [] (const size_t& sum, const Segment& s) {
                        return sum + s.insBuffers[worker_id][b].size();
                    });
                    
                    if (!pendingInserts)
                        continue;
                    
                    bands[b].requests.reserve(bands[b].requests.size() + pendingInserts);
                    
                    for(auto s : segments) {
                        auto& insBuf = s.insBuffers[worker_id][b];
                        bands[b].requests.insert(bands[b].requests.end(), insBuf.cbegin(), insBuf.cend());
                        insBuf.clear();
                    }

                    std::sort(bands[b].requests.begin(), bands[b].requests.end(), 
                              [] (const Request& a, const Request& b) {return a.start < b.start;});
                }
            }
        };
	
	
        vector<double> bandRadii = getBandRadii(n, R, seriesRatio);
        const unsigned int numBands = bandRadii.size()-1;
        bands.reserve(numBands);
        for(unsigned int b=0; b<bands.size(); b++) {
            bands.emplace_back( Band(bandRadii[b]) );
        }
        
        
        {
            Aux::Timer distrTimer;
            distrTimer.start();
        
            //Put points to bands
            for(index i=0; i<n; ++i) {
                const double& rad = radii[i];
                const double& theta = angles[i];
                
                const unsigned int bandIdx = std::distance(bandRadii.cbegin(), std::lower_bound(bandRadii.cbegin(), bandRadii.cend(), rad));

                const double e = std::exp(pt_radii);
                const double ie = 1.0 / e;
                const double cosh = 0.5*(e - ie);
                const double sinh = 0.5*(e + ie);
                
                bands[bandIdx].insertPoint(i, theta, cosh, sinh);
                
                // issue requests
                for(unsigned int b=bandIdx; b < numBands; ++b) {
                    bands[b].issueRequest(i, theta, cosh, sinh, coshR);
                }
            }
            distrTimer.stop();
            INFO("Distribution took ", distrTimer.elapsedMilliseconds(), " ms");
            
            Aux::Timer sortTimer;
            sortTimer.start();
            #pragma omp parallel for
            for(unsigned int b=0; b<bands.size(); b++) {
                bands[b].complete();
            }
            sortTimer.stop();
            
            INFO("Sorting of requests took ", sortTimer.elapsedMilliseconds(), " ms");
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
}
}


int main() {
    std::cout << "simple demonstration of NetworKit as a library\n";

    const uint64_t points = 10000;
    const double avgdeg = 50.0;
    const double alpha = 3.0;
    const double T = 0.0;
    
    NetworKit::HyperbolicGenerator gen(points, avgdeg, alpha, T);
    
    const double R = gen.getR();
    const auto pointsData = gen.generatePoints(points, R, alpha, T);
 
    for(unsigned int m=3; m; --m) {
        for(unsigned int r=0; r<3; r++) {
            Aux::Timer timer;
            timer.start();

            std::string key;
            switch(m) {
                case 1:
                    key = "Org";
                    gen.generateColdOrig(pointsData.first, pointsData.second, R);
                    break;

                case 2:
                    key = "Opt1";
                    gen.generateColdOrig(pointsData.first, pointsData.second, R);
                    break;

                case 3:
                    key = "Opt2";
                    gen.generateColdOrig(pointsData.first, pointsData.second, R);
                    break;
                    
                default:
                    abort();
            }
            timer.stop();
            std::cout << "Method: " << key << " took " << timer.elapsedMilliseconds() << "ms" << std::endl;
        }
    }
    

    return 0;

}

