#include <iostream>
#include <NetworKit/generators/HyperbolicGenerator.h>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "n avgdeg alpha" << std::endl;
        return -1;
    }
    
    uint64_t n = atoll(argv[1]);
    double avgdeg = atof(argv[2]);
    double alpha = atof(argv[3]);
    
    
    const double T = 0.0;
    
    NetworKit::HyperbolicGenerator gen(n, avgdeg, alpha, T);
 
    const double R = gen.getR();
    std::cout << "n=" << n << " avgdeg=" << avgdeg << " alpha=" << alpha << " R=" << R << std::endl;

    const auto points = gen.generatePoints(n, R, alpha, T);
 
    auto bandRadii = NetworKit::HyperbolicGenerator::getBandRadii(n, R);
    std::vector<uint64_t> counts(bandRadii.size());
    
    for(const auto& r : points.second) {
        auto lower = std::lower_bound(bandRadii.begin(), bandRadii.end(), r);
        counts.at( std::distance(bandRadii.begin(), lower) )++;
    }
    
    auto prob = [&] (double r) {return (std::cosh(r * alpha) - 1.0) / (std::cosh(R * alpha) - 1.0);};
    
    
    uint64_t sum = 0;
    for(unsigned int i=0; i < counts.size(); i++) {
        double expect = 0.0;
        if (i)
            expect = prob(bandRadii[i]) - prob(bandRadii[i-1]);
        
        std::cout << std::setw(4) << i 
                  << std::setw(10) << bandRadii[i] 
                  << " | "

                  << std::setw(10) << counts[i]
                  << std::setw(10) << std::round(expect * n)
                  
                  << " | "
                  
                  << std::setw(15) << (double(counts[i]) / n)
                  << std::setw(15) << (expect)

                  << " | "
                  
                  << std::setw(20) << n
                  << std::setw(8) << R
                  << std::setw(3) << alpha
        << std::endl;
        sum += counts[i];
    }
    
    assert(sum == n);
    
    
    return 0;
}
    

