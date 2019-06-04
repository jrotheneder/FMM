#ifndef FMM_UTILITY_H
#define FMM_UTILITY_H

#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 

#include "vector.hpp"

namespace fmm {

typedef std::complex<double> Complex; 

template<std::size_t d>
std::vector<PointSource_<d>> readSourcesFromFile(std::string filename) {

    std::ifstream ifile(filename); 
    std::string line;   

    std::vector<fmm::PointSource_<d>> sources; 
    std::array<double, d> pos; 
    double charge; 

    while(std::getline(ifile, line)) {

        std::istringstream ss(line);

        for(unsigned i = 0; i < d; ++i) {
            ss >> pos[i];   
        }
        ss >> charge;

        sources.push_back({{pos}, charge});   
    }

    return sources; 
}

template<typename T>
std::vector<T> vectorFromFile(std::string filename) {

    std::ifstream ifile(filename); 
    std::vector<T> outv; 

    T temp; 
    while(ifile >> temp) {
        outv.push_back(temp); 
    }

    return outv; 
}

template<typename T>
void vectorFromFile(std::string filename, T* dest_array) {

    std::ifstream ifile(filename); 
    T temp; 
    unsigned index = 0; 

    while(ifile >> temp) {
        dest_array[index++] = temp; 
    }
}

} // namespace fmm
#endif 
