#ifndef FMM_UTILITY_H
#define FMM_UTILITY_H

#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 
#include <complex> 

#include <gsl/gsl_sf_gamma.h>  

#include "vector.hpp"

namespace fmm {

typedef std::complex<double> Complex; 

// Implements the function A_n^m := (-1)^n / sqrt((n-m)! (n+m)!) [(5.23), 1]
double A_coeff(const int n, const int m) {

    if(n < 0) {
        throw std::logic_error("Argument n in A_coeff is expected "
            "to be nonnegative."); 
    }

    return (n % 2 ? -1 : 1) / std::sqrt(gsl_sf_fact(n-m) * gsl_sf_fact(n+m));
}

// Implements the function i^(|k| - |m| - |k-m|) for k, m integers
double sign_fun1(const int k, const int m) {

    const int exponent = std::abs(k) - std::abs(m) - std::abs(k-m); 
    switch(std::abs(exponent) % 4) {
        case 0 : return 1; 
        case 2 : return -1;  
    }

    throw std::logic_error("Exponent in sign_fun1() is not expected to "
        " be odd. Got input: k = " + std::to_string(k) + ", m = " 
        + std::to_string(m)+ ", exponent is " + std::to_string(exponent) + "\n"
        ); 

}

// Implements the function i^(|k-m|-|k|-|m|) for k, m integers
double sign_fun2(const int k, const int m) {

    const int exponent = std::abs(k-m) - std::abs(k) - std::abs(m); 
    switch(std::abs(exponent) % 4) {
        case 0  : return 1; 
        case 2  : return -1;  
    }

    throw std::logic_error("Exponent in sign_fun2() is not expected to "
        " be odd. Got input: k = " + std::to_string(k) + ", m = " 
        + std::to_string(m)+ ", exponent is " + std::to_string(exponent) + "\n"
        ); 

}

// Implements the function i^(|m|-|m-k|-|k|) for k, m integers
double sign_fun3(const int k, const int m) {

    const int exponent = std::abs(m) - std::abs(m-k) - std::abs(k); 
    switch(std::abs(exponent) % 4) {
        case 0 : return 1; 
        case 2 : return -1;  
    }

    throw std::logic_error("Exponent in sign_fun3() is not expected to "
        " be odd. Got input: k = " + std::to_string(k) + ", m = " 
        + std::to_string(m)+ ", exponent is " + std::to_string(exponent) + "\n"
        ); 
}

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

template<std::size_t d>
void sourceLocationsToFile(std::vector<PointSource_<d>> sources, 
        std::string filename, std::string sep = "\n") {

    std::ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific;

    for(auto s : sources) {
        for(unsigned i = 0; i < d; ++i) {
            ofile << s[i] << "\t"; 
        }
        ofile << "\n"; 
    }

    ofile.close();
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
