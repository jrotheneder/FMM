#ifndef FMM_UTILITY_H
#define FMM_UTILITY_H

#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "vector.hpp"

namespace fmm {

typedef std::complex<double> Complex; 
static const double PI = boost::math::constants::pi<double>();
constexpr auto sphericalHarmonicY = boost::math::spherical_harmonic<double, double>;
constexpr auto factorial = boost::math::factorial<double>;

// TODO remove
Complex YLM(unsigned l, int m, double theta, double phi) {
//  return std::sqrt(4 * PI / (2 *l + 1)) * (m > 0 ? sphericalHarmonicY(l,m,theta,phi) 
//               : std::conj(sphericalHarmonicY(l,-m,theta,phi)));
    if(m >= 0) {
        return std::sqrt(4 * PI / (2 *l + 1)) * sphericalHarmonicY(l, m, theta, phi);
    }
    else {
        return std::sqrt(4 * PI / (2 *l + 1)) * (m % 2 ? -1. : 1.) * 
            sphericalHarmonicY(l, m, theta, phi);
    }
}

Complex YLM_deriv_theta(unsigned l, int m, double theta, double phi) {
    using namespace std::complex_literals;

    if(m >= 0) {
        return -std::sqrt(4 * PI / (2 *l + 1)) * 
            (
            std::sqrt((l + m) * (l - m + 1)) * std::exp(1i * phi) 
                * sphericalHarmonicY(l, m - 1, theta, phi)
            + (double)m * sphericalHarmonicY(l, m, theta, phi) / std::tan(theta)
            );
    }
    else {
        return (m % 2 ? -1. : 1.) * -std::sqrt(4 * PI / (2 *l + 1)) * 
            (
            std::sqrt((l + m) * (l - m + 1)) * std::exp(1i * phi) 
            * sphericalHarmonicY(l, m - 1, theta, phi)
            + (double)m * sphericalHarmonicY(l, m, theta, phi) / std::tan(theta)
            );
    }
}

Complex YLM_deriv_phi(unsigned l, int m, double theta, double phi) {
    using namespace std::complex_literals;

    return YLM(l, m, theta, phi) * (double)m * 1i; 
}

// TODO use boost's binomial / lookup tables for this
double binomial(std::size_t n, std::size_t k) { 
    return 1 / ((n+1) * std::beta(n-k+1, k+1)); 
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
