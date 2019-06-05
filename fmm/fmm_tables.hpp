#ifndef FMM_TABLES_H
#define FMM_TABLES_H

#include <vector> 
#include <complex> 
#include <cmath> 

#include <gsl/gsl_sf_legendre.h> 
#include <gsl/gsl_sf_gamma.h>  

#include "fmm_general.hpp"

namespace fmm {
namespace tables {

template<typename T> struct PowTable {

    const T& x;            // base
    const unsigned nmax;   // max. exponent
    std::vector<T> table;  // stores [1, x, ... x^nmax] 

    PowTable(const T& x, unsigned nmax): x(x), nmax(nmax), table(nmax+1)  {

        T x_pow = 1; 
        for(unsigned n = 0; n <= nmax; ++n) {
            table[n] = x_pow; 
            x_pow *= x; 
        }
    }

    const T& operator()(unsigned n) const { return table[n]; }
};

struct BinomialTable {

    int max_order = 0;
    std::vector<double> table; 

    BinomialTable(int order = 0) { 
        refresh(order); 
        max_order = order; 
    }

    void refresh(int order) {

        if(order <= max_order) { return; }

        max_order = order;
        table.resize(((2*max_order+1)*(2*max_order+2))/2);

        for(int n = 0; n <= max_order; ++n) {
            for(int k = 0; k <= n; ++k) {
                table[getTableIndex(n,k)] = gsl_sf_choose(n,k);
            }
        }
    }

    unsigned getTableIndex(unsigned n, unsigned k) const {
        return (n*(n+1))/2 + k; 
    }

    const double& operator()(int n, int k) const {
        return table[getTableIndex(n,k)];
    }

};

struct YlmTable {

    const unsigned lmax; 
    const double theta; 
    const double phi; 

    std::vector<Complex> table;
     
    YlmTable(unsigned lmax, double theta, double phi): lmax(lmax), 
        theta(theta), phi(phi), table((lmax+1)*(lmax+1))  {

        using namespace std::complex_literals;
            
        double* result = new double[gsl_sf_legendre_array_n(lmax)];  
        gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SCHMIDT, lmax, 
            std::cos(theta), -1, result); 

        for(unsigned l = 0; l <= lmax; ++l) {

            //Handle m = 0 differently due to the Schmidt semi-normalization
            //convention used with gsl (comes the closest to what we want)
            table[getTableIndex(l, 0)] = result[gsl_sf_legendre_array_index(l, 0)];

            Complex exp_phi = std::exp(1i * phi); 
            Complex exp_m_phi = exp_phi; 

            for(unsigned m = 1; m <= l; ++m) {
                Complex Y_lm = result[gsl_sf_legendre_array_index(l, m)]
                     * exp_m_phi / std::sqrt(2);

                //If m = -|m| < 0, our convention implies M_l^{-|m|} = (M_l^|m|)^*
                table[getTableIndex(l, m)] = Y_lm; 
                table[getTableIndex(l, -m)] = std::conj(Y_lm); 

                exp_m_phi *= exp_phi;
            }
        }

        delete[] result;
    }

    unsigned getTableIndex(unsigned l, int m) const {
        return l*(l+1) + m; 
    }

    const Complex& operator()(unsigned l, int m) const {
        return table[getTableIndex(l, m)];
    }
};

struct YlmDerivTable {

    const unsigned lmax; 
    const double theta; 
    const double phi; 

    std::vector<Complex> ylm_table; // Y_l^m(theta, phi)
    std::vector<Complex> ylm_dtheta_table; // d/dθ Y_l^m(theta, phi)
    std::vector<Complex> exp_table; // e^(im *phi), m >=0 at position m.
     
    YlmDerivTable(unsigned lmax, double theta, double phi): lmax(lmax), 
        theta(theta), phi(phi), ylm_table((lmax+1)*(lmax+1)), 
        ylm_dtheta_table((lmax+1)*(lmax+1)), exp_table(lmax+1)    {

        using namespace std::complex_literals;
            
        Complex exp_phi = std::exp(1i * phi); 
        Complex exp_m_phi = 1; 

        for(unsigned l = 0; l <= lmax; ++l) {
            exp_table[l] = exp_m_phi;   
            exp_m_phi *= exp_phi; 
        }

        const unsigned gsl_array_size = gsl_sf_legendre_array_n(lmax); 
        double* result = new double[2 * gsl_array_size];  
        gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SCHMIDT, lmax, 
            std::cos(theta), -1, result, result + gsl_array_size);  
        
        for(unsigned l = 0; l <= lmax; ++l) {

            //Handle m = 0 differently due to the Schmidt semi-normalization
            //convention used with gsl (comes the closest to what we want)
            ylm_table[getTableIndex(l, 0)] 
                = result[gsl_sf_legendre_array_index(l, 0)];
            ylm_dtheta_table[getTableIndex(l, 0)] 
                = result[gsl_array_size + gsl_sf_legendre_array_index(l, 0)];


            for(unsigned m = 1; m <= l; ++m) {
                Complex Y_lm = result[gsl_sf_legendre_array_index(l, m)]
                     * exp_table[m] / std::sqrt(2);
                Complex Y_lm_dtheta 
                    = result[gsl_array_size + gsl_sf_legendre_array_index(l, m)]
                     * exp_table[m] / std::sqrt(2);

                //If m = -|m| < 0, our convention implies M_l^{-|m|} = (M_l^|m|)^*
                ylm_table[getTableIndex(l, m)] = Y_lm; 
                ylm_table[getTableIndex(l, -m)] = std::conj(Y_lm);  
                ylm_dtheta_table[getTableIndex(l, m)] = Y_lm_dtheta; 
                ylm_dtheta_table[getTableIndex(l, -m)] = std::conj(Y_lm_dtheta); 

            }
        }

        delete[] result;
    }

    unsigned getTableIndex(unsigned l, int m) const { return l*(l+1) + m; }

    const Complex& Y(unsigned l, int m) const {
        return ylm_table[getTableIndex(l,m)];    
    }

    const Complex& dYdtheta(unsigned l, int m) const {
        return ylm_dtheta_table[getTableIndex(l,m)];    
    }

    const Complex dYdphi(unsigned l, int m) const {

        using namespace std::complex_literals;
        return 1i * (double)m * ylm_table[getTableIndex(l,m)];    
    }
};

struct AlmTable {

    int max_order = 0;
    std::vector<double> table; 

    AlmTable(int order = 0) { 
        refresh(order); 
        max_order = order; 
    }

    void refresh(int order) {

        if(order <= max_order) { return; }

        max_order = order;
        table.resize((2*max_order+1)*(2*max_order+1));

        for(int l = 0; l <= 2*max_order; ++l) {
            for(int m = -l; m <= l; ++m) {
                table[getTableIndex(l,m)] = A_coeff(l, m);
            }
        }
    }

    unsigned getTableIndex(unsigned l, int m) const { return l*(l+1) + m; }

    const double& operator()(unsigned l, int m) const { 
        return table[getTableIndex(l,m) ]; 
    }

};

struct SignTable {

    const std::function <double (int, int)> signFunction;
    int max_order = 0;
    std::vector<double> table; 

    int offset; 

    SignTable(const std::function <double (int, int)> signFunction, 
            int order = 0): signFunction(signFunction) {

        refresh(order); 
        max_order = order; 
    }

    void refresh(int order) {

        if(order <= max_order) { return; }

        max_order = order;
        offset = 2 * max_order * (max_order + 1); 
        table.resize((2*max_order+1)*(2*max_order+1));

        for(int l = -max_order; l <= max_order; ++l) {
            for(int m = -max_order; m <= max_order; ++m) {
                table[getTableIndex(l,m)] = signFunction(l, m);
            }
        }
    }

    unsigned getTableIndex(int l, int m) const {
        return offset + l*(2*max_order+1) + m; 
    }

    const double& operator()(int l, int m) const {
        return table[getTableIndex(l,m)];
    }

};

} // namespace tables
} // namespace fmm

#endif
