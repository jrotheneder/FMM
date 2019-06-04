#ifndef SERIES_EXPANSION_H
#define SERIES_EXPANSION_H

#include <stdexcept> 

#include <gsl/gsl_sf_legendre.h> 
#include <gsl/gsl_sf_gamma.h>  

#include "fmm_utility.hpp"

namespace fmm {

template<typename Vector, typename Source, std::size_t d>
struct SeriesExpansion {};


/******************************************************************************/
/*                      2D Series Expansion Implementation                    */
/******************************************************************************/
template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 2> {

    struct BinomialTable;  // table of binomial coefficients

    template<typename T>
    struct PowTable;       // table of powers of a given quantity of type T

    unsigned order; 
    Complex center; // Center of the expansion
    std::vector<Complex> coefficients; 
    static BinomialTable binomial_table; 

    SeriesExpansion(): order(), center(), coefficients() {}
    SeriesExpansion(const Vector& center, unsigned order);

    SeriesExpansion& operator+=(const SeriesExpansion& rhs);
    Complex& operator()(unsigned n); // Access coefficient a_n
    const Complex& operator()(unsigned n) const;

    virtual double evaluatePotential(const Vector& eval_point) const = 0;
    virtual Vector evaluateForcefield(const Vector& eval_point) const = 0;

    template<typename V, typename S>
    friend std::ostream& operator<<(std::ostream& o, const SeriesExpansion& me);

    virtual ~SeriesExpansion() {}
};

template<typename Vector, typename Source>
template<typename T>
struct SeriesExpansion<Vector, Source, 2>::PowTable {

    const T& x;            // base
    const unsigned nmax;   // we store [1, x, ... x^nmax] 
    std::vector<T> table; 

    PowTable(const T& x, unsigned nmax): x(x), nmax(nmax), table(nmax+1)  {

        T x_pow = 1; 
        for(unsigned n = 0; n <= nmax; ++n) {
            table[n] = x_pow; 
            x_pow *= x; 
        }
    }

    const T& operator()(unsigned n) const { return table[n]; }
};

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 2>::BinomialTable
    SeriesExpansion<Vector, Source, 2>::binomial_table;

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 2>::BinomialTable {

    unsigned max_order = 0;
    std::vector<double> table; 

    BinomialTable(unsigned order = 0) { 
        refresh(order); 
        max_order = order; 
    }

    void refresh(unsigned order) {

        if(order <= max_order) { return; }

        max_order = order;
        table.resize(((2*max_order+1)*(2*max_order+2))/2);

        for(unsigned n = 0; n <= max_order; ++n) {
            for(unsigned k = 0; k <= n; ++k) {
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


template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 2>::SeriesExpansion(const Vector& center, 
        unsigned order): order(order), center({center[0], center[1]}), 
        coefficients(order + 1) {

    if(order > binomial_table.max_order) { binomial_table.refresh(order); }
}

template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 2>& SeriesExpansion<Vector, Source, 2>::
        operator+=(const SeriesExpansion& rhs) {

    if(coefficients.size() != rhs.coefficients.size() 
            || this->center != rhs.center) {

        throw std::logic_error("Series Expansions can only be added if "
            "they are of the same order and w.r.t. the same origin."); 
    }

    std::transform (
        coefficients.begin(), coefficients.end(), rhs.coefficients.begin(), 
        coefficients.begin(), std::plus<Complex>()
    );

    return *this;
}

template<typename Vector, typename Source> 
Complex& SeriesExpansion<Vector, Source, 2>::operator()(unsigned n) {

    return coefficients[n]; 
}

template<typename Vector, typename Source> 
const Complex& SeriesExpansion<Vector, Source, 2>::operator()(unsigned n) 
        const {

    return coefficients[n]; 
}

template<typename Vector, typename Source> 
std::ostream& operator<<(std::ostream& o, 
        const SeriesExpansion<Vector, Source, 2>& se) {

    o << "n\tRe a_n\tIm a_n\n";
    for(int n = 0; n <= se.order; ++n) { 
        o << n << "\t" << se.coefficients[n].real() << "\t"
            << se.coefficients[n].imag() << "\n";
    }
    o << "\n";

   return o; 
}

/******************************************************************************/
/*                      3D Series Expansion Implementation                    */
/******************************************************************************/
template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3> {

    // Lookup tables for frequently needed quantities
    struct YlmTable;  // spherical harmonics
    struct YlmDerivTable;  // derivatives of spherical harmonics
    struct AlmTable;  // A_l^m := (-1)^l/sqrt((l-m)!(l+m)!)
    struct SignTable; // sign patterns like I^(|m|-|k-m|-|k|) arising in M2M, M2L, L2L

    static AlmTable alm_table; 
    static SignTable sign_fun1_table; 
    static SignTable sign_fun2_table; 
    static SignTable sign_fun3_table; 

    unsigned order; 
    Vector center; // Center of the expansion
    std::vector<Complex> coefficients; // ME coefficients

    SeriesExpansion(): order(), center(), coefficients() {}
    SeriesExpansion(const Vector& center, unsigned order);

    SeriesExpansion& operator+=(const SeriesExpansion& rhs);
    Complex& operator()(unsigned n, int m); // Access coefficients M_n^m
    const Complex& operator()(unsigned n, int m) const;

    virtual double evaluatePotential(const Vector& eval_point) const = 0;
    virtual Vector evaluateForcefield(const Vector& eval_point) const = 0;

    template<typename V, typename S>
    friend std::ostream& operator<<(std::ostream& o, const SeriesExpansion& me);

    // Various coefficients needed in the M2M, M2L, L2L operators
    static double A_coeff(const int n, const int m);
    static double A_coeff_next(const int n, const int m, double A_nm); 
    static double sign_fun1(const int k, const int m);
    static double sign_fun2(const int k, const int m);
    static double sign_fun3(const int k, const int m);

    virtual ~SeriesExpansion() {}

};

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 3>::AlmTable
    SeriesExpansion<Vector, Source, 3>::alm_table;

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 3>::SignTable
    SeriesExpansion<Vector, Source, 3>::sign_fun1_table{sign_fun1, 0};

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 3>::SignTable
    SeriesExpansion<Vector, Source, 3>::sign_fun2_table{sign_fun2, 0};

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 3>::SignTable
    SeriesExpansion<Vector, Source, 3>::sign_fun3_table{sign_fun3, 0};

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::YlmTable {

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

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::YlmDerivTable {

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

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::AlmTable {

    unsigned max_order = 0;
    std::vector<double> table; 

    AlmTable(unsigned order = 0) { 
        refresh(order); 
        max_order = order; 
    }

    void refresh(unsigned order) {

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

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::SignTable {

    const std::function <double (int, int)> signFunction;
    unsigned max_order = 0;
    std::vector<double> table; 

    int offset; 

    SignTable(const std::function <double (int, int)> signFunction, 
            unsigned order = 0): signFunction(signFunction) {

        refresh(order); 
        max_order = order; 
    }

    void refresh(unsigned order) {

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

template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 3>::SeriesExpansion(const Vector& center, 
        unsigned order): order(order), center(center), 
        coefficients((1+order)*(1+order)) {
    
    if(order > alm_table.max_order) { 
        alm_table.refresh(order); 
        sign_fun1_table.refresh(order); 
        sign_fun2_table.refresh(order); 
        sign_fun3_table.refresh(order); 
    }
}

template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 3>& SeriesExpansion<Vector, Source, 3>::
        operator+=(const SeriesExpansion& rhs) {

    if(coefficients.size() != rhs.coefficients.size() 
            || this->center != rhs.center) {

        throw std::logic_error("Series Expansions can only be added if "
            "they are of the same order and w.r.t. the same origin."); 
    }

    std::transform (
        coefficients.begin(), coefficients.end(), rhs.coefficients.begin(), 
        coefficients.begin(), std::plus<Complex>()
    );

    return *this;
}

template<typename Vector, typename Source> 
Complex& SeriesExpansion<Vector, Source, 3>::operator()(unsigned n, int m) {

    // There are n*n coefficients M_k^m with k < n, => coefficient M_n^{-n} starts
    // at index n*n, coefficient M_n^m at n*n + n + m

    return coefficients[n*(n+1) + m]; 
}

template<typename Vector, typename Source> 
const Complex& SeriesExpansion<Vector, Source, 3>::operator()(unsigned n, int m) 
        const {

    return coefficients[n*(n+1) + m]; 
}

template<typename Vector, typename Source> 
std::ostream& operator<<(std::ostream& o, 
        const SeriesExpansion<Vector, Source, 3>& se) {

    unsigned coeff_index = 0; 

    o << "n\tm\tRe M_l^m\tIm M_l^m\n";

    for(int n = 0; n <= se.order; ++n) { 
        for(int m = -n; m <= n; ++m) {
            o << n << "\t" << m << "\t" 
                << se.coefficients[coeff_index].real() << "\t\t"
                << se.coefficients[coeff_index].imag() << "\n";
            ++coeff_index;
        }
    }
    o << "\n";
   return o; 
}


// Implements the function A_n^m := (-1)^n / sqrt((n-m)! (n+m)!) [(5.23), 1]
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::A_coeff(const int n, 
        const int m) {

    if(n < 0) {
        throw std::logic_error("Argument n in A_coeff is expected "
            "to be nonnegative."); 
    }

    return (n % 2 ? -1 : 1) / std::sqrt(gsl_sf_fact(n-m) * gsl_sf_fact(n+m));
}

// Implements the function i^(|k| - |m| - |k-m|) for k, m integers
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::sign_fun1(
        const int k, const int m) {

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
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::sign_fun2(
        const int k, const int m) {

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
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::sign_fun3(
        const int k, const int m) {

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

} // namespace fmm

#endif
