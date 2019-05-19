#ifndef SERIES_EXPANSION_H
#define SERIES_EXPANSION_H

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <gsl/gsl_sf_legendre.h> 

namespace fmm {

typedef std::complex<double> Complex; 
static const double PI = boost::math::constants::pi<double>();
constexpr auto sphericalHarmonicY = boost::math::spherical_harmonic<double, double>;
constexpr auto factorial = boost::math::factorial<double>;

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

// std::beta requires C++17
double binomial(std::size_t n, std::size_t k) { // TODO consider a lookup table if slow
    return 1 / ((n+1) * std::beta(n-k+1, k+1)); // TODO investigate accuracy
                                                // TODO signal if overflow
}

template<typename Vector, typename Source, std::size_t d>
struct SeriesExpansion {};

// TODO: unify 2d expansion interface
template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 2> {
};

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3> {


    // Lookup tables for frequently needed values
    struct YlmTable;  // spherical harmonics
    struct AlmTable;  // A_l^m := (-1)^l/sqrt((l-m)!(l+m)!)
    struct SignTable; // for sign patterns like I^(|m|-|k-m|-|k|) [aka sign_fun2(k, m)]

    static AlmTable alm_table; 
    static SignTable sign_fun1_table; 
    static SignTable sign_fun2_table; 
    static SignTable sign_fun3_table; 

    int order; 
    Vector center; // Center of the expansion
    std::vector<Complex> coefficients; // ME coefficients

    SeriesExpansion(): order(), center(), coefficients() {}
    SeriesExpansion(const Vector& center, int order);

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
     
    // TODO clean up
    YlmTable(unsigned lmax, double theta, double phi): lmax(lmax), 
//      theta(theta), phi(phi), table(((lmax + 1)*(lmax + 2))/2)  {
        theta(theta), phi(phi), table((lmax+1)*(lmax+1))  {

        using namespace std::complex_literals;
            
        double* result = new double[gsl_sf_legendre_array_n(lmax)];  
        gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SCHMIDT, lmax, 
            std::cos(theta), -1, result); 

        /*
        unsigned coeff_index = 0; // index of next coefficient to be computed
        for(std::size_t l = 0; l <= lmax; ++l) {

            //Handle m = 0 differently due to the Schmidt semi-normalization
            //convention used with gsl (comes the closest to what we want)
            table[coeff_index++] = result[gsl_sf_legendre_array_index(l, 0)];

            Complex exp_phi = std::exp(1i * phi); 
            Complex exp_m_phi = 1; 

            for(std::size_t m = 1; m <= l; ++m) {
                 table[coeff_index++] = result[gsl_sf_legendre_array_index(l, m)]
                     * (exp_m_phi *= exp_phi) / std::sqrt(2);
            }
        }
        */
        for(std::size_t l = 0; l <= lmax; ++l) {

            //Handle m = 0 differently due to the Schmidt semi-normalization
            //convention used with gsl (comes the closest to what we want)
            table[getTableIndex(l, 0)] = result[gsl_sf_legendre_array_index(l, 0)];

            Complex exp_phi = std::exp(1i * phi); 
            Complex exp_m_phi = exp_phi; 

            for(std::size_t m = 1; m <= l; ++m) {
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
        // There are l*l  coefficients M_k^m with k < l and m >= 0, 
        // => coefficient M_l^{-l} is at index l*l, coefficient 
        // M_l^m at l*l + l + m. 
        return l*(l+1) + m; 
    }

    const Complex& operator()(unsigned l, int m) const {
        //return table.at(getTableIndex(l, m));
        return table[getTableIndex(l, m)];
    }
};

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::AlmTable {

    int max_order;
    std::vector<double> table; 

    AlmTable(int order = 0): max_order(order) {
        refresh(order); 
    }

    void refresh(int order) {

        max_order = 2 * order;
        table.resize((max_order+1)*(max_order+1));

        unsigned coeff_index = 0;
        for(int l = 0; l <= max_order; ++l) {
            for(int m = -l; m <= l; ++m) {
                table[coeff_index++] = A_coeff(l, m);
            }
        }
    }

    const double& operator()(int l, int m) const {
        //return table.at(l*(l+1) + m);
        return table[l*(l+1) + m];
    }

};

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::SignTable {

    const std::function <double (int, int)> signFunction;
    int max_order;
    std::vector<double> table; 

    int offset; 

    SignTable(const std::function <double (int, int)> signFunction, 
            int order = 0): signFunction(signFunction), max_order(order) {

        refresh(order); 
    }

    void refresh(int order) {

        max_order = order;
        offset = 2 * max_order * (max_order + 1); 
        table.resize((2*max_order+1)*(2*max_order+1));

        unsigned coeff_index = 0;
        for(int l = -max_order; l <= max_order; ++l) {
            for(int m = -max_order; m <= max_order; ++m) {
                table[coeff_index++] = signFunction(l, m);
            }
        }
    }

    const double& operator()(int l, int m) const {
        //return table.at(offset + l*(2*max_order+1) + m);
        return table[offset + l*(2*max_order+1) + m];
    }

};

template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 3>::SeriesExpansion(const Vector& center, 
        int order): order(order), center(center), 
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

    // TODO exception this
    assert(coefficients.size() == rhs.coefficients.size());
    assert(this->center == rhs.center);

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

    //return coefficients.at(n*(n+1) + m); 
    return coefficients[n*(n+1) + m]; 
}

template<typename Vector, typename Source> 
const Complex& SeriesExpansion<Vector, Source, 3>::operator()(unsigned n, int m) 
        const {

    //return coefficients.at(n*(n+1) + m); 
    return coefficients[n*(n+1) + m]; 
}

template<typename Vector, typename Source> 
std::ostream& operator<<(std::ostream& o, SeriesExpansion<Vector, Source, 3>& me) {

    unsigned coeff_index = 0; // index of next coefficient to be computed

    o << "n\tm\tRe M_l^m\tIm M_l^m\n";

    for(int n = 0; n <= me.order; ++n) { 
        for(int m = -n; m <= n; ++m) {
            o << n << "\t" << m << "\t" 
                << me.coefficients[coeff_index].real() << "\t\t"
                << me.coefficients[coeff_index].imag() << "\n";
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

    return (n % 2 ? -1 : 1) / std::sqrt(factorial(n-m) * factorial(n+m));
}


// TODO depreacted in favor of lookup tables
// Implements the recurrence relation A_n^{m+1} = std::sqrt((n-m) / (n-m+1)) * A_n^m 
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::A_coeff_next(const int n, 
        const int m, double A_nm) {

    return std::sqrt((double)(n-m) / (n+m+1)) * A_nm;
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
