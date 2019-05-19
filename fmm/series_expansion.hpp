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

    struct YlmTable;
    struct AlmTable;

    static AlmTable alm_table; 

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

    static double A_coeff(const int n, const int m);
    static double A_coeff_next(const int n, const int m, double A_nm); 

    virtual ~SeriesExpansion() {}

};

template<typename Vector, typename Source>
typename SeriesExpansion<Vector, Source, 3>::AlmTable 
    SeriesExpansion<Vector, Source, 3>::alm_table{};

template<typename Vector, typename Source>
struct SeriesExpansion<Vector, Source, 3>::YlmTable {

    const unsigned lmax; 
    const double theta; 
    const double phi; 

    std::vector<Complex> table;
     
    YlmTable(unsigned lmax, double theta, double phi): lmax(lmax), 
        theta(theta), phi(phi), table(((lmax + 1)*(lmax + 2))/2)  {

        using namespace std::complex_literals;
            
        double* result = new double[gsl_sf_legendre_array_n(lmax)];  
        gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SCHMIDT, lmax, 
            std::cos(theta), -1, result); 

        // TODO: consider also caching values for negative m 
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

        delete[] result;
    }

    Complex operator()(unsigned l, int m) const {

        // There are 1/2 * l * (l + 1)  coefficients M_k^m with k < l and m >= 0, 
        // => coefficient M_l^{-l} is at index 1/2 * l * (l + 1), coefficient 
        // M_l^m at 1/2 * l * (l + 1) + m. If m = -|m| < 0, the convention used in [1] 
        // implies that M_l^{-|m|} = (M_l^|m|)^*

        if(m >= 0) { return table.at((l*(l + 1))/2 + m); }
        else { return std::conj(table.at((l*(l + 1))/2 - m)); }
        //if(m >= 0) { return table[(l*(l + 1))/2 + m]; }
        //else { return std::conj(table[(l*(l + 1))/2 - m]); }
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
            double A_lm = A_coeff(l, -l);  
            for(int m = -l; m <= l; ++m) {
                table[coeff_index++] = A_lm;
                A_lm = A_coeff_next(l, m, A_lm);  
            }
        }
    }

    const double& operator()(unsigned l, int m) const {
        return table.at(l*(l+1) + m);
    }

};

template<typename Vector, typename Source>
SeriesExpansion<Vector, Source, 3>::SeriesExpansion(const Vector& center, 
        int order): order(order), center(center), 
        coefficients((1+order)*(1+order)) {
    
    if(order > alm_table.max_order) { alm_table.refresh(order); }
      
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

    return coefficients.at(n*(n+1) + m); 
    //return coefficients[n*(n+1) + m]; 
}

template<typename Vector, typename Source> 
const Complex& SeriesExpansion<Vector, Source, 3>::operator()(unsigned n, int m) 
        const {

    return coefficients.at(n*(n+1) + m); 
    //return coefficients[n*(n+1) + m]; 
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


// Implements the recurrence relation A_n^{m+1} = std::sqrt((n-m) / (n-m+1)) * A_n^m 
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::A_coeff_next(const int n, 
        const int m, double A_nm) {

//  // TODO remove
//  if(n < 0) {
//      throw std::logic_error("Argument n in A_coeff is expected "
//         "to be nonnegative."); 
//  }

    return std::sqrt((double)(n-m) / (n+m+1)) * A_nm;
}

} // namespace fmm

#endif
