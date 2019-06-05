#ifndef SERIES_EXPANSION_H
#define SERIES_EXPANSION_H

#include <stdexcept> 

#include <gsl/gsl_sf_legendre.h> 
#include <gsl/gsl_sf_gamma.h>  

#include "fmm_tables.hpp"

namespace fmm {

template<std::size_t d>
struct SeriesExpansion {
    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};
 

/******************************************************************************/
/*                      2D Series Expansion Implementation                    */
/******************************************************************************/
template<> struct SeriesExpansion<2> {

    int order; // Order of the expansion
    Complex center; // Center of the expansion
    std::vector<Complex> coefficients; 
    static tables::BinomialTable binomial_table; 

    SeriesExpansion(): order(), center(), coefficients() {}
    SeriesExpansion(const Vector_<2>& center, int order);

    SeriesExpansion& operator+=(const SeriesExpansion& rhs);
    Complex& operator()(unsigned n); // Access coefficient a_n
    const Complex& operator()(unsigned n) const;

    virtual double evaluatePotential(const Vector_<2>& eval_point) const = 0;
    virtual Vector_<2> evaluateForcefield(const Vector_<2>& eval_point) const = 0;

    friend std::ostream& operator<<(std::ostream& o, const SeriesExpansion& me);

    virtual ~SeriesExpansion() {}
};

typename tables::BinomialTable SeriesExpansion<2>::binomial_table;

SeriesExpansion<2>::SeriesExpansion(const Vector_<2>& center, 
        int order): order(order), center({center[0], center[1]}), 
        coefficients(order + 1) {

    if(order > binomial_table.max_order) { binomial_table.refresh(order); }
}


SeriesExpansion<2>& SeriesExpansion<2>::
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

 
Complex& SeriesExpansion<2>::operator()(unsigned n) {

    return coefficients[n]; 
}

 
const Complex& SeriesExpansion<2>::operator()(unsigned n) 
        const {

    return coefficients[n]; 
}

 
std::ostream& operator<<(std::ostream& o, const SeriesExpansion<2>& se) {

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

template<> struct SeriesExpansion<3> {

    static tables::AlmTable alm_table; 
    static tables::SignTable sign_fun1_table; 
    static tables::SignTable sign_fun2_table; 
    static tables::SignTable sign_fun3_table; 

    int order; // Order of the expansion
    Vector_<3> center;  // Center of the expansion
    std::vector<Complex> coefficients; // ME coefficients

    SeriesExpansion(): order(), center(), coefficients() {}
    SeriesExpansion(const Vector_<3>& center, int order);

    SeriesExpansion& operator+=(const SeriesExpansion& rhs);
    Complex& operator()(unsigned n, int m); // Access coefficients M_n^m
    const Complex& operator()(unsigned n, int m) const;

    virtual double evaluatePotential(const Vector_<3>& eval_point) const = 0;
    virtual Vector_<3> evaluateForcefield(const Vector_<3>& eval_point) const = 0;

    template<typename V, typename S>
    friend std::ostream& operator<<(std::ostream& o, const SeriesExpansion& me);

    virtual ~SeriesExpansion() {}

};

typename tables::AlmTable SeriesExpansion<3>::alm_table; 
typename tables::SignTable SeriesExpansion<3>::sign_fun1_table{sign_fun1, 0};
typename tables::SignTable SeriesExpansion<3>::sign_fun2_table{sign_fun2, 0};
typename tables::SignTable SeriesExpansion<3>::sign_fun3_table{sign_fun3, 0};



SeriesExpansion<3>::SeriesExpansion(const Vector_<3>& center, int order): 
        order(order), center(center), coefficients((1+order)*(1+order)) {
    
    if(order > alm_table.max_order) { 
        alm_table.refresh(order); 
        sign_fun1_table.refresh(order); 
        sign_fun2_table.refresh(order); 
        sign_fun3_table.refresh(order); 
    }
}

SeriesExpansion<3>& SeriesExpansion<3>::operator+=(const SeriesExpansion& rhs) {

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

Complex& SeriesExpansion<3>::operator()(unsigned n, int m) {

    // There are n*n coefficients M_k^m with k < n, => coefficient M_n^{-n} starts
    // at index n*n, coefficient M_n^m at n*n + n + m

    return coefficients[n*(n+1) + m]; 
}
 
const Complex& SeriesExpansion<3>::operator()(unsigned n, int m) const {

    return coefficients[n*(n+1) + m]; 
}

std::ostream& operator<<(std::ostream& o, const SeriesExpansion<3>& se) {

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

} // namespace fmm

#endif
