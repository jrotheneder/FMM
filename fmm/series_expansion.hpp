#ifndef SERIES_EXPANSION_H
#define SERIES_EXPANSION_H

namespace fmm {

typedef std::complex<double> Complex; 
static const double PI = boost::math::constants::pi<double>();
constexpr auto sphericalHarmonicY = boost::math::spherical_harmonic<double, double>;
constexpr auto factorial = boost::math::factorial<double>;

Complex YLM(unsigned l, int m, double theta, double phi) {
//  return std::sqrt(4 * PI / (2 *l + 1)) * (m > 0 ? sphericalHarmonicY(l,m,theta,phi) 
//               : std::conj(sphericalHarmonicY(l,-m,theta,phi)));
    return std::sqrt(4 * PI / (2 *l + 1)) * (m > 0 ? sphericalHarmonicY(l,m,theta,phi) 
                 : (m % 2 ? -1. : 1.) * sphericalHarmonicY(l,m,theta,phi));
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

    int order; 
    Vector center; // Center of the expansion
    std::vector<Complex> coefficients; // ME coefficients

    SeriesExpansion() {}
    SeriesExpansion(const Vector& center, int order): order(order), center(center),
        coefficients((1+order)*(1+order)) {}

    SeriesExpansion& operator+=(const SeriesExpansion& rhs);
    Complex& operator()(unsigned n, int m); // Access coefficients M_n^m
    const Complex& operator()(unsigned n, int m) const;

    virtual double evaluatePotential(const Vector& eval_point) const = 0;
    virtual Vector evaluateForcefield(const Vector& eval_point) const = 0;

    template<typename V, typename S>
    friend std::ostream& operator<<(std::ostream& o, const SeriesExpansion& me);

    static double A_coeff(const int n, const int m);

    virtual ~SeriesExpansion() {}

};

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

    return coefficients.at(n*(n+1) + m); // TODO possibly exch. at() with []
}

template<typename Vector, typename Source> 
const Complex& SeriesExpansion<Vector, Source, 3>::operator()(unsigned n, int m) 
        const {

    return coefficients.at(n*(n+1) + m); // TODO possibly exch. at() with []
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


// Implements the function (-1)^n / sqrt((n-m)! (n+m)!) [(5.23), 1]
template<typename Vector, typename Source> 
double SeriesExpansion<Vector, Source, 3>::A_coeff(const int n, const int m) {

    // TODO remove
    if(n < 0) {
        throw std::logic_error("Argument n in A_coeff is expected "
            "to be nonnegative."); 
    }

    return (n % 2 ? -1 : 1) / std::sqrt(factorial(n-m) * factorial(n+m));
}

} // namespace fmm

#endif
