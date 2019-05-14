#include <complex> 
#include <type_traits> 
#include <numeric> 
#include <algorithm> 

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/constants/constants.hpp>

namespace fmm {

static const double PI = boost::math::constants::pi<double>();

// std::beta requires C++17
double binomial(std::size_t n, std::size_t k) { // TODO consider a lookup table if slow
    return 1 / ((n+1) * std::beta(n-k+1, k+1)); // TODO investigate accuracy
                                                // TODO signal if overflow
}

template<typename Vector, typename Source, std::size_t d>
struct MultipoleExpansion {

    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source> 
struct MultipoleExpansion<Vector, Source, 2> {

    using Complex = std::complex<double>;
    using ME = MultipoleExpansion; 

    std::vector<Complex> coefficients; // ME coefficients
    Complex center; // Center of the expansion
    double Q; // total source strength of all sources contained in box

    MultipoleExpansion() {}
    MultipoleExpansion(const Vector& center_vec, std::size_t order, 
            std::vector<Source>& sources): coefficients(order), 
            center({center_vec[0], center_vec[1]}), Q(0) {
        
        // Compute series expansion coefficients a_1 through a_order: 
        for(std::size_t i = 0; i < sources.size(); ++i) {

            Complex z{sources[i].position[0], sources[i].position[1]};
            Complex z_rel = z - center; // Express z in box-local coords

            for(std::size_t j = 1; j <= order; ++j) {
                coefficients[j-1] -=  // TODO: consider better pow
                    sources[i].q * std::pow(z_rel, j) / (double)j; 
            }
            Q += sources[i].q;
        }
    }

    // Construct expansion from other expansions by shifting
    MultipoleExpansion(const Vector& center_vec, std::vector<ME*>& expansions): 
            coefficients(expansions[0]->coefficients.size()), 
            center({center_vec[0], center_vec[1]}), Q(0) {
        
        for(ME* me : expansions) {
            Q += me->Q;     

            assert(me->coefficients.size() == coefficients.size()); // TODO remove
        
            Complex shift_vec = me->center - this->center; 
            std::vector<Complex> shifted_coefficients 
                = me->shift(shift_vec);

            std::transform (coefficients.begin(), coefficients.end(), 
                shifted_coefficients.begin(), coefficients.begin(), 
                std::plus<Complex>()
            );
        }
    }

    // Shift is the vector (complex number) from the new to the old center.
    std::vector<Complex> shift(const Complex& shift) const {

        std::vector<Complex> shifted_coefficients(coefficients.size()); 

        // TODO: precompute pows, better binomial if slow
        for(std::size_t l = 1; l <= shifted_coefficients.size() - 1; ++l) {
            shifted_coefficients[l-1] = -Q * pow(shift, l) / (double)l; 
            for(std::size_t k = 1; k <= l; ++k) {
                shifted_coefficients[l-1] += coefficients[k-1] * pow(shift, l-k) 
                    * binomial(l-1, k-1);  // [(4.15), 1]  
            }
        }

        return shifted_coefficients; 
    }

    double evaluatePotential(const Vector& eval_point) const {

        Complex z{eval_point[0], eval_point[1]}; // get complex repr.
        Complex z_rel = z - center; 

        Complex result = Q * log(z_rel); 
        
        // TODO Horner scheme this
        for(int j = coefficients.size(); j >= 1; --j) { // a_j is stored in coeff[j-1]  
            result += coefficients[j-1] / pow(z_rel, j);
        }

        // Result contains the eval. of the (inverse) complex power series for 
        // φ = q log(z-z0) (i.e. the complex 2d gravitational potential) summed 
        // over all source locations z0. 
        // Return -1 times this for the electrostatic potential
        return -result.real(); 
    } 

    Vector evaluateForcefield(const Vector& eval_point) const { 

        Complex z{eval_point[0], eval_point[1]}; // get complex repr.
        Complex z_rel = z - center; 

        Complex result = Q / z_rel;

        // TODO Horner scheme this
        for(int j = coefficients.size(); j >= 1; --j) { // a_j is stored in coeff[j-1]  
            result -= (double)j * coefficients[j-1] / pow(z_rel, j+1);
        }

        // Result contains the eval. of *the derivative* φ' of the (inverse) complex 
        // power series for φ = q log(z-z0) (i.e. the complex 2d gravitational potential) 
        // summed over all source locations z0. The force field is then given by 
        // \vec F - ∇φ = (-∂x φ, -∂y φ) = (-Re φ, Im φ). For the electric field,
        // return -1 times that, i.e. (Re φ, -Im φ).
        return {{result.real(), -result.imag()}}; 
    }

};


// 3-D Implementation
template<typename Vector, typename Source> 
struct MultipoleExpansion<Vector, Source, 3> {

    using Complex = std::complex<double>;
    using ME = MultipoleExpansion; 
    static constexpr auto sphericalHarmonicY 
        = boost::math::spherical_harmonic<double, double>;

    int order; 
    std::vector<Complex> coefficients; // ME coefficients
    Vector center; // Center of the expansion

    MultipoleExpansion() {}
    MultipoleExpansion(const Vector& center_vec, int order, 
            std::vector<Source>& sources);
    MultipoleExpansion(const Vector& center_vec, std::vector<ME*>& expansions);

    std::vector<double> shift(const Vector& shift) const { return {}; }

    double evaluatePotential(const Vector& eval_point) const;
    Vector evaluateForcefield(const Vector& eval_point) const;

    friend std::ostream& operator<<(std::ostream& o, const ME& me) {

        unsigned coeff_index = 0; // index of next coefficient to be computed
        for(int n = 0; n <= me.order; ++n) { 
            for(int m = -n; m <= n; ++m) {
                o << n << "\t" << m << "\t" 
                    << me.coefficients[coeff_index].real() << "\t"
                    << me.coefficients[coeff_index].imag() << "\n"; 
                ++coeff_index;
            }
        }
        o << "\n";
       return o; 
    }
};

template<typename Vector, typename Source> 
MultipoleExpansion<Vector, Source, 3>::MultipoleExpansion(
        const Vector& center, int order, std::vector<Source>& sources): 
        order(order), coefficients((1+order)*(1+order)), center(center) {

    for(std::size_t i = 0; i < sources.size(); ++i) {

        Vector rel_pos = (sources[i].position - center).toSpherical(); 

        double r = rel_pos[0];  // Distance of source to center
        double theta = rel_pos[1];  
        double phi = rel_pos[2];  


        // compute M_n^m according to [1, (5.15)] However (!)  we use a diff. 
        // convention (Mathematica, Boost, Wikipedia) than [1] for the harmonics)

        // TODO smarter pow
        // TODO exploit possible symmetry in coefficients? 
        unsigned coeff_index = 0; // index of next coefficient to be computed
        for(int n = 0; n <= order; ++n) {  
            for(int m = -n; m <= n; ++m) {  
                                                
                coefficients[coeff_index++] +=  
                    sources[i].sourceStrength() * pow(r, n) * 
                         std::conj(sphericalHarmonicY(n, m, theta, phi));
                        
            }
        }
    }
}

template<typename Vector, typename Source> 
MultipoleExpansion<Vector, Source, 3>::MultipoleExpansion(const Vector& center,
        std::vector<ME*>& expansions): 
        order(expansions[0]->order), coefficients((order+1)*(order+1)),
        center(center) { 
}

template<typename Vector, typename Source> 
double MultipoleExpansion<Vector, Source, 3>::evaluatePotential(
        const Vector& eval_point) const { 

    Vector rel_pos = (eval_point - center).toSpherical(); 
    Complex pot = 0; 

    double r = rel_pos[0];  // Distance of source to center
    double theta = rel_pos[1];  
    double phi = rel_pos[2];  

    // TODO smarter pow
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int n = 0; n <= order; ++n) { 
        Complex delta_pot = 0; 
        for(int m = -n; m <= n; ++m) {
            delta_pot +=  coefficients[coeff_index++] / pow(r, n+1) 
                * sphericalHarmonicY(n, m, theta, phi); 
        }
        pot += 1./(2*n+1) * delta_pot;
    }

    pot *= 4 * PI; 

    return pot.real(); 
}


template<typename Vector, typename Source> 
Vector MultipoleExpansion<Vector, Source, 3>::evaluateForcefield(
        const Vector& eval_point) const { 

    return Vector{}; 
}

} // namespace fmm
