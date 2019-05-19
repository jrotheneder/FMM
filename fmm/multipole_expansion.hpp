#include <complex> 
#include <type_traits> 
#include <numeric> 
#include <algorithm> 
#include <stdexcept> 

#include "series_expansion.hpp"

#ifndef MULTIPOLE_EXPANSION_H
#define MULTIPOLE_EXPANSION_H

namespace fmm {

template<typename Vector, typename Source, std::size_t d>
struct MultipoleExpansion {

    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source> 
struct MultipoleExpansion<Vector, Source, 2> {

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
    MultipoleExpansion(const Vector& center_vec, std::vector<const ME*>& expansions): 
            coefficients(expansions[0]->coefficients.size()), 
            center({center_vec[0], center_vec[1]}), Q(0) {
        
        for(const ME* me : expansions) {
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
struct MultipoleExpansion<Vector, Source, 3>: SeriesExpansion<Vector, Source, 3> {

    using Super = SeriesExpansion<Vector, Source, 3>;
    using ME = MultipoleExpansion; 

    MultipoleExpansion() {}
    MultipoleExpansion(const Vector& center, int order, 
            std::vector<Source>& sources);
    MultipoleExpansion(const Vector& center, std::vector<const ME*>& expansions);

    MultipoleExpansion& operator+=(const MultipoleExpansion& rhs);

    std::vector<Complex> shift(const Vector& shift) const;

    double evaluatePotential(const Vector& eval_point) const override;
    Vector evaluateForcefield(const Vector& eval_point) const override;

    static double sign_fun1(const int k, const int m);

};

template<typename Vector, typename Source> 
MultipoleExpansion<Vector, Source, 3>::MultipoleExpansion(
        const Vector& center, int order, std::vector<Source>& sources): 
        Super(center, order) {

    for(std::size_t i = 0; i < sources.size(); ++i) {

        const auto [r, theta, phi]  // Glorious C++17
            = (sources[i].position - center).toSpherical().data(); 

        // TODO still correct?
        // compute M_n^m according to [1, (5.15)]. However (!)  we use a diff. 
        // convention (Mathematica, Boost, Wikipedia) than [1] for the harmonics)
        // which means that Y_n^{-m} from [1, (5.15)] is replaced by conjugate(Y_n^m).

        // TODO switch to more readable ME(n,m) if no perf. penality
        // TODO smarter pow
        // TODO exploit symmetry in coefficients? 
        unsigned coeff_index = 0; // index of next coefficient to be computed
        for(int n = 0; n <= order; ++n) {  
            for(int m = -n; m <= n; ++m) {  
                                                
                this->coefficients[coeff_index++] +=  
                    sources[i].sourceStrength() * pow(r, n) 
                         * YLM(n, -m, theta, phi);
            }
        }
    }
}

// TODO consider a constructor from shift_vec + ME& together with operator += 
// instead
template<typename Vector, typename Source> 
MultipoleExpansion<Vector, Source, 3>::MultipoleExpansion(const Vector& center,
        std::vector<const ME*>& expansions): 
        Super(center, expansions[0]->order) { 
        
    for(const ME* me : expansions) {

        assert(me->coefficients.size() == this->coefficients.size()); // TODO remove
    
        Vector shift_vec = me->center - this->center; 
        std::vector<Complex> shifted_coefficients = me->shift(shift_vec);

        std::transform (this->coefficients.begin(), this->coefficients.end(), 
            shifted_coefficients.begin(), this->coefficients.begin(), 
            std::plus<Complex>()
        );
    }
}

template<typename Vector, typename Source>
MultipoleExpansion<Vector, Source, 3>& MultipoleExpansion<Vector, Source, 3>::
        operator+=(const MultipoleExpansion& rhs) {

    Super::operator+=(rhs); 
    return *this;
}

template<typename Vector, typename Source> 
std::vector<Complex> MultipoleExpansion<Vector, Source, 3>::shift(
        const Vector& shift) const {

    std::vector<Complex> shifted_coefficients(this->coefficients.size()); 
    const ME& outgoing = *this;  // outgoing expansion

    const auto [r, theta, phi] = shift.toSpherical().data(); 
      
    // TODO smarter pow
    // TODO exploit symmetry in coefficients? 
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int j = 0; j <= this->order; ++j) {  
        for(int k = -j; k <= j; ++k) {  
                                            
            Complex M_jk = 0; // Multipole coeff. M_j^k of the shifted expansion

            for(int n = 0; n <= j; ++n) {
                for(int m = std::max(-n, n+k-j); m <= std::min(n, j+k-n); ++m) {
                    M_jk += outgoing(j-n, k-m) * sign_fun1(k, m) 
                        * this->A_coeff(n, m) * this->A_coeff(j-n, k-m) 
                        / this->A_coeff(j,k) * std::pow(r, n) 
                        * YLM(n, -m, theta, phi);
                }
            }

            shifted_coefficients[coeff_index++] = M_jk;
        }
    }

    return shifted_coefficients;
}

template<typename Vector, typename Source> 
double MultipoleExpansion<Vector, Source, 3>::evaluatePotential(
        const Vector& eval_point) const { 

    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 
    Complex pot = 0; 

    // TODO smarter pow // TODO switch to more readable ME(n,m) if no perf. penality
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int n = 0; n <= this->order; ++n) { 

        for(int m = -n; m <= n; ++m) {
            pot += this->coefficients[coeff_index++] / pow(r, n+1) 
                //* sphericalHarmonicY(n, m, theta, phi); 
                * YLM(n, m, theta, phi); 
        }
    }

    return pot.real(); 
}


template<typename Vector, typename Source> 
Vector MultipoleExpansion<Vector, Source, 3>::evaluateForcefield(
        const Vector& eval_point) const { 

    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 

    // Components of the gradient of the potential evaluated in spherical
    // coordinates (and w.r.t. the spherical coordinate basis) 
    double force_r = 0;     // \hat r component (radial) of the force
    double force_theta = 0; // \hat theta component (polar) of the force
    double force_phi = 0;   // \hat phi component (azimuthal) of the force

    // TODO smarter pow 
    // TODO switch to more readable ME(n,m) if no perf. penality
    // TODO symmetry among coeff.
    // TODO double comp. of YLM (in and outside of d/dtheta), d/phi? 
    // TODO possible recursion relation for YLM / legendre? 
    unsigned coeff_index = 0; // index of next coefficient to be computed

    for(int n = 0; n <= this->order; ++n) { 

        double r_pow = 1/std::pow(r, n+2); 

        for(int m = -n; m <= n; ++m) {

            //Complex ylm_val = YLM(n, m, theta, phi) * r_pow;  
            const Complex M_nm = this->coefficients[coeff_index++];
            
            force_r -= r_pow * (double) (n + 1) * (M_nm * YLM(n, m, theta, phi)).real(); 
            force_theta += r_pow * (M_nm * YLM_deriv_theta(n, m, theta, phi)).real(); 
            force_phi += r_pow * (M_nm * YLM_deriv_phi(n, m, theta, phi)).real() 
                / std::sin(theta); 
        }
    }

    return -Vector{{force_r, force_theta, force_phi}}.toCartesianBasis(theta, phi); 
}

// Implements the function i^(|k|-|m| - |k-m|)
template<typename Vector, typename Source> 
double MultipoleExpansion<Vector, Source, 3>::sign_fun1(
        const int k, const int m) {

    //using namespace std::complex_literals;

    const int exponent = std::abs(k) - std::abs(m) - std::abs(k-m); 
    switch(std::abs(exponent) % 4) {
        case 0 : return 1; 
        case 2 : return -1;  
//      case 1 : return 1i; // These should never occur
//      case 3 : return -1i; 
    }

    throw std::logic_error("Exponent in sign_fun1() is not expected to "
        " be odd. Got input: k = " + std::to_string(k) + ", m = " 
        + std::to_string(m)+ ", exponent is " + std::to_string(exponent) + "\n"
        ); 

}


} // namespace fmm

# endif
