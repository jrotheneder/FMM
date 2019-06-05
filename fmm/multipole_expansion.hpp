#include <complex> 
#include <type_traits> 
#include <numeric> 
#include <algorithm> 
#include <stdexcept> 

#include "series_expansion.hpp"

#ifndef MULTIPOLE_EXPANSION_H
#define MULTIPOLE_EXPANSION_H

namespace fmm {

template<unsigned d> struct MultipoleExpansion: SeriesExpansion<d> {};


/******************************************************************************/
/*                      2D Multipole Expansion Implementation                    */
/******************************************************************************/
template<> struct MultipoleExpansion<2>: SeriesExpansion<2> {

    using Vector = Vector_<2>; 
    using PointSource = PointSource_<2>;
    using Super = SeriesExpansion<2>;
    using ME = MultipoleExpansion<2>;

    MultipoleExpansion(): Super() {}; 
    MultipoleExpansion(const Vector& center, unsigned order, 
            std::vector<PointSource>& sources);
    MultipoleExpansion(const Vector& center, std::vector<const ME*>& expansions);

    std::vector<Complex> shift(const Complex& shift) const;

    double evaluatePotential(const Vector& eval_point) const;
    Vector evaluateForcefield(const Vector& eval_point) const;

};

MultipoleExpansion<2>::MultipoleExpansion(const Vector_<2>& center, 
        unsigned order, std::vector<PointSource_<2>>& sources): Super(center, order) {
    
    // Compute series expansion coefficients a_0 through a_order: 
    for(std::size_t i = 0; i < sources.size(); ++i) {

        PointSource& src = sources[i];

        Complex z{src.position[0], src.position[1]};
        Complex z_rel = z - this->center; // Express z in box-local coordinates
        Complex z_rel_pow = z_rel; 

        (*this)(0) += src.sourceStrength(); 
        for(unsigned j = 1; j <= order; ++j) {
            (*this)(j) -=  src.sourceStrength() * z_rel_pow / (double)j; 
            z_rel_pow *= z_rel; 
        }
    }
}

// Construct expansion from other expansions by shifting
MultipoleExpansion<2>::MultipoleExpansion(const Vector& center, 
        std::vector<const ME*>& expansions): Super(center, expansions.at(0)->order) {
    
    for(const ME* me : expansions) {

        Complex shift_vec = me->center - this->center; 
        std::vector<Complex> shifted_coefficients 
            = me->shift(shift_vec);

        std::transform (this->coefficients.begin(), this->coefficients.end(), 
            shifted_coefficients.begin(), this->coefficients.begin(), 
            std::plus<Complex>()
        );
    }
}

// Shift is the vector (complex number) from the new to the old center.
std::vector<Complex>  MultipoleExpansion<2>::shift(const Complex& shift) const {

    std::vector<Complex> shifted_coefficients(this->order + 1); 

    const typename tables::BinomialTable& binomial_table = Super::binomial_table; 
    typename tables:: template PowTable<Complex> shift_pow_table(shift, this->order);

    double Q = (*this)(0).real(); 
    shifted_coefficients[0] = Q; 
    for(int l = 1; l <= this->order; ++l) {
        shifted_coefficients[l] = -Q * shift_pow_table(l)  / (double)l; 
        for(int k = 1; k <= l; ++k) {
            shifted_coefficients[l] += (*this)(k) * shift_pow_table(l-k) 
                * binomial_table(l-1, k-1);  // [(4.15), 1]  
        }
    }

    return shifted_coefficients; 
}

double MultipoleExpansion<2>::evaluatePotential(const Vector_<2>& eval_point) const {

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    Complex z_rel = z - this->center; 

    double Q = (*this)(0).real(); 
    Complex result = Q * log(z_rel); 
    
    Complex z_rel_inv_pow = 1./z_rel; 
    for(unsigned j = 1; j < this->coefficients.size(); ++j) { 
        result += (*this)(j) * z_rel_inv_pow; 
        z_rel_inv_pow /= z_rel; 
    }

    // Return -result.real() for the electrostatic potential, 
    // +result.real() for the gravitational potential. 
    return result.real(); 
} 

Vector_<2> MultipoleExpansion<2>::evaluateForcefield(const Vector_<2>& eval_point) 
        const { 

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    Complex z_rel = z - this->center; 
    double Q = (*this)(0).real(); 

    Complex result = Q / z_rel;

    Complex z_rel_inv_pow = 1./(z_rel * z_rel); 
    for(unsigned j = 1; j < this->coefficients.size(); ++j) { 
        result -= (double)j * (*this)(j) * z_rel_inv_pow;
        z_rel_inv_pow /= z_rel; 
    }

    // The gravitational field is given by {{-result.real(), result.imag()}}.
    // For the electric field, return {{result.real(), -result.imag()}}.
    return {{-result.real(), result.imag()}}; 
}

/******************************************************************************/
/*                      3D Multipole Expansion Implementation                 */
/******************************************************************************/
template<> struct MultipoleExpansion<3>: SeriesExpansion<3> {

    using Vector = Vector_<3>; 
    using PointSource = PointSource_<3>;
    using Super = SeriesExpansion<3>;
    using ME = MultipoleExpansion;  

    MultipoleExpansion() {}
    MultipoleExpansion(const Vector& center, int order, 
            std::vector<PointSource>& sources);
    MultipoleExpansion(const Vector& center, std::vector<const ME*>& expansions);

    MultipoleExpansion& operator+=(const MultipoleExpansion& rhs);

    std::vector<Complex> shift(const Vector& shift) const;

    double evaluatePotential(const Vector& eval_point) const override;
    Vector evaluateForcefield(const Vector& eval_point) const override;

};

MultipoleExpansion<3>::MultipoleExpansion(const Vector_<3>& center, int order, 
        std::vector<PointSource_<3>>& sources): Super(center, order) {

    ME& self = *this; 

    for(std::size_t i = 0; i < sources.size(); ++i) {

        const auto [r, theta, phi]  // Glorious C++17
            = (sources[i].position - center).toSpherical().data(); 

        // Precomputed values of Y_l^m(theta, phi) 
        const typename tables::YlmTable sphericalHarmonicY(this->order, theta, phi); 

        double r_pow = 1; 
        for(int n = 0; n <= order; ++n) {  
            for(int m = -n; m <= n; ++m) {  
                                                
                self(n,m) += sources[i].sourceStrength() * r_pow
                        * sphericalHarmonicY(n, -m);
            }
            r_pow *= r; 
        }
    }
}

MultipoleExpansion<3>::MultipoleExpansion(const Vector_<3>& center,
        std::vector<const ME*>& expansions): Super(center, expansions[0]->order) { 
        
    for(const ME* me : expansions) {

        Vector shift_vec = me->center - this->center; 
        std::vector<Complex> shifted_coefficients = me->shift(shift_vec);

        std::transform (this->coefficients.begin(), this->coefficients.end(), 
            shifted_coefficients.begin(), this->coefficients.begin(), 
            std::plus<Complex>()
        );
    }
}

MultipoleExpansion<3>& MultipoleExpansion<3>::operator+=(const MultipoleExpansion& rhs) {

    Super::operator+=(rhs); 
    return *this;
}

std::vector<Complex> MultipoleExpansion<3>::shift(const Vector& shift) const {

    std::vector<Complex> shifted_coefficients(this->coefficients.size()); 
    const ME& outgoing = *this;  // outgoing expansion

    const auto [r, theta, phi] = shift.toSpherical().data(); 

    // Precomputed values of Y_l^m(theta, phi), A_l^m, sign patterns & powers of r
    const typename tables::YlmTable sphericalHarmonicY(this->order, theta, phi); 
    typename tables::AlmTable& A = Super::alm_table;
    typename tables::SignTable& sign1 = Super::sign_fun1_table;
    typename tables::PowTable<double> r_pow_table(r, this->order);
      
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int j = 0; j <= this->order; ++j) {  
        for(int k = -j; k <= j; ++k) {  
                                            
            Complex M_jk = 0; // Multipole coeff. M_j^k of the shifted expansion

            for(int n = 0; n <= j; ++n) {

                Complex accumulant = 0; 

                for(int m = std::max(-n, n+k-j); m <= std::min(n, j+k-n); ++m) {
                    accumulant += sign1(k, m) * A(n, m) * A(j-n, k-m)     
                        * (outgoing(j-n, k-m) * sphericalHarmonicY(n, -m));
                }

                M_jk += accumulant * r_pow_table(n);
            }

            shifted_coefficients[coeff_index++] = M_jk / A(j,k);
        }
    }

    return shifted_coefficients;
}

double MultipoleExpansion<3>::evaluatePotential(const Vector_<3>& eval_point) const { 

    const ME& self = *this; 

    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 
    Complex pot = 0; 

    const typename tables::YlmTable sphericalHarmonicY(this->order, theta, phi); 

    double r_pow = 1/r; 
    for(int n = 0; n <= this->order; ++n) { 
        for(int m = -n; m <= n; ++m) {
            pot += self(n,m) * r_pow * sphericalHarmonicY(n, m); 
        }
        r_pow /= r; 
    }

    // +pot.real() for the gravitational potential, 
    // -pot.real() for the electrostatic potential
    // n.b. this distinction is handled within the FmmTree classes
    return pot.real(); 
}

Vector_<3> MultipoleExpansion<3>::evaluateForcefield(const Vector_<3>& eval_point) 
        const { 

    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 

    // Precomputed values of Y_l^m(theta, phi) and its derivativesj
    const typename tables::YlmDerivTable ylm_table(this->order, theta, phi); 

    // Components of the gradient of the potential evaluated in spherical
    // coordinates (and w.r.t. the spherical coordinate basis) 
    double force_r = 0;     // \hat r component (radial) of the force
    double force_theta = 0; // \hat theta component (polar) of the force
    double force_phi = 0;   // \hat phi component (azimuthal) of the force

    double r_pow = 1/(r*r); 

    for(int n = 0; n <= this->order; ++n) {  
        for(int m = -n; m <= n; ++m) {

            const Complex& M_nm = (*this)(n,m) ;
            
            force_r -= r_pow * (double) (n + 1) * (M_nm * ylm_table.Y(n, m)).real(); 
            force_theta += r_pow * (M_nm * ylm_table.dYdtheta(n, m)).real(); 
            force_phi += r_pow * (M_nm * ylm_table.dYdphi(n, m)).real() 
                / std::sin(theta); 
        }

        r_pow /= r; 
    }

    // +Vector{{force_r, force_theta, force_phi}} for the gravitational field, 
    // -Vector{{force_r, force_theta, force_phi}} for the electrostatic field
    // n.b. this distinction is handled within the FmmTree classes
    return Vector{{force_r, force_theta, force_phi}}.toCartesianBasis(theta, phi); 
}

} // namespace fmm

# endif
