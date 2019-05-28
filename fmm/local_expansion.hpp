#include "series_expansion.hpp"
#include "multipole_expansion.hpp"

#ifndef LOCAL_EXPANSION_H
#define LOCAL_EXPANSION_H

namespace fmm {

typedef std::complex<double> Complex; 

template<typename Vector, typename Source, std::size_t d>
struct LocalExpansion: SeriesExpansion<Vector, Source, d> {
    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};


/******************************************************************************/
/*                      2D Local Expansion Implementation                     */
/******************************************************************************/
template<typename Vector, typename Source>
struct LocalExpansion<Vector, Source, 2>: SeriesExpansion<Vector, Source, 2> {

    using Super = SeriesExpansion<Vector, Source, 2>;
    using ME = MultipoleExpansion<Vector, Source, 2>;
    using LE = LocalExpansion;

    LocalExpansion(): Super() {} // Empty default constructor
    LocalExpansion(const Vector& center, std::size_t order): Super(center, order) {}
    LocalExpansion(const Vector& center, std::vector<const ME*> expansions);
    LocalExpansion(const Vector& center, const ME& incoming);
    LocalExpansion(const Vector& center, const LE& incoming);

    std::vector<Complex> multipoleToLocal(const ME& incoming) const;
    std::vector<Complex> shift(const Complex shift) const;

    double evaluatePotential(const Vector& eval_point) const;
    Vector evaluateForcefield(const Vector& eval_point) const;

};

template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 2>::LocalExpansion(const Vector& center, 
        const ME& incoming): Super(center, incoming.order) {

    assert(incoming.coefficients.size() > 0); 

    this->coefficients = multipoleToLocal(incoming); 
}


template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 2>::LocalExpansion(const Vector& center, 
        std::vector<const ME*> expansions): Super(center, expansions.at(0)->order) {
    
    for(const ME* me : expansions) { 
        *this += LocalExpansion(center, *me); 
    }
}

template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 2>::LocalExpansion(const Vector& center, 
        const LocalExpansion& incoming): Super(center, incoming.order)  {
    
    assert(incoming.coefficients.size() > 0); 

    Complex shift_vec = incoming.center - this->center; 
    this->coefficients = incoming.shift(shift_vec); 

}

template<typename Vector, typename Source>
std::vector<Complex> LocalExpansion<Vector, Source, 2>::multipoleToLocal(
        const ME& incoming) const {

    Complex a_0 = incoming.coefficients[0];  // 0-th coeff of multipole exp.
    Complex z_0 = incoming.center - this->center; // ME center rel. to this->center
    
    // TODO precompute pows
    // Compute b_0:
    std::vector<Complex> coefficients(this->order + 1) ; 

    coefficients[0] += a_0 * log(-z_0); 
    for(int k = 1; k <= this->order; ++k) {
        double sign = k % 2 == 0 ? 1 : -1;  
        coefficients[0] += 
            sign * incoming.coefficients[k] / pow(z_0, k); // [(4.18), 1] 
    }
    // Compute b_l for 1 <= l <= order
    for(int l = 1; l <= this->order; ++l) {

        Complex b_l = -a_0/(double)l; 

        for(int k = 1; k < this->order; ++k) { // TODO reverse sum?
            double sign = k % 2 == 0 ? 1 : -1;  
            b_l += sign * incoming.coefficients[k]/pow(z_0,k) 
                * binomial(l+k-1, k-1); // [(4.19), 1]
        }

        b_l /= pow(z_0, l); // [(4.19), 1]
        coefficients[l] += b_l; 
    }

    return coefficients;
}


// Shift is the vector (complex number) from the new center to the old
// center.
template<typename Vector, typename Source>
std::vector<Complex> LocalExpansion<Vector, Source, 2>::shift(
        const Complex shift) const { // [(4.21), 1], shift === z0

    std::vector<Complex> shifted_coefficients{this->coefficients}; 

    for(int j = 0; j < this->order; ++j) {
        for(int k = this->order-j-1; k < this->order; ++k) {
            shifted_coefficients[k] -= shift * shifted_coefficients[k+1];    
        }
    }

    return shifted_coefficients;
}

template<typename Vector, typename Source>
double LocalExpansion<Vector, Source, 2>::evaluatePotential(
        const Vector& eval_point) const {

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    Complex z_rel = z - this->center; 

    Complex result{};
    Complex z_rel_pow = 1;

    for(std::size_t k = 0; k < this->coefficients.size(); ++k) {
        result += this->coefficients[k] * z_rel_pow;
        z_rel_pow *= z_rel; 
    }

    // Result contains the eval. of the local expansion of the multipole
    // expansion of the potential φ = q log(z-z0) (i.e. the 2d gravitational 
    // potential) summed over all source locations z0. 
    // Return -1 times this (i.e. -result.real()) for the electrostatic potential
    return -result.real(); 
} 


template<typename Vector, typename Source>
Vector LocalExpansion<Vector, Source, 2>::evaluateForcefield(
        const Vector& eval_point) const { 

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    Complex z_rel = z - this->center; 

    Complex result{};
    Complex z_rel_pow = 1;

    for(std::size_t k = 1; k < this->coefficients.size(); ++k) {
        result += (double)k * this->coefficients[k] * z_rel_pow;
        z_rel_pow *= z_rel; 
    }

    // Result contains the eval. of the local expansion of the multipole
    // expansion of the potential φ = q log(z-z0) (i.e. the 2d gravitational 
    // potential) summed over all source locations z0. 
    // Return -1 times this (i.e. -result.real()) for the electrostatic potential
    return {{result.real(), -result.imag()}}; 
}

/******************************************************************************/
/*                      3D Local Expansion Implementation                    */
/******************************************************************************/
template<typename Vector, typename Source>
struct LocalExpansion<Vector, Source, 3>: SeriesExpansion<Vector, Source, 3> {

    using Super = SeriesExpansion<Vector, Source, 3>;
    using ME = MultipoleExpansion<Vector, Source, 3>;
    using LE = LocalExpansion<Vector, Source, 3>;

    LocalExpansion() {} // Empty default constructor
    LocalExpansion(const Vector& center, std::size_t order);
    LocalExpansion(const Vector& center, std::vector<const ME*> expansions);
    LocalExpansion(const Vector& center, const ME& incoming);
    LocalExpansion(const Vector& center, const LE& incoming);

    LocalExpansion& operator+=(const LocalExpansion& rhs);

    std::vector<Complex> shift(const Vector& shift) const;

    double evaluatePotential(const Vector& eval_point) const override;
    Vector evaluateForcefield(const Vector& eval_point) const override;

    static double le_sign_fun2(const int k, const int m);
};


template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 3>::LocalExpansion(const Vector& center, 
        std::size_t order): Super(center, order) {}

// TODO consider deleting this constructor (need to change 2d imp. as well)
template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 3>::LocalExpansion(const Vector& center, 
        std::vector<const ME*> incoming): Super(center, incoming.at(0)->order) {

        // Convert the supplied multipole expansions to a local expansion
        for(const ME* me : incoming) {
            *this += LocalExpansion(center, *me);
        }
}

template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 3>::LocalExpansion(const Vector& center, 
        const ME& incoming): Super(center, incoming.order) {

    LE& self = *this; 
    const auto [r, theta, phi]  
        = (incoming.center - this->center).toSpherical().data(); 

    // Precompute powers of r:
    double* r_pow_table = new double[2 * this->order + 1];  
    r_pow_table[0] = r; 
    for(int i = 1; i <= 2 * this->order; ++i) { 
        r_pow_table[i] = r_pow_table[i-1] * r; 
    }

    // Precomputed values of Y_l^m(theta, phi) & A_l^m 
    const typename Super::YlmTable sphericalHarmonicY(2*this->order, theta, phi); 
    typename Super::template AlmTable& A = Super::alm_table;
    typename Super::template SignTable& sign2 = Super::sign_fun2_table;

    for(int j = 0; j <= this->order; ++j) {
        for(int k = -j; k <= j; ++k) {

            Complex L_jk = 0; // Coeff. L_j^k of the created local expansion 
            double A_jk = A(j,k) ; // A_j^k

            for(int n = 0; n <= this->order; ++n) {

                double sign = n % 2 ? -1 : 1; 
                Complex accumulant = 0; 

                for(int m = -n; m <= n; ++m) {
                    accumulant += sign2(k, m) * A(n,m) 
                        / (A(j+n, m-k) * r_pow_table[j+n])
                        * (incoming(n, m) * sphericalHarmonicY(j+n, m-k)); 
                }

                L_jk += sign * accumulant;
            }

            L_jk *= A_jk;
            self(j,k) = L_jk;
        }
    }

    delete[] r_pow_table;
}

template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 3>::LocalExpansion(const Vector& center, 
        const LE& incoming): Super(center, incoming.order) {

    Vector shift_vec = incoming.center - this->center; 
    this->coefficients = incoming.shift(shift_vec);
}

template<typename Vector, typename Source>
LocalExpansion<Vector, Source, 3>& LocalExpansion<Vector, Source, 3>::
        operator+=(const LocalExpansion& rhs) {

    Super::operator+=(rhs); 

    return *this;
}

// Shift is the vector from the new center to the old center.
template<typename Vector, typename Source>
std::vector<Complex> LocalExpansion<Vector, Source, 3>::shift(
        const Vector& shift) const {

    std::vector<Complex> shifted_coefficients(this->coefficients.size()); 
    const LE& outgoing = *this; // outgoing expansion

    const auto [r, theta, phi] = shift.toSpherical().data(); 

    // Precompute powers of r:
    double* r_pow_table = new double[this->order + 1];  
    r_pow_table[0] = 1; 
    for(int i = 1; i <= this->order; ++i) { 
        r_pow_table[i] = r_pow_table[i-1] * r; 
    }
      
    // Precomputed values of Y_l^m(theta, phi) & A_l^m 
    const typename Super::YlmTable sphericalHarmonicY(2 * this->order, theta, phi); 
    typename Super::template AlmTable& A = Super::alm_table;
    typename Super::template SignTable& sign3 = Super::sign_fun3_table;

    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int j = 0; j <= this->order; ++j) {  
        for(int k = -j; k <= j; ++k) {  
                                            
            Complex L_jk = 0; // Multipole coeff. L_j^k of the shifted expansion

            for(int n = j; n <= this->order; ++n) {

                double sign = (n + j) % 2 ? -1 : 1; 
                Complex accumulant = 0; 

                for(int m = std::max(-n, k+j-n); m <= std::min(n, k+n-j); ++m) {
                    accumulant += sign3(k, m) * A(n-j, m-k) / A(n,m)  
                        * (outgoing(n, m) * sphericalHarmonicY(n - j, m - k));
                }

                L_jk += sign * r_pow_table[n - j] * accumulant;

            }

            shifted_coefficients[coeff_index++] = A(j, k) * L_jk;
        }
    }

    delete[] r_pow_table; 

    return shifted_coefficients;
}

template<typename Vector, typename Source>
double LocalExpansion<Vector, Source, 3>::evaluatePotential(
        const Vector& eval_point) const {
 
    const LE& self = *this; 
    const auto [r, theta, phi] = (eval_point - self.center).toSpherical().data(); 
    const typename Super::YlmTable sphericalHarmonicY(self.order, theta, phi); 

    Complex pot = 0; 

    double r_pow = 1; 
    for(int n = 0; n <= self.order; ++n) { 
        for(int m = -n; m <= n; ++m) {
            pot += self(n,m) * r_pow * sphericalHarmonicY(n, m); 
        }
        r_pow *= r; 
    }

    return pot.real(); 
}

template<typename Vector, typename Source>
Vector LocalExpansion<Vector, Source, 3>::evaluateForcefield(
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
    unsigned coeff_index = 1; // index of next coefficient to be computed
    for(int n = 1; n <= this->order; ++n) { 

        double r_pow = std::pow(r, n-1); 

        for(int m = -n; m <= n; ++m) {

            //Complex ylm_val = YLM(n, m, theta, phi) * r_pow;  
            const Complex L_nm = this->coefficients[coeff_index++];
            
            force_r += r_pow * (double) n * (L_nm * YLM(n, m, theta, phi)).real(); 
            force_theta += r_pow * (L_nm * YLM_deriv_theta(n, m, theta, phi)).real(); 
            force_phi += r_pow * (L_nm * YLM_deriv_phi(n, m, theta, phi)).real() 
                / std::sin(theta); 
        }
    }

    return -Vector{{force_r, force_theta, force_phi}}.toCartesianBasis(theta, phi); 
}


// Implements the function i^(|k-m|-|k|-|m|) for k, m integers
template<typename Vector, typename Source> 
double LocalExpansion<Vector, Source, 3>::le_sign_fun2(
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

} // namespace fmm

#endif
