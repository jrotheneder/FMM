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

// 2-D Implementation
template<typename Vector, typename Source>
struct LocalExpansion<Vector, Source, 2>: SeriesExpansion<Vector, Source, 2> {

    using ME = MultipoleExpansion<Vector, Source, 2>;
    using LE = LocalExpansion;

    Complex center; // Center of the expansion
    std::vector<Complex> coefficients; // Local expansion coefficients

    LocalExpansion() {} // Empty default constructor
    LocalExpansion(Vector center_vec, std::size_t order): 
            center({center_vec[0], center_vec[1]}), coefficients(order+1) {}

    LocalExpansion(const Vector& center_vec, std::vector<const ME*> expansions):
            center({center_vec[0], center_vec[1]}) {
        
        assert(expansions.size() > 0); 

        std::size_t order = expansions[0]->coefficients.size();         
        coefficients.resize(order+1); //Store 0-th coeff expl. in LE


        // Convert the supplied multipole expansions to a local expansion
        // TODO this deserves its own function
        for(const ME* me : expansions) { 
            
            Complex a_0 = me->Q; // 0-th coeff of multipole exp.
            Complex z_0 = me->center - this->center; // ME center rel. to this->center
            
            // TODO precompute pows
            // Compute b_0:
            coefficients[0] += a_0 * log(-z_0); 
            for(std::size_t k = 1; k <= order; ++k) {
                double sign = k % 2 == 0 ? 1 : -1;  
                coefficients[0] += 
                    sign * me->coefficients[k-1] / pow(z_0, k); // [(4.18), 1] 
            }
            // Compute b_l for 1 <= l <= order
            for(std::size_t l = 1; l <= order; ++l) {

                Complex b_l = -a_0/(double)l; 

                for(std::size_t k = 1; k < order; ++k) { // TODO possibly reverse sum
                    double sign = k % 2 == 0 ? 1 : -1;  
                    b_l += sign * me->coefficients[k-1]/pow(z_0,k) 
                        * binomial(l+k-1, k-1); // [(4.19), 1]
                }

                b_l /= pow(z_0, l); // [(4.19), 1]
                coefficients[l] += b_l; 
            }
        }
    }

    LocalExpansion(const Vector& center_vec, const LocalExpansion& incoming): 
            center({center_vec[0], center_vec[1]}) {
        
        assert(incoming.coefficients.size() > 0); 

        Complex shift_vec = incoming.center - this->center; 
        this->coefficients = incoming.shift(shift_vec); 

    }

    LocalExpansion& operator+=(const LocalExpansion& rhs) {

        assert(coefficients.size() == rhs.coefficients.size());
        assert(this->center == rhs.center);

        std::transform (
            coefficients.begin(), coefficients.end(), rhs.coefficients.begin(), 
            coefficients.begin(), std::plus<Complex>()
        );

        return *this;
    }


    // Shift is the vector (complex number) from the new center to the old
    // center.
    std::vector<Complex> shift(const Complex shift) const { // [(4.21), 1], shift === z0

        std::size_t order = coefficients.size();
        std::vector<Complex> shifted_coefficients{coefficients}; 

        for(std::size_t j = 0; j < order-1; ++j) {
            for(std::size_t k = order-j-2; k < order-1; ++k) {
                shifted_coefficients[k] -= shift * shifted_coefficients[k+1];    
            }
        }

        return shifted_coefficients;
    };

    double evaluatePotential(const Vector& eval_point) const {

        Complex z{eval_point[0], eval_point[1]}; // get complex repr.
        Complex z_rel = z - center; 

        Complex result{};
        Complex z_rel_pow = 1;

        for(std::size_t k = 0; k < coefficients.size(); ++k) {
            result += coefficients[k] * z_rel_pow;
            z_rel_pow *= z_rel; 
        }

        // Result contains the eval. of the local expansion of the multipole
        // expansion of the potential φ = q log(z-z0) (i.e. the 2d gravitational 
        // potential) summed over all source locations z0. 
        // Return -1 times this (i.e. -result.real()) for the electrostatic potential
        return -result.real(); 
    } 

    Vector evaluateForcefield(const Vector& eval_point) const { 

        Complex z{eval_point[0], eval_point[1]}; // get complex repr.
        Complex z_rel = z - center; 

        Complex result{};
        Complex z_rel_pow = 1;

        for(std::size_t k = 1; k < coefficients.size(); ++k) {
            result += (double)k * coefficients[k] * z_rel_pow;
            z_rel_pow *= z_rel; 
        }

        // Result contains the eval. of the local expansion of the multipole
        // expansion of the potential φ = q log(z-z0) (i.e. the 2d gravitational 
        // potential) summed over all source locations z0. 
        // Return -1 times this (i.e. -result.real()) for the electrostatic potential
        return {{result.real(), -result.imag()}}; 
    }
};

// 3-D Implementation
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

    double sign_fun2(const int k, const int m) const;
    double sign_fun3(const int k, const int m) const;

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

    const auto [r, theta, phi]  
        = (incoming.center - this->center).toSpherical().data(); 

    const typename Super::YlmTable ylm(2 * this->order, theta, phi); 
    double* r_pow_table = new double[2 * this->order + 1];  

    r_pow_table[0] = r; 
    for(int i = 1; i <= 2 * this->order; ++i) { 
        r_pow_table[i] = r_pow_table[i-1] * r; 
    }

    //unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int j = 0; j <= this->order; ++j) {
        for(int k = -j; k <= j; ++k) {

            Complex L_jk = 0; // Coeff. L_j^k of the created local expansion 

////        double A_jk = Super::A_coeff(j, k); // A_j^k

            for(int n = 0; n <= this->order; ++n) {

                double sign = n % 2 ? -1 : 1; 
////            double A_jn_mk = Super::A_coeff(j+n, -n-k); // stores A_{j+n}^{m-k}
////            double A_nm = Super::A_coeff(n, -n);        // stores A_n^m

                for(int m = -n; m <= n; ++m) {
                    L_jk += incoming(n, m) * sign_fun2(k, m) * sign
                        //* Super::A_coeff(n, m) * A_jk
                        /// (Super::A_coeff(j+n, m-k) 
////                    * A_nm * A_jk / (A_jn_mk
                        * this->alm_table(n,m) * this->alm_table(j,k) 
                        / (this->alm_table(j+n, m-k) 
                        //* pow(r, j+n+1)) 
                        * r_pow_table[j+n])  
                        //* YLM(j+n, m-k, theta, phi); 
                        * ylm(j+n, m-k); 

////                A_nm = Super::A_coeff_next(n, m, A_nm);  
////                A_jn_mk = Super::A_coeff_next(j+n, m-k, A_jn_mk);  
                }
            }

            //this->coefficients[coeff_index++] = L_jk; //TODO clear
            (*this)(j,k) = L_jk;
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
    const LE& outgoing = *this;  // outgoing expansion

    const auto [r, theta, phi] = shift.toSpherical().data(); 
      
    // TODO smarter pow
    // TODO exploit symmetry in coefficients? 
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int j = 0; j <= this->order; ++j) {  
        for(int k = -j; k <= j; ++k) {  
                                            
            Complex L_jk = 0; // Multipole coeff. L_j^k of the shifted expansion

            for(int n = j; n <= this->order; ++n) {
                double sign = (n + j) % 2 ? -1 : 1; 
                for(int m = std::max(-n, k+j-n); m <= std::min(n, k+n-j); ++m) {
                    L_jk += outgoing(n, m) * sign_fun3(k, m) * sign
                        * this->A_coeff(n-j, m-k) * this->A_coeff(j, k) 
                        / this->A_coeff(n,m) * std::pow(r, n - j) 
                        * YLM(n - j, m - k, theta, phi);
                }
            }

            shifted_coefficients[coeff_index++] = L_jk;
        }
    }

    return shifted_coefficients;
}

template<typename Vector, typename Source>
double LocalExpansion<Vector, Source, 3>::evaluatePotential(
        const Vector& eval_point) const {
 
    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 
    Complex pot = 0; 

    // TODO smarter pow // TODO switch to more readable ME(n,m) if no perf. penality
    unsigned coeff_index = 0; // index of next coefficient to be computed
    for(int n = 0; n <= this->order; ++n) { 

        for(int m = -n; m <= n; ++m) {
//          std::cout << YLM(n, m, theta, phi) << "\n";
//          std::cout << this->coefficients[coeff_index] * pow(r, n) 
//              * YLM(n, m, theta, phi) << "\n";
            pot +=  this->coefficients[coeff_index++] * pow(r, n) 
                * YLM(n, m, theta, phi); 
        }
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

// Implements the function i^(|k-m|-|k|-|m|)
template<typename Vector, typename Source> 
double LocalExpansion<Vector, Source, 3>::sign_fun2(
        const int k, const int m) const {

    //using namespace std::complex_literals;

    const int exponent = std::abs(k-m) - std::abs(k) - std::abs(m); 
    switch(exponent % 4) {
        case 0  : return 1; 
        case 2  : return -1;  
        case -2 : return -1; 
    }

    throw std::logic_error("Exponent in sign_fun2() is not expected to "
        " be odd. Got input: k = " + std::to_string(k) + ", m = " 
        + std::to_string(m)+ ", exponent is " + std::to_string(exponent) + "\n"
        ); 

}

// Implements the function i^(|m|-|m-k|-|k|)
template<typename Vector, typename Source> 
double LocalExpansion<Vector, Source, 3>::sign_fun3(
        const int k, const int m) const {

    //using namespace std::complex_literals;

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
