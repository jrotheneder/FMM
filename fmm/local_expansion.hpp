#include "series_expansion.hpp"
#include "multipole_expansion.hpp"

#ifndef LOCAL_EXPANSION_H
#define LOCAL_EXPANSION_H

namespace fmm {

template<std::size_t d>
struct LocalExpansion: SeriesExpansion<d> {};


/******************************************************************************/
/*                      2D Local Expansion Implementation                     */
/******************************************************************************/
template<> struct LocalExpansion<2>: SeriesExpansion<2> {

    using Vector = Vector_<2>; 
    using PointSource = PointSource_<2>;
    using Super = SeriesExpansion<2>;

    using ME = MultipoleExpansion<2>; 
    using LE = LocalExpansion;

    LocalExpansion(): Super() {} // Empty default constructor
    LocalExpansion(const Vector& center, std::size_t order): Super(center, order) {}
    LocalExpansion(const Vector& center, std::size_t order, 
            std::vector<PointSource>& sources);
    LocalExpansion(const Vector& center, std::vector<const ME*> expansions);
    LocalExpansion(const Vector& center, const ME& incoming);
    LocalExpansion(const Vector& center, const LE& incoming);

    std::vector<Complex> multipoleToLocal(const ME& incoming) const;
    std::vector<Complex> shift(const Complex shift) const;

    double evaluatePotential(const Vector& eval_point) const;
    Vector evaluateForcefield(const Vector& eval_point) const;

};

LocalExpansion<2>::LocalExpansion(const Vector_<2>& center, 
        const ME& incoming): Super(center, incoming.order) {

    assert(incoming.order > 0); 

    this->coefficients = multipoleToLocal(incoming); 
}

LocalExpansion<2>::LocalExpansion(const Vector_<2>& center, std::size_t order, 
        std::vector<PointSource_<2>>& sources): Super(center, order) {
    
    // Compute series expansion coefficients a_0 through a_order for every
    // source 
    for(std::size_t i = 0; i < sources.size(); ++i) {

        PointSource& src = sources[i];

        Complex z{src.position[0], src.position[1]};
        Complex z_rel = z - this->center; // Express z in box-local coordinates
        Complex z_rel_inv_pow = 1./z_rel; 

        this->coefficients[0] += src.sourceStrength() * std::log(z_rel); 
        for(std::size_t j = 1; j <= order; ++j) {
            this->coefficients[j] -=  
                src.sourceStrength() * z_rel_inv_pow / (double)j; 
            z_rel_inv_pow /= z_rel; 
        }
    }
}

LocalExpansion<2>::LocalExpansion(const Vector_<2>& center, 
        std::vector<const ME*> expansions): Super(center, expansions.at(0)->order) {
    
    for(const ME* me : expansions) { 
        *this += LocalExpansion(center, *me); 
    }
}

LocalExpansion<2>::LocalExpansion(const Vector_<2>& center, 
        const LocalExpansion& incoming): Super(center, incoming.order)  {
    
    assert(incoming.order > 0); 

    Complex shift_vec = incoming.center - this->center; 
    this->coefficients = incoming.shift(shift_vec); 

}

std::vector<Complex> LocalExpansion<2>::multipoleToLocal(
        const ME& incoming) const {

    Complex z0 = incoming.center - this->center; // ME center rel. to this->center
    std::vector<Complex> coefficients(this->order + 1); 

    const tables::BinomialTable& binomial_table = Super::binomial_table; 
    const tables::PowTable<Complex> z0_inv_pow_table(1./z0, this->order);   
    

    coefficients[0] = incoming(0) * std::log(-z0); 
    for(int k = 1; k <= this->order; ++k) {
        double sign = k % 2 == 0 ? 1 : -1;  
        coefficients[0] += sign * incoming(k) * z0_inv_pow_table(k); // [(4.18), 1] 
    }

    // Compute b_l for 1 <= l <= order
    for(int l = 1; l <= this->order; ++l) {

        Complex b_l = -incoming(0)/(double)l; 

        for(int k = 1; k < this->order; ++k) { 
            double sign = k % 2 == 0 ? 1 : -1;  
            b_l += sign * incoming(k) * z0_inv_pow_table(k)   
                * binomial_table(l+k-1, k-1); // [(4.19), 1]
        }

        b_l *= z0_inv_pow_table(l); // [(4.19), 1]
        coefficients[l] = b_l; 
    }

    return coefficients;
}

// Shift is the vector (complex number) from the new center to the old
// center. In [(4.21), 1], shift === z0.
std::vector<Complex> LocalExpansion<2>::shift(const Complex shift) const { 
    
    std::vector<Complex> shifted_coefficients{this->coefficients}; 

    for(int j = 0; j < this->order; ++j) {
        for(int k = this->order-j-1; k < this->order; ++k) {
            shifted_coefficients[k] -= shift * shifted_coefficients[k+1];    
        }
    }

    return shifted_coefficients;
}

double LocalExpansion<2>::evaluatePotential(const Vector_<2>& eval_point) const {

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    const Complex z_rel = z - this->center; 

    Complex result{};
    Complex z_rel_pow = 1;

    for(std::size_t k = 0; k < this->coefficients.size(); ++k) {
        result += (*this)(k) * z_rel_pow;
        z_rel_pow *= z_rel; 
    }

    // Return -result.real() for the electrostatic potential, 
    // +result.real() for the gravitational potential. 
    return result.real(); 
} 


Vector_<2> LocalExpansion<2>::evaluateForcefield(
        const Vector_<2>& eval_point) const { 

    Complex z{eval_point[0], eval_point[1]}; // get complex repr.
    Complex z_rel = z - this->center; 

    Complex result{};
    Complex z_rel_pow = 1;

    for(std::size_t k = 1; k < this->coefficients.size(); ++k) {
        result += (double)k * (*this)(k) * z_rel_pow;
        z_rel_pow *= z_rel; 
    }
    // The gravitational field is given by {{-result.real(), result.imag()}}.
    // For the electric field, return {{result.real(), -result.imag()}}.
    return {{-result.real(), result.imag()}}; 
}

/******************************************************************************/
/*                      3D Local Expansion Implementation                    */
/******************************************************************************/
template<> struct LocalExpansion<3>: SeriesExpansion<3> {

    using Vector = Vector_<3>; 
    using PointSource = PointSource_<3>;
    using Super = SeriesExpansion<3>;
    using ME = MultipoleExpansion<3>;
    using LE = LocalExpansion<3>;

    LocalExpansion() {} // Empty default constructor
    LocalExpansion(const Vector& center, std::size_t order);
    LocalExpansion(const Vector& center, int order, std::vector<PointSource>& sources);
    LocalExpansion(const Vector& center, std::vector<const ME*> expansions);
    LocalExpansion(const Vector& center, const ME& incoming);
    LocalExpansion(const Vector& center, const LE& incoming);

    LocalExpansion& operator+=(const LocalExpansion& rhs);

    std::vector<Complex> shift(const Vector& shift) const;

    double evaluatePotential(const Vector& eval_point) const override;
    Vector evaluateForcefield(const Vector& eval_point) const override;

    static double le_sign_fun2(const int k, const int m);
};


LocalExpansion<3>::LocalExpansion(const Vector_<3>& center, 
        std::size_t order): Super(center, order) {}

LocalExpansion<3>::LocalExpansion(const Vector_<3>& center, int order, 
        std::vector<PointSource_<3>>& sources): Super(center, order) {

    LE& self = *this; 

    for(std::size_t i = 0; i < sources.size(); ++i) {

        const auto [r, theta, phi]  // Glorious C++17
            = (sources[i].position - center).toSpherical().data(); 

        // Precomputed values of Y_l^m(theta, phi) 
        const typename tables::YlmTable sphericalHarmonicY(this->order, theta, phi); 

        double r_pow_inv = 1/r; 
        for(int n = 0; n <= order; ++n) {  
            for(int m = -n; m <= n; ++m) {  
                                                
                self(n,m) += sources[i].sourceStrength() * r_pow_inv
                        * sphericalHarmonicY(n, -m);
            }
            r_pow_inv /= r; 
        }
    }
}

LocalExpansion<3>::LocalExpansion(const Vector_<3>& center, 
        std::vector<const ME*> incoming): Super(center, incoming.at(0)->order) {

    // Convert the supplied multipole expansions to a local expansion
    for(const ME* me : incoming) {
        *this += LocalExpansion(center, *me);
    }
}

LocalExpansion<3>::LocalExpansion(const Vector_<3>& center, const ME& incoming): 
        Super(center, incoming.order) {

    LE& self = *this; 
    const auto [r, theta, phi]  
        = (incoming.center - this->center).toSpherical().data(); 

    // Precomputed values of Y_l^m(theta, phi), A_l^m, sign patterns & powers of  r
    const tables::YlmTable sphericalHarmonicY(2 * this->order, theta, phi); 
    tables::AlmTable& A = Super::alm_table;
    tables::SignTable& sign2 = Super::sign_fun2_table;
    tables::PowTable<double> 
        r_inv_pow_table(1/r, 2*this->order+1);

    for(int j = 0; j <= this->order; ++j) {
        for(int k = -j; k <= j; ++k) {

            Complex L_jk = 0; // Coeff. L_j^k of the created local expansion 

            for(int n = 0; n <= this->order; ++n) {
                Complex accumulant = 0; 
                for(int m = -n; m <= n; ++m) {
                    const double A_coeff = sign2(k, m) * A(n,m) / A(j+n, m-k);
                    accumulant += A_coeff
                        * (incoming(n, m) * sphericalHarmonicY(j+n, m-k)); 
                }

                double sign = n % 2 ? -1 : 1; 
                L_jk += sign * r_inv_pow_table(j+n+1) * accumulant ;
            }

            L_jk *= A(j,k);
            self(j,k) = L_jk;
        }
    }
}

LocalExpansion<3>::LocalExpansion(const Vector_<3>& center, const LE& incoming): 
        Super(center, incoming.order) {

    Vector shift_vec = incoming.center - this->center; 
    this->coefficients = incoming.shift(shift_vec);
}

LocalExpansion<3>& LocalExpansion<3>::operator+=(const LocalExpansion& rhs) {

    Super::operator+=(rhs); 
    return *this;
}

// Shift is the vector from the new center to the old center.
std::vector<Complex> LocalExpansion<3>::shift(const Vector_<3>& shift) const {

    std::vector<Complex> shifted_coefficients(this->coefficients.size()); 
    const LE& outgoing = *this; // outgoing expansion

    const auto [r, theta, phi] = shift.toSpherical().data(); 

    // Precomputed values of Y_l^m(theta, phi), A_l^m, sign patterns & powers of r
    const typename tables::YlmTable sphericalHarmonicY(2 * this->order, theta, phi); 
    typename tables::AlmTable& A = Super::alm_table;
    typename tables::SignTable& sign3 = Super::sign_fun3_table;
    typename tables::PowTable<double> r_pow_table(r, this->order);

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

                L_jk += sign * r_pow_table(n - j) * accumulant;

            }

            shifted_coefficients[coeff_index++] = A(j, k) * L_jk;
        }
    }

    return shifted_coefficients;
}

double LocalExpansion<3>::evaluatePotential(const Vector_<3>& eval_point) 
        const {
 
    const LE& self = *this; 
    const auto [r, theta, phi] = (eval_point - self.center).toSpherical().data(); 
    const typename tables::YlmTable sphericalHarmonicY(self.order, theta, phi); 

    Complex pot = 0; 

    double r_pow = 1; 
    for(int n = 0; n <= self.order; ++n) { 
        for(int m = -n; m <= n; ++m) {
            pot += self(n,m) * r_pow * sphericalHarmonicY(n, m); 
        }
        r_pow *= r; 
    }

    // +pot.real() for the gravitational potential, 
    // -pot.real() for the electrostatic potential
    // n.b. this distinction is handled within the FmmTree classes
    return pot.real(); 
}

Vector_<3> LocalExpansion<3>::evaluateForcefield(const Vector_<3>& eval_point) 
        const { 

    const auto [r, theta, phi] = (eval_point - this->center).toSpherical().data(); 

    // Precomputed values of Y_l^m(theta, phi) and its derivativesj
    const typename tables::YlmDerivTable ylm_table(this->order, theta, phi); 

    // Components of the gradient of the potential evaluated in spherical
    // coordinates (and w.r.t. the spherical coordinate basis) 
    double force_r = 0;     // \hat r component (radial) of the force
    double force_theta = 0; // \hat theta component (polar) of the force
    double force_phi = 0;   // \hat phi component (azimuthal) of the force

    double r_pow = 1; 
    for(int n = 1; n <= this->order; ++n) { 

        for(int m = -n; m <= n; ++m) {

            const Complex L_nm = (*this)(n,m) ;

            force_r += r_pow * (double) n * (L_nm * ylm_table.Y(n,m)).real(); 
            force_theta += r_pow * (L_nm * ylm_table.dYdtheta(n, m)).real(); 
            force_phi += r_pow * (L_nm * ylm_table.dYdphi(n, m)).real() 
                / std::sin(theta); 
        }

        r_pow *= r; 
    }

    // +Vector{{force_r, force_theta, force_phi}} for the gravitational field, 
    // -Vector{{force_r, force_theta, force_phi}} for the electrostatic field
    // n.b. this distinction is handled within the FmmTree classes
    return Vector{{force_r, force_theta, force_phi}}.toCartesianBasis(theta, phi); 
}

} // namespace fmm

#endif
