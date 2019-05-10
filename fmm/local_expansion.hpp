namespace fmm {

template<typename Vector, typename Source, std::size_t d, typename = void>
struct LocalExpansion {
    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source, std::size_t d>
struct LocalExpansion<Vector, Source, d, typename std::enable_if<d==2>::type> {

    using Complex = std::complex<double>;
    using ME = MultipoleExpansion<Vector, Source, d>;

    std::vector<Complex> coefficients; // Local expansion coefficients
    Complex center; // Center of the expansion

    LocalExpansion() {}
    LocalExpansion(const Vector& center_vec, std::size_t order, 
            std::vector<ME*> me_expansions): coefficients(order + 1), 
            center({center_vec[0], center_vec[1]}) {


        // Convert the supplied multipiole expansions to a local expansion
        for(ME * me : me_expansions) { 
            
            assert(order == me->coefficients.size());

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

    // Shift is the vector (complex number) from the new center to the old
    // center.
    std::vector<Complex> shift(const Complex shift) { // [(4.21), 1], shift === zo

        std::size_t order = coefficients.size();
        std::vector<Complex> shifted_coefficients{coefficients}; 

        for(std::size_t j = 0; j < order-1; ++j) {
            for(std::size_t k = order-j-2; k < order-1; ++k) {
                shifted_coefficients[k] -= shift * shifted_coefficients[k+1];    
            }
        }

        return shifted_coefficients;
    };

    // TODO mark const
    double evaluatePotential(const Vector& eval_point) {

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

    //TODO mark const
    Vector evaluateForcefield(const Vector& eval_point) { 

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
template<typename Vector, typename Source, std::size_t d >
struct LocalExpansion<Vector, Source, d, typename std::enable_if<d==3>::type> {

    using ME = MultipoleExpansion<Vector, Source, d>;

    std::vector<double> coefficients; 

    LocalExpansion(ME me = {}) {} // TODO 
    std::vector<double> shift(Vector v);
};

} // namespace fmm
