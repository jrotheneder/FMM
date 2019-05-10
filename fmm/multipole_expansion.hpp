#include <complex> 
#include <type_traits> 
#include <numeric> 
#include <algorithm> 

namespace fmm {

// std::beta requires C++17
double binomial(std::size_t n, std::size_t k) { // TODO consider a lookup table if slow
    return 1 / ((n+1) * std::beta(n-k+1, k+1)); // TODO investigate accuracy
                                                // TODO signal if overflow
}

template<typename Vector, typename Source, std::size_t d, typename = void>
struct MultipoleExpansion {

    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source, std::size_t d> 
struct MultipoleExpansion<Vector, Source, d, typename std::enable_if<d==2>::type> {

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
            center({center_vec[0], center_vec[1]}), Q(0){
        
        for(ME* me : expansions) {
            Q += me->Q;     

            assert(me->coefficients.size() > 0);  // TODO remove in a bit
        
            Complex shift_vec = me->center - this->center; //TODO: check this!!
            std::vector<Complex> shifted_coefficients 
                = me->shift(shift_vec);

            std::transform (coefficients.begin(), coefficients.end(), 
                shifted_coefficients.begin(), coefficients.begin(), 
                std::plus<Complex>()
            );
        }
    }

    // TODO mark const 
    // Shift is the vector (complex number) from the new center to the old
    // center.
    std::vector<Complex> shift(const Complex& shift) {

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

    // TODO mark const
    double evaluatePotential(const Vector& eval_point) {

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

    //TODO mark const
    Vector evaluateForcefield(const Vector& eval_point) { 

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
template<typename Vector, typename Source, std::size_t d> 
struct MultipoleExpansion<Vector, Source, d, typename std::enable_if<d==3>::type> {

    std::vector<double> coefficients; 

    MultipoleExpansion(std::vector<Source> sources = {}) {} //TODO

    std::vector<double> shift(Vector v) { return {}; }
};


} // namespace fmm
