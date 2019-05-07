#include <complex> 
#include <type_traits> 

namespace fmm {

template<typename Vector, typename Source, std::size_t d, typename = void>
struct MultipoleExpansion {
    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source, std::size_t d> 
struct MultipoleExpansion<Vector, Source, d, typename std::enable_if<d==2>::type> {

    std::vector<complex<double>> coefficients; // ME coefficients
    std::complex<double> center;
    double Q_tot; // Total source strength

    MultipoleExpansion() {}
    MultipoleExpansion(std::vector<Source>& sources, const Vector& center_vec,
            std::size_t order): coefficients(order), 
            center({center_vec[0], center_vec[1]}), Q_tot(0) {
        
        // Compute series expansion coefficients a_1 through a_order: 
        for(std::size_t i = 0; i < sources.size(); ++i) {

            std::complex<double> z{sources[i].position[0], sources[i].position[1]};
            std::complex<double> z_rel = z - center;

            for(std::size_t j = 1; j <= order; ++j) {
                coefficients[j-1] -= 
                    sources[i].q * std::pow(z_rel, j) / (double)j; 
            }
            Q_tot += sources[i].q;
        }
    }

    MultipoleExpansion(std::vector<std::vector<Source>*>& source_vectors, 
            const Vector& center_vec, std::size_t order): coefficients(order), 
            center({center_vec[0], center_vec[1]}), Q_tot(0) {
        // Compute series expansion coefficients a_1 through a_order: 
        for(std::size_t k = 0; k < source_vectors.size(); ++k) {

            std::vector<Source>& sources = *source_vectors[k]; 

            for(int i = 0; i < sources.size(); ++i) {

                std::complex<double> z{sources[i].position[0], 
                    sources[i].position[1]};
                std::complex<double> z_rel = z - center;

                for(int j = 1; j <= order; ++j) {
                    coefficients[j-1] -= sources[i].q 
                        * std::pow(z_rel, j) / j; 
                }
                Q_tot += sources[i].q;
            }
        }
    }

    std::vector<double> shift(const Vector& shift_vector) { return {}; }

    double evaluatePotential(const Vector& eval_point) {

        std::complex<double> z{eval_point[0], eval_point[1]}; // get complex repr.
        std::complex<double> z_rel = z - center; 

        std::complex<double> result = Q_tot * log(z_rel); 
        
        // TODO Horner scheme this
        for(int j = coefficients.size(); j >= 1; --j) { // a_j is stored in coeff[j-1]  
            result += coefficients[j-1] / pow(z_rel, j);
        }

        // Result contains the eval. of the power series for φ = q log(z-z0) (i.e.
        // the 2d gravitational potential) summed over all source locations z0. 
        // Return -1 times this for the electrostatic potential
        return -result.real(); 
    } 

    Vector evaluateForcefield(const Vector& eval_point) {

        std::complex<double> z{eval_point[0], eval_point[1]}; // get complex repr.
        std::complex<double> z_rel = z - center; 

        std::complex<double> result = Q_tot / z_rel;

        // TODO Horner scheme this
        for(int j = coefficients.size(); j >= 1; --j) { // a_j is stored in coeff[j-1]  
            result -= (double)j * coefficients[j-1] / pow(z_rel, j+1);
        }

        // Result contains the eval. of *the derivative* φ' of the power series for 
        // φ = q log(z-z0) (i.e. the 2d gravitational potential) summed over all 
        // source locations z0. The force field is then given by 
        // \vec F = (Re φ, -Im φ).
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
