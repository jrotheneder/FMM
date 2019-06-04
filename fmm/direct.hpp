#include <algorithm> 
#include <functional> 
#include <stdexcept> 

#ifndef DIRECT_H
#define DIRECT_H

namespace fmm {
namespace fields {

template<typename Vector, typename Source, std::size_t d>
double electrostaticPotential(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {

    double r = (src.position - evaluation_point).norm(); 

    if constexpr(d == 2) { return -src.q * log(r + eps); }
    else { return src.q / (r + eps); }
}

template<typename Vector, typename Source, std::size_t d>
double gravitationalPotential(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {
    return -electrostaticPotential(src, evaluation_point, eps); 
}

template<typename Vector, typename Source, std::size_t d>
Vector electrostaticForce(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {

    Vector diff = evaluation_point - src.position; 
    double r = diff.norm(); 

    if constexpr(d == 2) { return src.q / (r*r + eps) * diff; }
    else { return src.q / (r*r*r + eps) * diff; }
}

template<typename Vector, typename Source, std::size_t d>
Vector gravitationalForce(const Source& src, const Vector& evaluation_point,
        const double eps = 0) {
    return - electrostaticForce(src, evaluation_point, eps); 
}

// safe implementations: these check whether source.position == evaluation_point
// and return zero in that case
template<typename Vector, typename Source, std::size_t d>
double safeElectrostaticPotential(const Source& src, const Vector& evaluation_point,
        const double eps = 0) {

    double r = (src.position - evaluation_point).norm(); 

    if(r == 0) {
        return 0; 
    /* throw std::runtime_error("Cannot evaluate field at source location."); */
    }
    
    return electrostaticPotential<Vector, Source, d>(src, evaluation_point, eps); 
}

template<typename Vector, typename Source, std::size_t d>
double safeGravitationalPotential(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {
    return -electrostaticPotential_safe(src, evaluation_point, eps); 
}

template<typename Vector, typename Source, std::size_t d>
Vector safeElectrostaticForce(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {

    Vector diff = evaluation_point - src.position; 
    double r = diff.norm(); 

    if(r == 0) {
        return 0; 
    /* throw std::runtime_error("Cannot evaluate field at source location."); */
    }
    
    return electrostaticForce<Vector, Source, d>(src, evaluation_point, eps);
}

template<typename Vector, typename Source, std::size_t d>
Vector safeGravitationalForce(const Source& src, const Vector& evaluation_point, 
        const double eps = 0) {
    return - electrostaticForce(src, evaluation_point, eps); 
}

} //namespace fields

// InteractionResult can be e.g. double, Vector for potentials and forces
template<typename Vector, typename Source, typename InteractionResult> 
InteractionResult evaluateInteraction(const std::vector<Source>& sources, 
        const Vector& evaluation_point, const double eps, 
        std::function <InteractionResult (const Source&, const Vector&, 
        const double)> interactionFunction) {

    InteractionResult res{}; // Default init must be '0' element of result type
    res = std::accumulate(sources.begin(), sources.end(), res, 
        [&interactionFunction, &evaluation_point, eps](InteractionResult& accumulant, 
        const Source& src) -> InteractionResult { 
            return accumulant + interactionFunction(src, evaluation_point, eps); 
        }
    );

    return res; 
}
} //namespace fmm

#endif
