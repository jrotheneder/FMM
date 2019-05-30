#include <algorithm> 
#include <functional> 
#include <stdexcept> 

#ifndef DIRECT_H
#define DIRECT_H

namespace fmm {
namespace fields {

template<typename Vector, typename Source, std::size_t d>
double electrostaticPotential(const Source& src, const Vector& evaluation_point) {

    double r = (src.position - evaluation_point).norm(); 

    if constexpr(d == 2) { return -src.q * log(r); }
    else { return src.q / r; }
}

template<typename Vector, typename Source, std::size_t d>
double gravitationalPotential(const Source& src, const Vector& evaluation_point) {
    return -electrostaticPotential(src, evaluation_point); 
}

template<typename Vector, typename Source, std::size_t d>
Vector electrostaticForce(const Source& src, const Vector& evaluation_point) {

    Vector diff = evaluation_point - src.position; 
    double r = diff.norm(); 

    if constexpr(d == 2) { return src.q / (r*r) * diff; }
    else { return src.q / (r*r*r) * diff; }
}

template<typename Vector, typename Source, std::size_t d>
Vector gravitationalForce(const Source& src, const Vector& evaluation_point) {
    return - electrostaticForce(src, evaluation_point); 
}

// safe implementations: these check whether source.position == evaluation_point
template<typename Vector, typename Source, std::size_t d>
double electrostaticPotential_s(const Source& src, const Vector& evaluation_point) {

    double r = (src.position - evaluation_point).norm(); 

    if(r == 0) {
        return 0; 
    /* throw std::runtime_error("Cannot evaluate field at source location."); */
    }

    if constexpr(d == 2) { return -src.q * log(r); }
    else { return src.q / r; }
}

template<typename Vector, typename Source, std::size_t d>
double gravitationalPotential_s(const Source& src, const Vector& evaluation_point) {
    return -electrostaticPotential_s(src, evaluation_point); 
}

template<typename Vector, typename Source, std::size_t d>
Vector electrostaticForce_s(const Source& src, const Vector& evaluation_point) {

    Vector diff = evaluation_point - src.position; 
    double r = diff.norm(); 

    if(r == 0) {
        return 0; 
    /* throw std::runtime_error("Cannot evaluate field at source location."); */
    }
    
    if constexpr(d == 2) { return src.q / (r*r) * diff; }
    { return src.q / (r*r*r) * diff; }
}

template<typename Vector, typename Source, std::size_t d>
Vector gravitationalForce_s(const Source& src, const Vector& evaluation_point) {
    return - electrostaticForce(src, evaluation_point); 
}

} //namespace fields

// InteractionResult can be e.g. double, Vector for potentials and forces
template<typename Vector, typename Source, typename InteractionResult> 
InteractionResult evaluateInteraction(std::vector<Source>& sources, 
        const Vector& evaluation_point, std::function <InteractionResult 
        (const Source&, const Vector&)> interactionFunction) {

    InteractionResult res{}; // Default init must be '0' element of result type
    res = std::accumulate(sources.begin(), sources.end(), res, 
        [&interactionFunction, &evaluation_point](InteractionResult& accumulant, 
        Source& src) -> InteractionResult { 
            return accumulant + interactionFunction(src, evaluation_point); 
        }
    );

    return res; 
}
} //namespace fmm

#endif
