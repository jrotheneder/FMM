#include <algorithm> 
#include <functional> 

namespace direct {

template<typename Vector, typename Source, std::size_t d>
double Electrostatic_Potential(const Source& src, const Vector& evaluation_point) {

    double r = (src.position - evaluation_point).norm(); 

    if constexpr(d == 2) { return -src.q * log(r); }
    else { return src.q / r; }
}

template<typename Vector, typename Source, std::size_t d>
double Gravitational_Potential(const Source& src, const Vector& evaluation_point) {
    return -Electrostatic_Potential(src, evaluation_point); 
}

template<typename Vector, typename Source, std::size_t d>
Vector Electrostatic_Force(const Source& src, const Vector& evaluation_point) {

    Vector diff = evaluation_point - src.position; 
    double r = diff.norm(); 

    if constexpr(d == 2) { return src.q / (r*r) * diff; }
    else { return src.q / (r*r*r) * diff; }
}

template<typename Vector, typename Source, std::size_t d>
Vector Gravitational_Force(const Source& src, const Vector& evaluation_point) {
    return - Electrostatic_Force(src, evaluation_point); 
}

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

} //namespace direct
