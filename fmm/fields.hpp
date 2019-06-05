#include <cmath> 

#ifndef FMM_FIELDS_H
#define FMM_FIELDS_H

namespace fmm {
namespace fields {

// Gravitational (if grav = true) or Coulomb (if grav = false) potential 
// by a single source in d \in {2,3} dimensions. With safe = true, 
// this function checks whether the evaluation point 
// coincides with the source location and returns zero in that case
template<std::size_t d, bool grav = true, bool safe = true>
double potential(const PointSource_<d>& src, const Vector_<d>& eval_point, 
        const double eps = 0) {

    double r = (src.position - eval_point).norm(); 

    if constexpr(safe) { if(r == 0) { return 0; }}

    double pot; 
    if constexpr(d == 2) { pot = src.q * std::log(r + eps); }
    else { pot = src.q / (r + eps); }

    if constexpr(grav) { return pot; }
    else { return -pot; }  // Coulomb potential
}

// Overload of the potential function for multiple sources
template<std::size_t d, bool grav = true, bool safe = true>
double potential(const std::vector<PointSource_<d>>& sources, 
        const Vector_<d>& eval_point, const double eps = 0) {

    double pot = 0; 

    for(std::size_t i = 0; i < sources.size(); ++i) {
        pot += potential<d, grav, safe>(sources[i], eval_point, eps);  
    }

    return pot; 
}

// Gravitational (if grav = true) or Coulomb (if grav = false) force field 
// by a single source in d \in {2,3} dimensions. With safe = true, this 
// function checks whether the evaluation point coincides with the source 
// location and returns the zero vector in that case
template<std::size_t d, bool grav = true, bool safe = true>
Vector_<d> forcefield(const PointSource_<d>& src, const Vector_<d>& eval_point, 
        const double eps = 0) {

    Vector_<d> diff = src.position - eval_point; 

    if constexpr(d == 2) {
        double r_sq = diff.norm_sq(); 
        if constexpr(safe) { if(r_sq == 0) { return Vector_<d>{}; }}
        diff *= src.q/(r_sq + eps); 
    } 
    else {
        double r = diff.norm(); 
        if constexpr(safe) { if(r == 0) { return Vector_<d>{}; }}
        diff *= src.q / (r*r*r + eps); 
    }

    if constexpr(grav) { return diff; }
    else { return -diff; }
}

// Overload of the force function for multiple sources. 
template<std::size_t d, bool grav = true, bool safe = true>
Vector_<d> forcefield(const std::vector<PointSource_<d>>& sources, 
        const Vector_<d>& eval_point, const double eps = 0) {

    Vector_<d> frc{};

    for(std::size_t i = 0; i < sources.size(); ++i) {
        frc += forcefield<d, grav, safe>(sources[i], eval_point, eps);  
    }

    return frc; 
}

} //namespace fields
} //namespace fmm

#endif
