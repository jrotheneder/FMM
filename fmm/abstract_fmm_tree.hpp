#ifndef ABSTRACT_FMM_TREE_H
#define ABSTRACT_FMM_TREE_H

#include "vector.hpp"
#include "multipole_expansion.hpp"
#include "local_expansion.hpp"
#include "fields.hpp"

namespace fmm {

template<std::size_t d, bool field_type> // true => Grav., false => Coulomb
class AbstractFmmTree {

protected:

    using Vector = Vector_<d>; 
    using PointSource = PointSource_<d>;

    using ME = MultipoleExpansion<d>;
    using LE = LocalExpansion<d>;

    const std::vector<PointSource>& sources;
    unsigned order;
    double eps; 
    double force_smoothing_eps; 

public: 
    AbstractFmmTree(std::vector<PointSource>& sources, double eps, 
            double force_smoothing_eps): sources(sources), eps(eps),
            force_smoothing_eps(force_smoothing_eps) {

        if constexpr(!(d == 2 || d == 3)) {
            throw std::logic_error("Implementations are only available for "
                "2 & 3 dimensions."); 
        } 

        // Determine expansion order
        double A = std::accumulate(sources.begin(), sources.end(), 0.0, 
            [](double acc, const PointSource& src) 
            { return acc + std::abs(src.sourceStrength()); }); 

        this->order = (ceil(log(A/eps) / log(2))); 
        //if constexpr(d == 2) {
        //  this->order = (ceil(log(A/eps) / log(2))); 
        //}
        //else if constexpr(d == 3) {
        //    this->order = (ceil(log(A/eps) / log(4./3)));  
        //}
    };

    virtual double evaluatePotential(const Vector& eval_point) const = 0; 
    virtual Vector evaluateForcefield(const Vector& eval_point) const = 0; 
    std::vector<double> evaluateParticlePotentialEnergies() const; 
    std::vector<Vector> evaluateParticleForces() const; 

    unsigned getOrder() { return order; }
    double getAccuracyEps() { return eps; }

    ~AbstractFmmTree() {}

protected:
    std::tuple<Vector, Vector> getDataRange() const;

    template<typename FmmNode> void localToLocal(FmmNode& node);
    template<typename FmmNode> void multipoleToLocal(FmmNode& node);
};

template<std::size_t d, bool field_type>
std::vector<double> AbstractFmmTree<d, field_type>::
        evaluateParticlePotentialEnergies() const {

    std::vector<double> potentials(sources.size()); 
    #pragma omp parallel for schedule(dynamic)
    for(std::size_t i = 0; i < sources.size(); ++i) {
        potentials[i] = sources[i].sourceStrength() 
            * evaluatePotential(sources[i].position);
    }

    return potentials; 
}

template<std::size_t d, bool field_type>
std::vector<Vector_<d>> AbstractFmmTree<d, field_type>::
        evaluateParticleForces() const {

    std::vector<Vector> forces(sources.size()); 

    #pragma omp parallel for schedule(dynamic)
    for(std::size_t i = 0; i < sources.size(); ++i) {
        forces[i] = sources[i].sourceStrength() 
            * evaluateForcefield(sources[i].position);
    }
    
    return forces; 
}

template<std::size_t d, bool field_type>
std::tuple<Vector_<d>, Vector_<d>> AbstractFmmTree<d, field_type>::getDataRange() const {

    Vector lower_bounds, upper_bounds; 
    lower_bounds.fill(HUGE_VAL);
    upper_bounds.fill(-HUGE_VAL);

    std::for_each(sources.begin(), sources.end(), [&](PointSource p) { 
        Vector& pos = p.position; 
        for(std::size_t i = 0; i < d; i++) {
            lower_bounds[i] = pos[i] < lower_bounds[i] ? pos[i] : lower_bounds[i];
            upper_bounds[i] = pos[i] > upper_bounds[i] ? pos[i] : upper_bounds[i];  
        }
    }); 

    // Expand bounding box s.t. all points lie completely within it (and not
    // on boundaries). This simplifies the index calculations in getLeafBoxIndices().
    double paddingFactor = 1E-5;
    lower_bounds = lower_bounds + paddingFactor * lower_bounds;
    upper_bounds = upper_bounds + paddingFactor * upper_bounds;

    return std::make_tuple(lower_bounds, upper_bounds);
}

template<std::size_t d, bool field_type>
template<typename FmmNode>
void AbstractFmmTree<d, field_type>::localToLocal(FmmNode& node) {
 
    FmmNode * parent = static_cast<FmmNode*>(node.parent);

    assert(parent != nullptr);
    assert(parent->local_expansion.coefficients.size() > 0);   

    node.local_expansion += LE(node.center, parent->local_expansion);
}


template<std::size_t d, bool field_type>
template<typename FmmNode>
void AbstractFmmTree<d, field_type>::multipoleToLocal(FmmNode& node) {
 
    std::vector<const ME*> incoming; 
    for(FmmNode* interaction_partner : node.interaction_list) {
        incoming.push_back(&interaction_partner->multipole_expansion);  
    }

    if(incoming.size() > 0) {
        node.local_expansion += LE(node.center, incoming); 
    }
}

} // namespace fmm

#endif
