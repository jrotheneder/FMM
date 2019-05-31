#ifndef ABSTRACT_FMM_TREE_H
#define ABSTRACT_FMM_TREE_H

#include "multipole_expansion.hpp"
#include "local_expansion.hpp"
#include "direct.hpp"

namespace fmm {

template<std::size_t d>
class AbstractFmmTree {

protected:

    using Vec = Vector<d>; 
    using Source = PointSource<d>;

    using AOT = AbstractOrthtree<Vec, d>;
    using Super = BalancedOrthtree<Vec, d>;
    using BaseNode = typename AOT::BaseNode; 
    using ME = MultipoleExpansion<Vec, Source, d>;
    using LE = LocalExpansion<Vec, Source, d>;


    static constexpr auto evalScalarInteraction 
        = evaluateInteraction<Vec, Source, double>;
    static constexpr auto evalVecInteraction 
        = evaluateInteraction<Vec, Source, Vec>;

    static constexpr auto potentialFunction
        = fields::electrostaticPotential<Vec, Source, d>;
    static constexpr auto safePotentialFunction
        = fields::electrostaticPotential_s<Vec, Source, d>;
    static constexpr auto forceFunction
        = fields::electrostaticForce<Vec, Source, d>;
    static constexpr auto safeForceFunction
        = fields::electrostaticForce_s<Vec, Source, d>;

    const double max_neighbour_distance = 1.1 * sqrt(d); 
    //in units of box_length + padding to avoid numerical issues

    std::vector<Source>& sources;
    std::size_t order;

public: 
    AbstractFmmTree(std::vector<Source>& sources): sources(sources) {};

    virtual double evaluatePotential(const Vec& eval_point) const = 0; 
    virtual Vec evaluateForcefield(const Vec& eval_point) const = 0; 
    std::vector<double> evaluateParticlePotentials() const; 
    std::vector<Vec> evaluateParticleForcefields() const; 

    virtual ~AbstractFmmTree() {}

protected:
    static std::tuple<Vec, Vec> getDataRange(const std::vector<Source>& sources);

    template<typename FmmNode>
    void localToLocal(FmmNode& node);
    template<typename FmmNode>
    void multipoleToLocal(FmmNode& node);
};

template<std::size_t d>
std::vector<double> AbstractFmmTree<d>::evaluateParticlePotentials() const {

    std::vector<double> potentials(sources.size()); 
    for(std::size_t i = 0; i < sources.size(); ++i) {
        potentials[i] = evaluatePotential(sources[i].position);
    }

    return potentials; 
}

template<std::size_t d>
std::vector<Vector<d>> AbstractFmmTree<d>::
        evaluateParticleForcefields() const {

    std::vector<Vec> forces(sources.size()); 
    for(std::size_t i = 0; i < sources.size(); ++i) {
        forces[i] = evaluateForcefield(sources[i].position);
    }
    
    return forces; 
}

template<std::size_t d>
std::tuple<Vector<d>, Vector<d>> AbstractFmmTree<d>::getDataRange(
        const std::vector<Source> & sources) {

    Vec lower_bounds, upper_bounds; 
    lower_bounds.fill(HUGE_VAL);
    upper_bounds.fill(-HUGE_VAL);

    std::for_each(sources.begin(), sources.end(), [&](Source p) { 
        Vec & pos = p.position; 
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

template<std::size_t d>
template<typename FmmNode>
void AbstractFmmTree<d>::localToLocal(FmmNode& node) {
 
    FmmNode * parent = static_cast<FmmNode*>(node.parent);

    assert(parent != nullptr);
    assert(parent->local_expansion.coefficients.size() > 0);   

    node.local_expansion += LE(node.center, parent->local_expansion);
}


template<std::size_t d>
template<typename FmmNode>
void AbstractFmmTree<d>::multipoleToLocal(FmmNode& node) {
 
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
