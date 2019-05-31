#ifndef ADAPTIVE_ORTHTREE_H
#define ADAPTIVE_ORTHTREE_H

#include "abstract_orthtree.hpp"

template<typename Vector, std::size_t d>
class AdaptiveOrthtree: public AbstractOrthtree<Vector, d> {

protected:
    using Super = AbstractOrthtree<Vector, d>;
    using Node = typename Super::BaseNode; 

public:
    AdaptiveOrthtree(): AbstractOrthtree<Vector, d>() {}
    
    static unsigned getOrthant(const Vector& center, const Vector& point); 

    virtual ~AdaptiveOrthtree() = 0; 
};

// Computes index of orthant that a point falls into w.r.t. a given center. 
// This method assumes a specific ordering of the orthants, i.e. it will break
// if AbstractOrthtree::getChildCenterDirections() is modified. 
template<typename Vector, std::size_t d>
unsigned AdaptiveOrthtree<Vector, d>::getOrthant(const Vector& center, 
        const Vector& point) {

    unsigned orthantIndex = 0; 
    for(unsigned i = 0; i < d; ++i) {
        orthantIndex <<= 1; 
        orthantIndex |= (point[i] >= center[i] ? 1 : 0);    
    }

    return orthantIndex; 
}

template<typename Vector, std::size_t d>
AdaptiveOrthtree<Vector, d>::~AdaptiveOrthtree() {}

#endif
