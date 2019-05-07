#ifndef BALANCED_ORTHTREE_H
#define BALANCED_ORTHTREE_H

#include "abstract_orthtree.hpp"

template<typename Vector, std::size_t d>
class BalancedOrthtree: public AbstractOrthtree<Vector, d> {

protected:
    using Super = AbstractOrthtree<Vector, d>;
    using Node = typename Super::Node; 

public:
    std::array<std::size_t, d> getLeafBoxIndices(Vector p) const;
    std::size_t getFlatIndex(const std::array<std::size_t, d>& indices) const;

    BalancedOrthtree(): AbstractOrthtree<Vector, d>() {}
    virtual ~BalancedOrthtree() = 0; 
};

//
// Given a point, determines indices (i_1,i_2,...) of leaf that this point belongs
// to, where leaves are indexed according to the part of they domain they
// own (smallest x coordinate -> i_1 = 0, smallest y coordinate -> i_2 = 0 etc...) 
template<typename Vector, std::size_t d>
std::array<std::size_t, d> BalancedOrthtree<Vector, d>::getLeafBoxIndices(
        Vector point) const {

    Node * root = this->root;
    double box_length = root->box_length; 

    Vector diagonal;
    diagonal.fill(1);

    Vector lower_left_box_corner = root->center - box_length/2 * diagonal; 
    Vector ratios = (point - lower_left_box_corner) * (1 / box_length);

    const std::size_t n_boxes_per_dim = pow(2, this->height); 

    std::array<std::size_t, d> indices; 
    for(std::size_t j = 0; j < d; j++) {
        indices[j] = floor(ratios[j] * n_boxes_per_dim);
    }
    
    return indices;
}


template<typename Vector, std::size_t d>
std::size_t BalancedOrthtree<Vector, d>::getFlatIndex(
        const std::array<std::size_t, d>& indices) const {
    
    const std::size_t n_boxes_per_dim = pow(2, this->height);   

    std::size_t flat_index = 0; 
    std::size_t offset_factor = 1; 

    for(std::size_t i = 0; i < d; ++i) {
        flat_index += offset_factor * indices[d-i-1];  
        offset_factor *= n_boxes_per_dim; 
    }

    return flat_index;
}

template<typename Vector, std::size_t d>
BalancedOrthtree<Vector, d>::~BalancedOrthtree() {}

#endif
