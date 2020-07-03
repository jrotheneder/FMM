#ifndef BALANCED_ORTHTREE_H
#define BALANCED_ORTHTREE_H

#include "abstract_orthtree.hpp"

template<typename Vector, std::size_t d>
class BalancedOrthtree: public AbstractOrthtree<Vector, d> {

protected:
    using Super = AbstractOrthtree<Vector, d>;
    using Node = typename Super::BaseNode; 

public:
    std::array<uint32_t, d> getLeafBoxIndices(Vector p) const;
    uint32_t getFlatIndex(const std::array<uint32_t, d>& indices) const;
    unsigned long long getMortonIndex(const std::array<uint32_t, d>& indices) const;

    BalancedOrthtree(): AbstractOrthtree<Vector, d>() {}
    virtual ~BalancedOrthtree() = 0; 
};

// Given a point, determines indices (i_1,i_2,...) of leaf that this point belongs
// to, where leaves are indexed according to the part of they domain they
// own (smallest x coordinate -> i_1 = 0, smallest y coordinate -> i_2 = 0 etc...) 
template<typename Vector, std::size_t d>
std::array<uint32_t, d> BalancedOrthtree<Vector, d>::
        getLeafBoxIndices(Vector point) const {

    Node * root = this->root;
    double box_length = root->box_length; 

    Vector diagonal(1);

    Vector lower_left_box_corner = root->center - box_length/2 * diagonal; 
    Vector ratios = (point - lower_left_box_corner) * (1 / box_length);

    uint32_t n_boxes_per_dim = pow(2, this->height); 
    std::array<uint32_t, d> indices; 

    for(unsigned i = 0; i < d; i++) {
        indices[i] = floor(ratios[i] * n_boxes_per_dim);

        if(indices[i] >= n_boxes_per_dim) {
            throw std::runtime_error("Illegal index: point lies outside of "
                "tree region."); 
        }
    }
    
    return indices;
}

template<typename Vector, std::size_t d>
uint32_t BalancedOrthtree<Vector, d>::
        getFlatIndex(const std::array<uint32_t, d>& indices) const {
    
    const uint32_t n_boxes_per_dim = pow(2, this->height);   

    uint32_t flat_index = 0; 
    uint32_t offset_factor = 1; 

    for(unsigned i = 0; i < d; ++i) {
        flat_index += offset_factor * indices[d-i-1];  
        offset_factor *= n_boxes_per_dim; 
    }

    return flat_index;
}

// Compute 64-bit morton index from d-tuple of 32-bit unsigned indices 
// Interleaving order x (least significant), y, z(most significant) 
// is used here.
template<typename Vector, std::size_t d>
unsigned long long BalancedOrthtree<Vector, d>::
        getMortonIndex(const std::array<uint32_t, d>& indices) const {

    using Morton_t = unsigned long long;  // Datatype of returned morton index
    
    const unsigned morton_index_size = sizeof(Morton_t) * 8;  
    const unsigned max_bits_per_index = morton_index_size/d; 

    // (2^(d*h) boxes at height h require 2^h indices per dimension for indexing
    // which in turn require h bits each => the individual indices may only use 
    // h bits, else we overflow. 
    const unsigned bits_used = this->height; 

    if(bits_used >= max_bits_per_index) {
        throw std::runtime_error("ULL size insufficient for given tree."); 
    }; 

    std::bitset<morton_index_size> morton_index; 
    std::bitset<morton_index_size> helper_bitset; 

    for(unsigned i = 0; i < d; ++i) {
        helper_bitset <<= bits_used; 
        helper_bitset |= std::bitset<morton_index_size>(indices[i]); 
    }

    unsigned index_counter = 0; 
    for(unsigned i = 0; i < bits_used; ++i) {
        for(unsigned j = 0; j < d; ++j) {
            morton_index[index_counter++] = helper_bitset[j * bits_used + i];  
        }
         
    }

    return morton_index.to_ullong(); // throws if overflow
}

template<typename Vector, std::size_t d>
BalancedOrthtree<Vector, d>::~BalancedOrthtree() {}

#endif
