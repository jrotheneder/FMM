#ifndef ORTHTREE_H
#define ORTHTREE_H

#include <vector> 
#include <cmath> 
#include <algorithm> 
#include <array> 
#include <queue> 
#include <stack> 
#include <fstream> 
#include <iostream> 
#include <string> 
#include <sstream> 
#include <bitset> 
#include <functional> 
#include <cassert> 

#include "debugging.hpp" 
#include "vector.hpp"

template<typename Vector, std::size_t d>
class AbstractOrthtree {
public:
    
    struct BaseNode; 

    BaseNode * root; 

    std::size_t height;
    static const int n_children = 1 << d; 
    const std::array<Vector, n_children> child_center_directions; 

    AbstractOrthtree(): child_center_directions(getChildCenterDirections()) {};

    template<typename Node>
    void traverseBFSCore(const std::function <void(Node *)>& processNode);

    std::size_t getHeight() const { return this->height; }
    Vector getCenter() const { return this->root->center; }
    double getBoxLength() const { return this->root->box_length; }

    virtual void toFile() = 0;

    static std::array<Vector, n_children> getChildCenterDirections();
    static std::tuple<Vector, Vector> getDataRange(const std::vector<Vector> & points);

    virtual ~AbstractOrthtree() = 0; 

};


template<typename Vector, std::size_t d>
struct AbstractOrthtree<Vector, d>::BaseNode {

    Vector center;         
    double box_length;     // This nodes's bounding box length
    std::size_t depth;     // Depth of the node in the tree

    BaseNode * parent; 
    std::array<BaseNode*, n_children> children; 

    BaseNode() {} // Empty default constructor
    BaseNode(Vector center, double box_length, std::size_t depth, 
        BaseNode * parent = nullptr): center(center), box_length(box_length), 
        depth(depth), parent(parent), children() {}

    bool adjacent(BaseNode* other) const; 
    bool contains(const Vector& point) const; 

    virtual ~BaseNode() {}; // Implement tailored destructor in deriving classes
}; 

template<typename Vector, std::size_t d>
bool AbstractOrthtree<Vector, d>::BaseNode::adjacent(BaseNode* other) const {

    Vector ones(1); 
    Vector self_lower = this->center - this->box_length/2 * ones; 
    Vector self_upper = this->center + this->box_length/2 * ones; 
    Vector other_lower = other->center - other->box_length/2 * ones; 
    Vector other_upper = other->center + other->box_length/2 * ones; 

    // We check whether the regions owned by two nodes intersect by checking 
    // overlap in each dimension separately. In 1D, two lines (x1, x2), (y1, y2)
    // do not intersect if x2 < y1 || x1 > y2 or equiv. if y1-x2 > 0 && x1 - y2
    // > 0. To allow for some numerical error, we instead check whether 
    // y1 - x2 - 1/2 length(min(box1, box2)) > 0 etc, which is possible to due
    // to regular layout of boxes in this setting.

    Vector tolerance(std::min(this->box_length, other->box_length)/2);
    Vector diffs1 = other_lower - self_upper - tolerance; 
    Vector diffs2 = self_lower - other_upper - tolerance; 

    for(unsigned i = 0; i < d; ++i) {

        if(diffs1[i] > 0 || diffs2[i] > 0) { return false; }
    }

    return true; 
}


template<typename Vector, std::size_t d>
bool AbstractOrthtree<Vector, d>::BaseNode::contains(const Vector& point) const {

    Vector ones(1); 
    Vector lower = this->center - this->box_length/2 * ones; 
    Vector upper = this->center + this->box_length/2 * ones; 

    for(unsigned i = 0; i < d; ++i) {
        if(point[i] < lower[i] || point[i] > upper[i]) { return false; }
    }

    return true; 
}

template<typename Vector, std::size_t d>
template<typename Node>
void AbstractOrthtree<Vector, d>::traverseBFSCore(const std::function 
        <void(Node *)>& processNode) {

    std::queue<Node*> node_queue({static_cast<Node*>(this->root)}); 

    while(!node_queue.empty()) {

        Node * current_node = node_queue.front();
        node_queue.pop(); 

        processNode(current_node); 

        for(int i = 0; i < n_children; i++) {
            if(current_node->children[i] != nullptr)  {
                node_queue.push(static_cast<Node*>(current_node->children[i]));
            }
        }
    }
}

template<typename Vector, std::size_t d>  
std::array<Vector, AbstractOrthtree<Vector, d>::n_children> 
        AbstractOrthtree<Vector, d>::getChildCenterDirections() { 

    std::array<Vector, AbstractOrthtree<Vector, d>::n_children> 
        child_center_directions;

    //Use a length d bit vector to represent any of the n_children = 2^d 
    //possible directions in which the centers of a given node's children 
    //may lie. A 0 in the d-j-1-st position of the bit vector corresponds to 
    //a direction vector whose j-th component is -1; a 1 in the same position 
    //to a direction vector whose j-th component is +1. 
    //E.g. bitset 01 goes to {-1,1}. 

    for(std::size_t i = 0; i < AbstractOrthtree<Vector, d>::n_children; ++i) {

        std::bitset<d> direction_bits = std::bitset<d>(i); 

        Vector v;  
        for(std::size_t j = 0; j < d; ++j) {
            v[j] = direction_bits[d-j-1] ? 1 : -1;   
        }

        child_center_directions[i] = v;  
    }

    // Sample outputs: 
    // 2D: {-1, -1}, {-1, 1}, {1, -1}, {1, 1}
    // 3D: {-1, -1, -1}, {-1, -1, 1 }, {-1, 1, -1}, {-1, 1, 1}, 
    //     {1, -1, -1}, {1, -1, 1 }, {1, 1, -1}, {1, 1, 1}

    return child_center_directions; 
}

// Return tuple of (x_min, x_max, y_min, y_max, ...) of set of input vectors
// (x, y, ...)
template<typename Vector, std::size_t d>
std::tuple<Vector, Vector> AbstractOrthtree<Vector, d>::getDataRange(
        const std::vector<Vector> & points) {

    Vector lower_bounds, upper_bounds; 
    lower_bounds.fill(HUGE_VAL);
    upper_bounds.fill(-HUGE_VAL);

    std::for_each(points.begin(), points.end(), [&](Vector p) { 
        for(std::size_t i = 0; i < d; i++) {
            lower_bounds[i] = p[i] < lower_bounds[i] ? p[i] : lower_bounds[i];  
            upper_bounds[i] = p[i] > upper_bounds[i] ? p[i] : upper_bounds[i];  
        }
    }); 

    // Expand bounding box s.t. all points lie completely within it (and not
    // on boundaries. This simplifies the index calculations in getLeafBoxIndices().

    double paddingFactor = 1E-5;
    lower_bounds = lower_bounds + paddingFactor * lower_bounds;
    upper_bounds = upper_bounds + paddingFactor * upper_bounds;


    return std::make_tuple(lower_bounds, upper_bounds);
}

template<typename Vector, std::size_t d>
AbstractOrthtree<Vector, d>::~AbstractOrthtree() {}

#endif
