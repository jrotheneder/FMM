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
    static const int n_children = pow(2,d); 
    const std::array<Vector, n_children> child_center_directions; 

    AbstractOrthtree(): child_center_directions(getChildCenterDirections()) {};

    void traverseBFSCore(const std::function <void(BaseNode *)>& processNode);

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

    virtual ~BaseNode() {}; // Implement tailored destructor in deriving classes
}; 


template<typename Vector, std::size_t d>
void AbstractOrthtree<Vector, d>::traverseBFSCore(const std::function 
        <void(BaseNode *)>& processNode) {

    std::queue<BaseNode*> node_queue({this->root}); 
    while(!node_queue.empty()) {

        BaseNode * current_node = node_queue.front();
        node_queue.pop(); 

        processNode(current_node); 

        for(int i = 0; i < n_children; i++) {
            if(current_node->children[i] != nullptr)  {
                node_queue.push(current_node->children[i]);
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

    vecToFile(child_center_directions, child_center_directions.size(), "test.dat"); 
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


    return make_tuple(lower_bounds, upper_bounds);
}

template<typename Vector, std::size_t d>
AbstractOrthtree<Vector, d>::~AbstractOrthtree() {}

#endif
