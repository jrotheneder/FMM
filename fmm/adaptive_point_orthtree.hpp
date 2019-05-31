#ifndef ADAPTIVE_POINT_ORTHTREE_H
#define ADAPTIVE_POINT_ORTHTREE_H

#include "adaptive_orthtree.hpp"

template<typename Vector, std::size_t d>
class AdaptivePointOrthtree: public AdaptiveOrthtree<Vector, d> {

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = AdaptiveOrthtree<Vector, d>;
    using BaseNode = typename AOT::BaseNode; 
    using AbstractOrthtree<Vector, d>::height;

public:
    struct Node; 
    struct Leaf; 

    AdaptivePointOrthtree(std::vector<Vector> points, std::size_t s);

    static std::array<std::vector<Vector>, AOT::n_children> refineOrthant(
        const Vector& center, const std::vector<Vector>& sources);
    void splitNode(Node* node, std::vector<Vector> points, 
            std::size_t max_items_per_leaf);

    void toFile() override; 

    ~AdaptivePointOrthtree() { delete this->root; }

};

template<typename Vector, std::size_t d>
struct AdaptivePointOrthtree<Vector, d>::Node: BaseNode {

    Node(Vector center, double box_length, std::size_t depth, 
        Node * parent): BaseNode(center, box_length, depth, parent) {}

    virtual ~Node() { for(auto child : this->children) { delete child; } }
};

template<typename Vector, std::size_t d>
struct AdaptivePointOrthtree<Vector, d>::Leaf: Node {

    std::vector<Vector> points;

    Leaf(Vector center, double box_length, std::size_t depth, 
        Node * parent, std::vector<Vector> points): 
        Node(center, box_length, depth, parent), points(points) {}

    virtual ~Leaf() {}
};

template<typename Vector, std::size_t d>
AdaptivePointOrthtree<Vector, d>::AdaptivePointOrthtree(std::vector<Vector> points, 
        std::size_t max_items_per_leaf): Super() {

    // Determine bounding box lenghts and center    

    auto [lower_bounds, upper_bounds] = AOT::getDataRange(points);
    auto extents = Vector{upper_bounds - lower_bounds}.data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    std::size_t child_depth = 0; 
    this->height = 0; 
    this->root = new Node(center, box_length, child_depth, nullptr); 

    splitNode(static_cast<Node*>(this->root), points, max_items_per_leaf); 
}

template<typename Vector, std::size_t d>
std::array<std::vector<Vector>, AbstractOrthtree<Vector, d>::n_children> 
        AdaptivePointOrthtree<Vector, d>::refineOrthant(const Vector& center, 
        const std::vector<Vector>& sources) {

    std::array<std::vector<Vector>, AOT::n_children> child_sources; 

    for(auto& source : sources) {
        unsigned index = Super::getOrthant(center, source); 
        child_sources[index].push_back(source);   
    }
    
    return child_sources; 
}

template<typename Vector, std::size_t d>
void AdaptivePointOrthtree<Vector, d>::splitNode(Node* node, 
        std::vector<Vector> points, std::size_t max_items_per_leaf) {

    Vector parent_center = node->center; 
    double child_box_length = node->box_length/2; 
    unsigned child_depth = node->depth + 1;

    if(child_depth >= this->height) {
        this->height = child_depth; 
    }

    auto child_sources = refineOrthant(node->center, points); 

    for(std::size_t i = 0; i < AOT::n_children; ++i) {

        Vector child_center = parent_center + 
            child_box_length/2 * AOT::child_center_directions[i];

        if(child_sources[i].size() > max_items_per_leaf) { // child needs refinement

            Node * child = new Node(child_center, child_box_length, 
                        child_depth, node);
            node->children[i] = child; 
            splitNode(child, child_sources[i], max_items_per_leaf); 
        }
        else {
            node->children[i] = new Leaf(child_center, child_box_length, 
                        child_depth, node, child_sources[i]);  
        }
    }
}

template<typename Vector, std::size_t d>
void AdaptivePointOrthtree<Vector, d>::toFile() {

    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";

    std::ofstream geometry_file, data_file; 
    geometry_file.open(geometry_filename);
    data_file.open(data_filename);

    std::size_t n_node = 0;

    AOT::traverseBFSCore([&](BaseNode * current) {

        Vector center = current->center;
        double box_length = current->box_length;
        bool is_leaf = current->children[0] == nullptr;

        geometry_file << n_node << ", " << is_leaf << "," << current->depth 
            << ", " << box_length;
        for(auto coord : center.data()) { geometry_file << ", " << coord; }
        geometry_file << std::endl;

        // For leaf nodes, write contained points to file:
        if(current->children[0] == nullptr) { 

            Leaf * leaf = static_cast<Leaf*>(current);
            std::vector<Vector>& points = leaf->points; 

            data_file << n_node;
            for(Vector& p : points) {
                for(double coord : p.data()) {
                    data_file << ", " << coord;         
                }
            }
            data_file << "\n";
        }
        ++n_node; 
    });
    
    geometry_file.close();
    data_file.close();
}

#endif
