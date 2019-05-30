#ifndef BALANCED_POINT_ORTHTREE_H
#define BALANCED_POINT_ORTHTREE_H

#include "balanced_orthtree.hpp"

template<typename Vector, std::size_t d>
class BalancedPointOrthtree: public BalancedOrthtree<Vector, d> {

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = BalancedOrthtree<Vector, d>;
    using BaseNode = typename AOT::BaseNode; 
    using AbstractOrthtree<Vector, d>::height;

public:
    struct Node; 
    struct Leaf; 

    BalancedPointOrthtree(std::vector<Vector> points, std::size_t s);

    void toFile() override; 

    ~BalancedPointOrthtree() { delete this->root; }

};

template<typename Vector, std::size_t d>
struct BalancedPointOrthtree<Vector, d>::Node: BaseNode {

    Node(Vector center, double box_length, std::size_t depth, 
        Node * parent): BaseNode(center, box_length, depth, parent) {}

    virtual ~Node() { for(auto child : this->children) { delete child; } }
};

template<typename Vector, std::size_t d>
struct BalancedPointOrthtree<Vector, d>::Leaf: Node {

    std::vector<Vector> * points;

    Leaf(Vector center, double box_length, std::size_t depth, 
        Node * parent): Node(center, box_length, depth, parent), points() {}
    virtual ~Leaf() { delete points; }
};

template<typename Vector, std::size_t d>
BalancedPointOrthtree<Vector, d>::BalancedPointOrthtree(std::vector<Vector> points, 
        std::size_t items_per_leaf): BalancedOrthtree<Vector, d>() {

    // Determine tree height, bounding box lenghts and center    
    this->height = ceil(log((double)points.size()/items_per_leaf) / 
            log(AOT::n_children));

    auto [lower_bounds, upper_bounds] = AOT::getDataRange(points);
    auto extents = Vector{upper_bounds - lower_bounds}.data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());

    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    std::size_t child_depth = 0; 

    this->root = new Node(center, box_length, child_depth, nullptr); 
    std::queue<Node*> node_queue({static_cast<Node*>(this->root)}); 

    while(child_depth < this->height) {

        std::size_t n_nodes = node_queue.size(); 
        assert(n_nodes == pow(AOT::n_children, child_depth));
        
        ++child_depth;

        for(std::size_t i = 0; i < n_nodes; i++) {

            Node * parent = node_queue.front();
            node_queue.pop(); 

            Vector parent_center = parent->center;
            double child_box_length = parent->box_length/2;

            for(int j = 0; j < AOT::n_children; j++) {

                Node * child;
                Vector child_center = parent_center + 
                    child_box_length/2 * AOT::child_center_directions[j];

                if(child_depth < this->height) {
                    child = new Node(child_center, child_box_length, 
                        child_depth, parent); 
                }
                else { // Treat leaves seperately
                    child = new Leaf(child_center, child_box_length, 
                        child_depth, parent); 
                }

                node_queue.push(child);
                parent->children[j] = child; 
            }
        }
    }

    // Distribute points to leaves
    std::size_t n_leaves = node_queue.size(); 
    assert(n_leaves == pow(AOT::n_children, child_depth));

    std::vector<Vector> ** leaf_vectors = new std::vector<Vector>*[n_leaves];
    for(std::size_t i = 0; i < n_leaves; i++) {
        leaf_vectors[i] = new std::vector<Vector>;
    }

    for(std::size_t k = 0; k < points.size(); k++) {
        std::array<uint32_t, d> indices = Super::getLeafBoxIndices(points[k]);
        assert(Super::getFlatIndex(indices) < n_leaves);
        leaf_vectors[Super::getFlatIndex(indices)]->push_back(points[k]);
    }

    for(std::size_t k = 0; k < n_leaves; k++) {

        Leaf * leaf = static_cast<Leaf*>(node_queue.front());
        node_queue.pop(); 

        std::array<uint32_t, d> indices = Super::getLeafBoxIndices(leaf->center);
        leaf->points = leaf_vectors[Super::getFlatIndex(indices)];
    }
     
    delete[] leaf_vectors; 
}

template<typename Vector, std::size_t d>
void BalancedPointOrthtree<Vector, d>::toFile() {

    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";

    std::ofstream geometry_file, data_file; 
    geometry_file.open(geometry_filename);
    data_file.open(data_filename);

    std::size_t n_node = 0;

    AOT::traverseBFSCore([&](BaseNode * current) {

        Vector center = current->center;
        double box_length = current->box_length;

        geometry_file << n_node << ", " << current->depth << ", " << box_length;
        for(auto coord : center.data()) { geometry_file << ", " << coord; }
        geometry_file << std::endl;

        // For leaf nodes, write contained points to file:
        if(current->children[0] == nullptr) { 

            Leaf * leaf = static_cast<Leaf*>(current);
            std::vector<Vector> * points = leaf->points; 

            data_file << n_node;
            for(Vector& p : *points) {
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
