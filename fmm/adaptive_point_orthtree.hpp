#ifndef ADAPTIVE_POINT_ORTHTREE_H
#define ADAPTIVE_POINT_ORTHTREE_H

#include "adaptive_orthtree.hpp"

template<typename Vector, std::size_t d>
class AdaptivePointOrthtree: public AdaptiveOrthtree<Vector, d> {

protected: 

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = AdaptiveOrthtree<Vector, d>;
    using BaseNode = typename AOT::BaseNode; 
    using AbstractOrthtree<Vector, d>::height;

    struct Node; 
    struct Leaf; 

    std::vector<std::vector<Node*>> nodes;  // nodes at depth i stored in nodes[i]
    std::vector<std::vector<Leaf*>> leaves; // leaves at depth i stored in leaves[i]

public:

    AdaptivePointOrthtree(std::vector<Vector> points, std::size_t s);

    static std::array<std::vector<Vector>, AOT::n_children> refineOrthant(
        const Vector& center, const std::vector<Vector>& sources);
    void buildTree(Node* node, std::vector<Vector>& points, 
            std::size_t max_items_per_leaf);

    void toFile() override; 

    ~AdaptivePointOrthtree() { delete this->root; }

};

template<typename Vector, std::size_t d>
struct AdaptivePointOrthtree<Vector, d>::Node: BaseNode {

    Node(Vector center, double box_length, std::size_t depth, 
        Node * parent): BaseNode(center, box_length, depth, parent) {}

    virtual ~Node() { 
        for(auto child : this->children) { 
            if(child != nullptr) {
                delete child; 
            }
        } 
    }
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
        std::size_t max_items_per_leaf): Super(), leaves(1) {

    // Determine bounding box lenghts and center    

    auto [lower_bounds, upper_bounds] = AOT::getDataRange(points);
    auto extents = Vector{upper_bounds - lower_bounds}.data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    // TODO handle case where root is the only leaf
    std::size_t child_depth = 0; 
    this->height = 0; 

    Node * root = new Node(center, box_length, child_depth, nullptr);
    this->nodes.push_back(std::vector<Node*>{root});
    this->root = root; 

    buildTree(this->nodes[0][0], points, max_items_per_leaf); 
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
void AdaptivePointOrthtree<Vector, d>::buildTree(Node* node, 
        std::vector<Vector>& points, std::size_t max_items_per_leaf) {

    unsigned depth = 0; 
    bool refine = true; 

    //sources of nodes at depth which still need to be assigned to 
    //leaves further down in the tree. Unfortunately we incur quite a few copies
    //of large vectors here, but storage requirements are limited to 2N + overhead, 
    //and we are cpu bound anyways 
    std::vector<std::vector<Vector>> node_sources{points}; 
    // temporary storage for sources assigned to the nodes on the next level
    std::vector<std::vector<Vector>> temp_node_sources{}; 

    while(refine) { // construct the tree level by level

        ++depth;
        nodes.push_back(std::vector<Node*>()); 
        leaves.push_back(std::vector<Leaf*>()); 

        assert(temp_node_sources.empty());

        for(unsigned i = 0; i < nodes[depth-1].size(); ++i) {

            Node* node = nodes[depth-1][i];

            Vector parent_center = node->center; 
            double child_box_length = node->box_length/2; 
            unsigned child_depth = node->depth + 1;

            assert(child_depth == depth); 

            auto child_sources = refineOrthant(node->center, node_sources[i]); 

            for(std::size_t i = 0; i < AOT::n_children; ++i) {

                Vector child_center = parent_center + 
                    child_box_length/2 * AOT::child_center_directions[i];

                if(child_sources[i].size() > max_items_per_leaf) { // refine

                    Node* child = new Node(child_center, child_box_length, 
                                child_depth, node);
                    node->children[i] = child;

                    nodes[depth].push_back(child);
                    temp_node_sources.push_back(child_sources[i]);

                }
                else {
                    Leaf* child = new Leaf(child_center, child_box_length, 
                        child_depth, node, child_sources[i]);
                    node->children[i] = child;

                    leaves[depth].push_back(child);

                }
            }
        }

        node_sources.swap(temp_node_sources); 
        temp_node_sources.clear(); 
        
        assert(node_sources.size() == nodes[depth].size());
        refine = !nodes[depth].empty(); 
    }
    this->height = depth; 
}

template<typename Vector, std::size_t d>
void AdaptivePointOrthtree<Vector, d>::toFile() {

    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";

    std::ofstream geometry_file, data_file; 
    geometry_file.open(geometry_filename);
    data_file.open(data_filename);

    std::size_t n_node = 0;

    this-> template traverseBFSCore<BaseNode>([&](BaseNode * current) {

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
