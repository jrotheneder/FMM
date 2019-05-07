#ifndef FMM_TREE_H
#define FMM_TREE_H

#include "balanced_orthtree.hpp"
#include "multipole_expansion.hpp"
#include "local_expansion.hpp"

namespace fmm {

template<typename Vector, typename Source, std::size_t d>
class BalancedFmmTree: public BalancedOrthtree<Vector, d> {


    using AOT = AbstractOrthtree<Vector, d>;
    using Super = BalancedOrthtree<Vector, d>;

    using BaseNode = typename AOT::Node; 

    struct FmmNode; 
    struct FmmLeaf; 

    std::vector<Source>& sources;

public: 

    static constexpr double max_neighbour_distance = 1.1 * sqrt(d); 
    //in units of box_length + padding to avoid numerical issues

    BalancedFmmTree(std::vector<Source>& sources, std::size_t items_per_cell, 
            double eps);

    std::vector<double> computePotentials(); 
    std::vector<Vector> computeForcefields(); 

    void traverseBFSCore(const std::function <void(FmmNode *)>& processNode);
    void toFile() override;

    static std::tuple<Vector, Vector> getDataRange(
        const std::vector<Source>& sources);
    void expandNode(std::queue<FmmNode*>& leaf_queue, FmmNode* node);
    void distributeSources(std::queue<FmmNode*>& leaf_queue, 
        const std::vector<Source>& sources);
    void computeNodeNeighbourhood(FmmNode* node); 


    ~BalancedFmmTree() { delete this->root; }
};

template<typename Vector, typename Source, std::size_t d>
struct BalancedFmmTree<Vector, Source, d>::FmmNode: AbstractOrthtree<Vector, d>::Node {

protected:
    using Super = typename AbstractOrthtree<Vector, d>::Node; 
    using ME = MultipoleExpansion<Vector, Source, d>;
    using LE = LocalExpansion<Vector, Source, d>;
public:

    ME multipole_expansion; 
    LE local_expansion;

    std::vector<FmmNode*> interaction_list;
    std::vector<FmmNode*> near_neighbours;

    FmmNode(Vector center, double box_length, std::size_t depth, 
        FmmNode * parent = nullptr): Super(center, box_length, depth, parent) {};

    virtual ~FmmNode() {};
};

template<typename Vector, typename Source, std::size_t d>
struct BalancedFmmTree<Vector, Source, d>::FmmLeaf: 
        BalancedFmmTree<Vector, Source, d>::FmmNode {

    using Super = BalancedFmmTree<Vector, Source, d>::FmmNode;
    using BaseNode = typename Super::Node; 

    std::vector<Source> * sources;

    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent = nullptr): Super(center, box_length, depth, parent),
        sources() {};

    virtual ~FmmLeaf() { delete sources; }
};

template<typename Vector, typename Source, std::size_t d>
BalancedFmmTree<Vector, Source, d>::BalancedFmmTree(std::vector<Source>& sources, 
        std::size_t items_per_cell, double eps): BalancedOrthtree<Vector, d>(), 
        sources(sources) { 

    if(d > 3lu) { 
        throw std::runtime_error("Determining near neighbours from "
            "euclidean center separation is only valid for d <= 3 "
            "dimensions, modify implementation."); 
    }

    // Determine tree height, bounding box lenghts and center 
    this->height = ceil(
        log((double)sources.size()/items_per_cell) / log(AOT::n_children)
    );

    Vector lower_bounds, upper_bounds;  
    std::tie(lower_bounds, upper_bounds) = getDataRange(sources); 
    auto extents = (upper_bounds - lower_bounds).data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    this->root = new FmmNode(center, box_length, 0, nullptr); 
    std::queue<FmmNode*> leaf_queue; 

    traverseBFSCore( // Init children for every node down to the desired height
        [this, &leaf_queue] (FmmNode * node) -> void 
        { this->expandNode(leaf_queue, node); }
    ); 

    distributeSources(leaf_queue, sources); // Distribute sources to leaves

    traverseBFSCore( // Compute neighbourhoods of all nodes
        [this](FmmNode* node) -> void { this->computeNodeNeighbourhood(node); }
    ); 
}

/* N log N variant for now */
template<typename Vector, typename Source, std::size_t d>
std::vector<double> BalancedFmmTree<Vector, Source, d>::computePotentials() {
    
    std::vector<double> potentials(5);
    return {}; 
}

/* N log N variant for now */
template<typename Vector, typename Source, std::size_t d>
std::vector<Vector> BalancedFmmTree<Vector, Source, d>::computeForcefields() {
    return {}; 
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::traverseBFSCore(
        const std::function <void(FmmNode*)>& processNode) {

    const std::function <void(BaseNode*)>& processNodeAdaptor = [&processNode]
        (BaseNode* node) { processNode(static_cast<FmmNode*>(node)); }; 

    AOT::traverseBFSCore(processNodeAdaptor); 
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::toFile() {

    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";
    std::string neighbours_filename = "neighbours.dat";
    std::string interaction_filename = "interactions.dat";

    ofstream geometry_file, data_file, neighbours_file, interaction_file; 

    geometry_file.open(geometry_filename);
    data_file.open(data_filename);
    neighbours_file.open(neighbours_filename);
    interaction_file.open(interaction_filename);

    std::size_t n_node = 0;

    traverseBFSCore([&](FmmNode* node) {

        Vector& center = node->center;
        std::size_t depth = node->depth;
        double box_length = node->box_length;

        // 1. log information about the tree itself
        geometry_file << n_node << ", " << depth << ", " << box_length;
        for(auto coord : center.data()) { geometry_file << ", " << coord; }
        geometry_file << std::endl;

        // 2. write to file the sources contained in the various leaves
        if(node->children[0] == nullptr) { 

            FmmLeaf *leaf = static_cast<FmmLeaf*>(node);
            std::vector<Source> &sources = *leaf->sources; 

            data_file << n_node;
            for(const Source &s : sources) {
                for(double coord : s.position.data()) {
                    data_file << ", " << coord;         
                }
                data_file << ", " << s.sourceStrength(); 
            }
            data_file << "\n";
        }

        // 3., 4. for each node, log near neighbour and interaction lists
        std::stringstream head; 
        head << n_node << ", " << depth << ", " << box_length;
        for(auto coord : center.data()) { head << ", " << coord; }

        neighbours_file << head.str(); 
        interaction_file << head.str(); 

        for(auto neighbour : node->near_neighbours) {
            for(auto coord : neighbour->center.data()) {
                neighbours_file << ", " << coord; 
            }
        }
        neighbours_file << std::endl;
        
        for(auto partner : node->interaction_list) {
            for(auto coord : partner->center.data()) {
                interaction_file << ", " << coord; 
            }
        }
        interaction_file << std::endl;

        ++n_node;

    }); 
    
    neighbours_file.close();
    interaction_file.close(); 
}

// Computes near neighbour list and interaction list of a node, assuming that 
// these lists have already been computed for the parent node.
template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::computeNodeNeighbourhood(FmmNode* node) {

    auto& interaction_list = node->interaction_list;
    auto& near_neighbours = node->near_neighbours;

    FmmNode* parent = static_cast<FmmNode*>(node->parent);
    
    // if node is root, node is node's only NN: 
    if(parent == nullptr) { 
        near_neighbours.push_back(node); 
        return; 
    }

    Vector& center = node->center;
    double node_max_neighbour_distance = max_neighbour_distance * node->box_length;

    //If node is not root, divide children of node's parent's NNs (including
    //parent itself) into NNs of node and interaction list of node

    for(auto parent_nn : parent->near_neighbours) {
        for(auto child : parent_nn->children) {
            double distance = (center - child->center).norm();  
            FmmNode* child_ptr = static_cast<FmmNode*>(child);
            if(distance < node_max_neighbour_distance) { // child is near neighbour
                near_neighbours.push_back(child_ptr);  
            }
            else {
                interaction_list.push_back(child_ptr); 
                assert(distance > 1.99 * node->box_length); 
            }
        }
    }
} 

template<typename Vector, typename Source, std::size_t d>
std::tuple<Vector, Vector> BalancedFmmTree<Vector, Source, d>::getDataRange(
        const std::vector<Source> & sources) {

    Vector lower_bounds, upper_bounds; 
    lower_bounds.fill(HUGE_VAL);
    upper_bounds.fill(-HUGE_VAL);

    std::for_each(sources.begin(), sources.end(), [&](Source p) { 
        Vector & pos = p.position; 
        for(std::size_t i = 0; i < d; i++) {
            lower_bounds[i] = pos[i] < lower_bounds[i] ? pos[i] : lower_bounds[i];
            upper_bounds[i] = pos[i] > upper_bounds[i] ? pos[i] : upper_bounds[i];  
        }
    }); 

    // Expand bounding box s.t. all points lie completely within it (and not
    // on boundaries). This simplifies the index calculations in getLeafBoxIndices().
    double paddingFactor = 1E-5;
    lower_bounds = lower_bounds - paddingFactor * lower_bounds;
    upper_bounds = upper_bounds + paddingFactor * upper_bounds;

    return make_tuple(lower_bounds, upper_bounds);
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::expandNode(
        std::queue<FmmNode*>& leaf_queue, FmmNode* node) {

    Vector parent_center = node->center;
    double child_box_length = node->box_length/2;
    std::size_t child_depth = node->depth + 1; 

    if(child_depth < this->height) {
        for(std::size_t j = 0; j < AOT::n_children; ++j) {

            Vector child_center = parent_center + 
                child_box_length/2 * AOT::child_center_directions[j];

            FmmNode* child = new FmmNode(child_center, child_box_length, 
                node->depth + 1, node);

            node->children[j] = child; 
        }
    } 
    else if(child_depth == this->height)  {
        for(std::size_t j = 0; j < AOT::n_children; ++j) {

            Vector child_center = parent_center + 
                child_box_length/2 * AOT::child_center_directions[j];

            FmmNode* child = new FmmLeaf(child_center, child_box_length, 
                node->depth + 1, node);


            leaf_queue.push(child); 
            node->children[j] = child; 
        }
    }
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::distributeSources(
        std::queue<FmmNode*>& leaf_queue, const std::vector<Source>& sources) {
     
    std::size_t n_leaves = leaf_queue.size(); 
    assert(n_leaves == pow(AOT::n_children, this->height));

    std::vector<Source> ** leaf_vectors = new std::vector<Source>*[n_leaves];
    for(std::size_t i = 0; i < n_leaves; i++) {
        leaf_vectors[i] = new std::vector<Source>;
    }

    for(std::size_t k = 0; k < sources.size(); k++) {
        std::array<size_t, d> indices = 
            Super::getLeafBoxIndices(sources[k].position);
        assert(Super::getFlatIndex(indices) < n_leaves); 
        leaf_vectors[Super::getFlatIndex(indices)]->push_back(sources[k]);
    }

    for(std::size_t k = 0; k < n_leaves; k++) {

        FmmLeaf* leaf = static_cast<FmmLeaf*>(leaf_queue.front());
        leaf_queue.pop(); 

        std::array<size_t, d> indices = Super::getLeafBoxIndices(leaf->center);
        leaf->sources = leaf_vectors[Super::getFlatIndex(indices)]; 
    }
     
    delete[] leaf_vectors; 
}

} // namespace fmm

#endif
