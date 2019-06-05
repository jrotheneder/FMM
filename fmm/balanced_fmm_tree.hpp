#ifndef BALANCED_FMM_TREE_H
#define BALANCED_FMM_TREE_H

#include <iomanip> 
#include <filesystem> 

#include "balanced_orthtree.hpp"
#include "abstract_fmm_tree.hpp"
#include "multipole_expansion.hpp"
#include "local_expansion.hpp"
#include "direct.hpp"

namespace fmm {

template<std::size_t d>
class BalancedFmmTree: public BalancedOrthtree<Vector_<d>, d>, 
    public AbstractFmmTree<d> {

protected:

    struct FmmNode; 
    struct FmmLeaf; 

    using Vector = Vector_<d>; 
    using PointSource = PointSource_<d>;

    using AOT = AbstractOrthtree<Vector, d>;
    using AFMMT = AbstractFmmTree<d>;
    using BaseNode = typename AOT::BaseNode; 
    using ME = MultipoleExpansion<d>;
    using LE = LocalExpansion<d>;

    FmmNode* nodes;     
    FmmLeaf* leaves; 

    std::size_t n_nodes;
    std::size_t n_leaves;

public: 
    BalancedFmmTree(std::vector<PointSource>& sources, std::size_t items_per_cell, 
            double eps, double force_smoothing_eps);

    double evaluatePotential(const Vector& eval_point) const; 
    Vector evaluateForcefield(const Vector& eval_point) const; 

    void toFile() override;

    ~BalancedFmmTree() { 
        delete[] this->nodes;
        delete[] this->leaves;
    }

private:
    void expandNode(FmmNode* node, std::size_t& node_offset, 
        std::size_t& leaf_offset);
    void distributePointSources();
    void computeNodeNeighbourhood(FmmNode* node); 
};

template<std::size_t d>
struct BalancedFmmTree<d>::FmmNode: BaseNode {

    MultipoleExpansion<d> multipole_expansion; 
    LocalExpansion<d> local_expansion;

    std::vector<FmmNode*> interaction_list;
    std::vector<FmmNode*> near_neighbours;

    FmmNode() {} // Empty default constructor
    FmmNode(Vector center, double box_length, std::size_t depth, 
        FmmNode * parent, std::size_t expansion_order): 
        BaseNode(center, box_length, depth, parent), 
        local_expansion(center, expansion_order) {};

    virtual ~FmmNode() {};
};

template<std::size_t d>
struct BalancedFmmTree<d>::FmmLeaf: BalancedFmmTree<d>::FmmNode {

    using Super = BalancedFmmTree<d>::FmmNode;

    std::vector<PointSource>* sources;

    FmmLeaf() {};
    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent, std::size_t expansion_order): Super(center, 
        box_length, depth, parent, expansion_order), sources() {};

    virtual ~FmmLeaf() { delete sources; }
};

template<std::size_t d>
BalancedFmmTree<d>::BalancedFmmTree(std::vector<PointSource>& sources, 
        std::size_t sources_per_cell, double eps, double force_smoothing_eps): 
        BalancedOrthtree<Vector, d>(), 
        AbstractFmmTree<d>(sources, force_smoothing_eps) { 

    // Determine expansion order, tree height, numbers of nodes and leaves
    this->order = (ceil(log(1/eps) / log(2))), 
    this->height = ceil(log((double)sources.size()/sources_per_cell) 
        / log(AOT::n_children));

    n_nodes  = pow(AOT::n_children, this->height) / (AOT::n_children - 1);
    n_leaves = pow(AOT::n_children, this->height);

    // Determine bounding box lenghts and center 
    auto [lower_bounds, upper_bounds] = AbstractFmmTree<d>::getDataRange(); 
    auto extents = (upper_bounds - lower_bounds).data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    if(this->height == 0) { // Short circuit construction

        this->nodes = nullptr; 
        this->leaves = new FmmLeaf[n_leaves /* == 1 */];   
        this->leaves[0] = FmmLeaf(center, box_length, 0, nullptr, this->order);
        this->computeNodeNeighbourhood(&this->leaves[0]);

        std::vector<PointSource>* leaf_sources 
            = new std::vector<PointSource>(this->sources);
        this->leaves[0].sources = leaf_sources;  
        this->root = this->leaves; 

        return;
    }

    this->nodes = new FmmNode[n_nodes];   
    this->leaves = new FmmLeaf[n_leaves];   

    std::size_t node_offset = 0; // Indicates the next free location in this->nodes
    std::size_t leaf_offset = 0; // Indicates the next free location in this->leaves

    this->nodes[node_offset++] = FmmNode(center, box_length, 0, nullptr, this->order); 
    this->root = this->nodes;  

    // Set up node and leaf structure:
    this->template traverseBFSCore<FmmNode>(
        [this, &node_offset, &leaf_offset] (FmmNode * node) 
        { this->expandNode(node, node_offset, leaf_offset); }
    ); 

    distributePointSources(); // Assign sources to leaves

    // Compute neighbourhoods of all nodes
    this->template traverseBFSCore<FmmNode>( 
        [this](FmmNode* node) -> void { this->computeNodeNeighbourhood(node); }
    ); 

    // Upward pass: Form expansions at leaves and shift to parents

    #pragma omp parallel for
    for(std::size_t i = 0; i < n_leaves; ++i) { 

        FmmLeaf& leaf = this->leaves[i];
        leaf.multipole_expansion = ME(leaf.center, this->order, *leaf.sources);

    }

    for(unsigned depth = this->height-1; depth > 1; --depth) {

        std::size_t n_nodes_at_depth = std::pow(AOT::n_children, depth); 
        std::size_t offset = (std::pow(AOT::n_children, depth) - 1) 
            / (AOT::n_children - 1);
         
        #pragma omp parallel for 
        for(std::size_t i = 0; i < n_nodes_at_depth; i++) {

            FmmNode& node = nodes[offset + i]; 
            std::vector<const ME*> children_expansions;

            for(BaseNode* child : node.children) {
                children_expansions.push_back(
                        &(static_cast<FmmNode*>(child)->multipole_expansion)
                ); 
            }

            node.multipole_expansion = ME(node.center, children_expansions); 
        }
    }
    
    // Downward pass: Convert multipole expansions to local expansions
    for(std::size_t depth = 2; depth < this->height; ++depth) {

        std::size_t n_nodes_at_depth = std::pow(AOT::n_children, depth); 
        std::size_t offset = (std::pow(AOT::n_children, depth) - 1) 
            / (AOT::n_children - 1);

        #pragma omp parallel for
        for(std::size_t i = 0; i < n_nodes_at_depth; i++) {

            FmmNode& current_node = nodes[offset + i]; 

            //auto t1 = std::chrono::high_resolution_clock::now();

            this->localToLocal(current_node);
            this->multipoleToLocal(current_node);

            //auto t2 = std::chrono::high_resolution_clock::now();
            //std::cout << offset + i << "\t" <<  chrono_duration(t2-t1)  << "\n";

        }
    }

    #pragma omp parallel for 
    for(std::size_t leaf_index = 0; leaf_index < n_leaves; ++leaf_index) {

        FmmNode& current_node = leaves[leaf_index]; 

//      auto t1 = std::chrono::high_resolution_clock::now();

        this->localToLocal(current_node);
        this->multipoleToLocal(current_node);

//      auto t2 = std::chrono::high_resolution_clock::now();
//      std::cout << leaf_index << "\t" <<  chrono_duration(t2-t1)  << "\n";
    }
}

template<std::size_t d>
double BalancedFmmTree<d>::
        evaluatePotential(const Vector& eval_point) const {

    uint64_t leaf_index = this->getMortonIndex(this->getLeafBoxIndices(eval_point)); 
    FmmLeaf& containing_leaf = leaves[leaf_index];  

    // Contributions from local expansion of containing leaf
    double pot = containing_leaf.local_expansion.evaluatePotential(eval_point); 

    // Contributions from sources in near neighbour leaves (eval. directly).
    // The list of NNs includes the containing_leaf; this one is handled
    // separately with a potential which checks whether eval_point coincides
    // with a source location (in which case an exception is thrown)
    for(auto leaf : containing_leaf.near_neighbours) {
        if(leaf != &containing_leaf) {
            pot += AFMMT::evalScalarInteraction(*(static_cast<FmmLeaf*>(leaf)->sources), 
                eval_point, this->force_smoothing_eps, AFMMT::potentialFunction);
        }
        else {
            pot += AFMMT::evalScalarInteraction(*containing_leaf.sources, 
                eval_point, this->force_smoothing_eps, AFMMT::safePotentialFunction);
        }
    }
    
    return pot; 
}

template<std::size_t d>
Vector_<d> BalancedFmmTree<d>::evaluateForcefield(const Vector_<d>& eval_point) 
        const {

    uint64_t leaf_index = this->getMortonIndex(this->getLeafBoxIndices(eval_point)); 
    FmmLeaf& containing_leaf = leaves[leaf_index];  

    // Contributions from local expansion of containing leaf
    Vector force_vec = containing_leaf.local_expansion.evaluateForcefield(eval_point); 

    // Contributions from sources in near neighbour leaves (eval. directly).
    // The list of NNs includes the containing_leaf; this one is handled
    // separately with a force which checks whether eval_point coincides
    // with a source location (in which case an exception is thrown) 
    for(auto leaf : containing_leaf.near_neighbours) {
        if(leaf != &containing_leaf) {
            force_vec += AFMMT::evalVectorInteraction(
                *(static_cast<FmmLeaf*>(leaf)->sources), eval_point, 
                this->force_smoothing_eps, AFMMT::forceFunction);
        }
        else {
            force_vec += AFMMT::evalVectorInteraction(*containing_leaf.sources, 
                eval_point, this->force_smoothing_eps, AFMMT::safeForceFunction);
        }
    }

    return force_vec;
}

template<std::size_t d>
void BalancedFmmTree<d>::toFile() {

    namespace fs =  std::filesystem;
    std::string logs_dir = "./logs"; 
    fs::create_directory(fs::path(logs_dir)); 
    
    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";
    std::string neighbours_filename = "neighbour_lists.dat";
    std::string interaction_filename = "interaction_lists.dat";

    std::ofstream geometry_file, data_file, neighbours_file, interaction_file; 

    geometry_file.open(logs_dir + "/" + geometry_filename);
    data_file.open(logs_dir + "/" + data_filename);
    neighbours_file.open(logs_dir + "/" + neighbours_filename);
    interaction_file.open(logs_dir + "/" + interaction_filename);

    std::size_t n_node = 0;

    this->template traverseBFSCore<FmmNode>([&](FmmNode* node) {

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
            std::vector<PointSource> &sources = *leaf->sources; 

            data_file << n_node 
                << std::setprecision(std::numeric_limits<double>::digits10) 
                << std::scientific;

            for(const PointSource &s : sources) {
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
}

// Computes near neighbour list and interaction list of a node, assuming that 
// these lists have already been computed for the parent node.
template<std::size_t d>
void BalancedFmmTree<d>::computeNodeNeighbourhood(FmmNode* node) {

    auto& interaction_list = node->interaction_list;
    auto& near_neighbours = node->near_neighbours;

    FmmNode* parent = static_cast<FmmNode*>(node->parent);
    
    // if node is root, node is node's only NN: 
    if(parent == nullptr) { 
        near_neighbours.push_back(node); 
        return; 
    }

    Vector& center = node->center;

    //If node is not root, divide children of node's parent's NNs (including
    //parent itself) into NNs of node and interaction list of node

    for(auto parent_nn : parent->near_neighbours) {
        for(auto child : parent_nn->children) {
            double distance = (center - child->center).norm();  
            FmmNode* child_ptr = static_cast<FmmNode*>(child);
            if(node->adjacent(child_ptr)) { // child is near neighbour
            //if(distance < this->max_neighbour_distance * node->box_length) {
            //child is near neighbour
                near_neighbours.push_back(child_ptr);  
            }
            else {
                interaction_list.push_back(child_ptr); 

                if(distance < 1.99 * node->box_length) {
//                  std::cerr << "Distance is " << distance 
//                      << ", vs 1.99 * node->box_length = " 
//                      << 1.99 * node->box_length << "\n";
//                  std::cerr << "Depth is " << node->depth << "\n";
//                  std::cerr << node->center << "\t" << node->box_length << "\n";
//                  std::cerr << child->center << "\t" << child->box_length << "\n";
                    throw std::logic_error("Node too close to be part of the "
                            "interaction list.");
                } 
            }
        }
    }
} 

template<std::size_t d>
void BalancedFmmTree<d>::expandNode(FmmNode* node, 
        std::size_t& node_offset, std::size_t& leaf_offset) {

    Vector parent_center = node->center;
    double child_box_length = node->box_length/2;
    std::size_t child_depth = node->depth + 1; 

    if(child_depth < this->height) { // Child at depth < h is node
        for(std::size_t j = 0; j < AOT::n_children; ++j) {

            Vector child_center = parent_center + 
                child_box_length/2 * AOT::child_center_directions[j];

            FmmNode* child = this->nodes + node_offset++; 
            *child = FmmNode(child_center, child_box_length, 
                    node->depth+1, node, this->order);

            node->children[j] = child; 
        }
    }
    else if(child_depth == this->height) { // Child at depth h is leaf

        for(std::size_t j = 0; j < AOT::n_children; ++j) {

            Vector child_center = parent_center + 
                child_box_length/2 * AOT::child_center_directions[j];

            FmmLeaf* child = this->leaves + leaf_offset++; 
            *child = FmmLeaf(child_center, child_box_length, 
                node->depth + 1, node, this->order);

            node->children[j] = child; 
        }
    }
}

template<std::size_t d>
void BalancedFmmTree<d>::distributePointSources() {
     
    std::vector<PointSource> ** leaf_source_vectors = new std::vector<PointSource>*[n_leaves];
    for(std::size_t i = 0; i < n_leaves; i++) {
        leaf_source_vectors[i] = new std::vector<PointSource>;
    }

    // For each source, determine the corresponding leaf. This is equivalent to the
    // first step in bucket sorting
    for(std::size_t k = 0; k < this->sources.size(); k++) {
        auto indices = BalancedOrthtree<Vector, d>::getLeafBoxIndices(
            this->sources[k].position);
        leaf_source_vectors[this->getFlatIndex(indices)]->push_back(this->sources[k]);
    }

    // Assign the 'buckets' of source to their respective leaves.     
    // std::size_t source_index = 0; //Next pos. in this->sources to write to
    for(std::size_t k = 0; k < n_leaves; k++) {

        FmmLeaf& leaf = leaves[k]; 
        auto indices = this->getLeafBoxIndices(leaf.center);
        leaf.sources = leaf_source_vectors[this->getFlatIndex(indices)]; 
    }
     
    delete[] leaf_source_vectors; 
}

} // namespace fmm

#endif
