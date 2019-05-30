#ifndef FMM_TREE_H
#define FMM_TREE_H

#include <iomanip> 

#include "balanced_orthtree.hpp"
#include "multipole_expansion.hpp"
#include "local_expansion.hpp"
#include "direct.hpp"

namespace fmm {

template<typename Vector, typename Source, std::size_t d>
class BalancedFmmTree: public BalancedOrthtree<Vector, d> {

protected:

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = BalancedOrthtree<Vector, d>;
    using BaseNode = typename AOT::BaseNode; 
    using ME = MultipoleExpansion<Vector, Source, d>;
    using LE = LocalExpansion<Vector, Source, d>;

    static constexpr auto evalScalarInteraction 
        = evaluateInteraction<Vector, Source, double>;
    static constexpr auto evalVectorInteraction 
        = evaluateInteraction<Vector, Source, Vector>;

    static constexpr auto potentialFunction
        = fields::electrostaticPotential<Vector, Source, d>;
    static constexpr auto safePotentialFunction
        = fields::electrostaticPotential_s<Vector, Source, d>;
    static constexpr auto forceFunction
        = fields::electrostaticForce<Vector, Source, d>;
    static constexpr auto safeForceFunction
        = fields::electrostaticForce_s<Vector, Source, d>;

    struct FmmNode; 
    struct FmmLeaf; 

    std::vector<Source>& sources;

    FmmNode* nodes;     
    FmmLeaf* leaves; 

    std::size_t order;
    std::size_t n_nodes;
    std::size_t n_leaves;

    const double max_neighbour_distance = 1.1 * sqrt(d); 
    //in units of box_length + padding to avoid numerical issues

public: 
    BalancedFmmTree(std::vector<Source>& sources, std::size_t items_per_cell, 
            double eps);

    double evaluatePotential(const Vector& eval_point) const; 
    Vector evaluateForcefield(const Vector& eval_point) const; 
    std::vector<double> evaluateParticlePotentials() const; 
    std::vector<Vector> evaluateParticleForcefields() const; 

    void traverseBFSCore(const std::function <void(FmmNode *)>& processNode);
    void toFile() override;

    ~BalancedFmmTree() { 
        delete[] this->nodes;
        delete[] this->leaves;
    }

private:
    static std::tuple<Vector, Vector> getDataRange(
        const std::vector<Source>& sources);
    void expandNode(FmmNode* node, std::size_t& node_offset, 
        std::size_t& leaf_offset);
    void distributeSources();
    void computeNodeNeighbourhood(FmmNode* node); 

    void localToLocal(FmmNode& node);
    void multipoleToLocal(FmmNode& node);
};

template<typename Vector, typename Source, std::size_t d>
struct BalancedFmmTree<Vector, Source, d>::FmmNode: BaseNode {

    using ME = MultipoleExpansion<Vector, Source, d>;
    using LE = LocalExpansion<Vector, Source, d>;

    ME multipole_expansion; 
    LE local_expansion;

    std::vector<FmmNode*> interaction_list;
    std::vector<FmmNode*> near_neighbours;

    FmmNode() {} // Empty default constructor
    FmmNode(Vector center, double box_length, std::size_t depth, 
        FmmNode * parent, std::size_t expansion_order): 
        BaseNode(center, box_length, depth, parent), 
        local_expansion(center, expansion_order) {};

    virtual ~FmmNode() {};
};

template<typename Vector, typename Source, std::size_t d>
struct BalancedFmmTree<Vector, Source, d>::FmmLeaf: 
        BalancedFmmTree<Vector, Source, d>::FmmNode {

    using Super = BalancedFmmTree<Vector, Source, d>::FmmNode;

    std::vector<Source> * sources;

    FmmLeaf() {};
    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent, std::size_t expansion_order): Super(center, 
        box_length, depth, parent, expansion_order), sources() {};

    virtual ~FmmLeaf() { delete sources; }
};

template<typename Vector, typename Source, std::size_t d>
BalancedFmmTree<Vector, Source, d>::BalancedFmmTree(std::vector<Source>& sources, 
        std::size_t sources_per_cell, double eps): BalancedOrthtree<Vector, d>(), 
        sources(sources) { 

    //We compute nearest neighbour & interaction lists by euclidean distance,
    //this breaks down above three dimensions
    static_assert(d <= 3lu); 

    // Determine expansion order, tree height, numbers of nodes and leaves
    order = (ceil(log(1/eps) / log(2))), 
    this->height = ceil(log((double)sources.size()/sources_per_cell) 
        / log(AOT::n_children));

    n_nodes  = pow(AOT::n_children, this->height) / (AOT::n_children - 1);
    n_leaves = pow(AOT::n_children, this->height);

    // Determine bounding box lenghts and center 
    auto [lower_bounds, upper_bounds] = getDataRange(sources); 
    auto extents = (upper_bounds - lower_bounds).data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    std::cout << "Tree has height " << this->height << ", " << n_nodes << 
        " nodes, " << n_leaves << " leaves, order is p = " << order << "\n";

    if(this->height == 0) { // Short circuit construction

        this->nodes = nullptr; 
        this->leaves = new FmmLeaf[n_leaves /* == 1 */];   
        this->leaves[0] = FmmLeaf(center, box_length, 0, nullptr, order);

        std::vector<Source>* leaf_sources = new std::vector<Source>(sources);
        this->leaves[0].sources = leaf_sources;  
        this->root = this->leaves; 

        return;
    }

    this->nodes = new FmmNode[n_nodes];   
    this->leaves = new FmmLeaf[n_leaves];   

    std::size_t node_offset = 0; // Indicates the next free location in this->nodes
    std::size_t leaf_offset = 0; // Indicates the next free location in this->leaves

    this->nodes[node_offset++] = FmmNode(center, box_length, 0, nullptr, order); 
    this->root = this->nodes;  

    // Set up node and leaf structure:
    traverseBFSCore(
        [this, &node_offset, &leaf_offset] (FmmNode * node) 
        { this->expandNode(node, node_offset, leaf_offset); }
    ); 

    distributeSources(); // Assign sources to leaves

    // Compute neighbourhoods of all nodes
    traverseBFSCore( 
        [this](FmmNode* node) -> void { this->computeNodeNeighbourhood(node); }
    ); 

    // Upward pass: Form expansions at leaves and shift to parents

    #pragma omp parallel for
    for(std::size_t i = 0; i < n_leaves; ++i) { 

        FmmLeaf& leaf = this->leaves[i];
        leaf.multipole_expansion = ME(leaf.center, order, *leaf.sources);

    }

    for(int64_t depth = this->height-1; depth >= 0; --depth) {

        std::size_t n_nodes_at_depth = std::pow(AOT::n_children, depth); 
        std::size_t offset = (std::pow(AOT::n_children, depth) - 1) 
            / (AOT::n_children - 1);
         
        // TODO possible if clause to parallelize only if enough grain available
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
    
    /*
    for(int64_t i = n_nodes - 1; i >= 0; --i) { // TODO: parallelizable within levels

        FmmNode& node = this->nodes[i]; 
        std::vector<const ME*> children_expansions;

        for(BaseNode* child : node.children) {
            children_expansions.push_back(
                    &(static_cast<FmmNode*>(child)->multipole_expansion)
            ); 
        }

        node.multipole_expansion = ME(node.center, children_expansions); 

    }
    */
    
    // Downward pass: Convert multipole expansions to local expansions
    for(std::size_t depth = 2; depth < this->height; ++depth) {

        std::size_t n_nodes_at_depth = std::pow(AOT::n_children, depth); 
        std::size_t offset = (std::pow(AOT::n_children, depth) - 1) 
            / (AOT::n_children - 1);

        #pragma omp parallel for
        for(std::size_t i = 0; i < n_nodes_at_depth; i++) {

            FmmNode& current_node = nodes[offset + i]; 

            // Parent LEs and neighbour MEs have to have been built by now
            assert(current_node.multipole_expansion.coefficients.size() > 0); 
            assert(current_node.local_expansion.coefficients.size() > 0);   

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

        // Parent LEs and neighbour MEs have to have been built by now
        assert(current_node.multipole_expansion.coefficients.size() > 0); 
        assert(current_node.local_expansion.coefficients.size() > 0);   

//      auto t1 = std::chrono::high_resolution_clock::now();

        this->localToLocal(current_node);
        this->multipoleToLocal(current_node);

//      auto t2 = std::chrono::high_resolution_clock::now();
//      std::cout << leaf_index << "\t" <<  chrono_duration(t2-t1)  << "\n";
    }
}

template<typename Vector, typename Source, std::size_t d>
double BalancedFmmTree<Vector, Source, d>::
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
            pot += evalScalarInteraction(*(static_cast<FmmLeaf*>(leaf)->sources), 
                eval_point, potentialFunction);
        }
        else {
            pot += evalScalarInteraction(*containing_leaf.sources, 
                    eval_point, safePotentialFunction);
        }
    }
    
    return pot; 
}

template<typename Vector, typename Source, std::size_t d>
Vector BalancedFmmTree<Vector, Source, d>::
        evaluateForcefield(const Vector& eval_point) const {

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
            force_vec += evalVectorInteraction(
                *(static_cast<FmmLeaf*>(leaf)->sources), eval_point, forceFunction);
        }
        else {
            force_vec += evalVectorInteraction(*containing_leaf.sources, 
                    eval_point, safeForceFunction);
        }
    }

    return force_vec;
}

template<typename Vector, typename Source, std::size_t d>
std::vector<double> BalancedFmmTree<Vector, Source, d>::
        evaluateParticlePotentials() const {

    std::vector<double> potentials(sources.size()); 
    for(std::size_t i = 0; i < sources.size(); ++i) {
        potentials[i] = evaluatePotential(sources[i].position);
    }

    return potentials; 
}

template<typename Vector, typename Source, std::size_t d>
std::vector<Vector> BalancedFmmTree<Vector, Source, d>::
        evaluateParticleForcefields() const {

    std::vector<Vector> forces(sources.size()); 
    for(std::size_t i = 0; i < sources.size(); ++i) {
        forces[i] = evaluateForcefield(sources[i].position);
    }
    
    return forces; 
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

    std::ofstream geometry_file, data_file, neighbours_file, interaction_file; 

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

            data_file << n_node 
                << std::setprecision(std::numeric_limits<double>::digits10) 
                << std::scientific;

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
    lower_bounds = lower_bounds + paddingFactor * lower_bounds;
    upper_bounds = upper_bounds + paddingFactor * upper_bounds;

    return std::make_tuple(lower_bounds, upper_bounds);
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::expandNode(FmmNode* node, 
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
                    node->depth+1, node, order);

            node->children[j] = child; 
        }
    }
    else if(child_depth == this->height) { // Child at depth h is leaf

        for(std::size_t j = 0; j < AOT::n_children; ++j) {

            Vector child_center = parent_center + 
                child_box_length/2 * AOT::child_center_directions[j];

            FmmLeaf* child = this->leaves + leaf_offset++; 
            *child = FmmLeaf(child_center, child_box_length, 
                node->depth + 1, node, order);

            node->children[j] = child; 
        }
    }
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::distributeSources() {
     
    std::vector<Source> ** leaf_source_vectors = new std::vector<Source>*[n_leaves];
    for(std::size_t i = 0; i < n_leaves; i++) {
        leaf_source_vectors[i] = new std::vector<Source>;
    }

    // For each source, determine the corresponding leaf. This is equivalent to the
    // first step in bucket sorting
    for(std::size_t k = 0; k < sources.size(); k++) {
        auto indices = Super::getLeafBoxIndices(sources[k].position);
        leaf_source_vectors[this->getFlatIndex(indices)]->push_back(sources[k]);
    }

    // Assign the 'buckets' of source to their respective leaves.     
    // std::size_t source_index = 0; //Next pos. in sources to write to
    for(std::size_t k = 0; k < n_leaves; k++) {

        FmmLeaf& leaf = leaves[k]; 
        auto indices = this->getLeafBoxIndices(leaf.center);
        leaf.sources = leaf_source_vectors[this->getFlatIndex(indices)]; 

        // Store all
        // buckets contiguously in the original sources array. Since the leaves
        // array is in BFS order (or equivalently, Morton sorted), this is
        // equivalent to sorting the buckets w.r.t. the Morton index of their leaves
        //
        // TODO could use this for implicit association of sources and leaves
        // via offsets in sources instead of explicitly storing copies of
        // sources in leaves

//      for(auto source : *leaf.sources) {
//          sources[source_index++] = source; 
//      }
    }
     
    delete[] leaf_source_vectors; 
}

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::localToLocal(FmmNode& node) {
 
    FmmNode * parent = static_cast<FmmNode*>(node.parent);

    assert(parent != nullptr);
    assert(parent->local_expansion.coefficients.size() > 0);   

    node.local_expansion += LE(node.center, parent->local_expansion);
}


template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::multipoleToLocal(FmmNode& node) {
 
    std::vector<const ME*> incoming; 
    for(FmmNode* interaction_partner : node.interaction_list) {
        incoming.push_back(&interaction_partner->multipole_expansion);  
    }

    if(incoming.size() > 0) {
        node.local_expansion += LE(node.center, incoming); 
    }
}

} // namespace fmm

#endif
