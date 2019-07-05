#ifndef ADAPTIVE_FMM_TREE_H
#define ADAPTIVE_FMM_TREE_H

#include <iomanip> 
#include <filesystem> 

#include "adaptive_orthtree.hpp"
#include "abstract_fmm_tree.hpp"

namespace fmm {

template<std::size_t d, bool field_type = true> // true => Grav., false => Coulomb
class AdaptiveFmmTree: public AdaptiveOrthtree<Vector_<d>, d>, 
    public AbstractFmmTree<d, field_type> {

protected:

    struct FmmNode; 
    struct FmmLeaf; 

    using Vector = Vector_<d>; 
    using PointSource = PointSource_<d>;

    using AOT = AbstractOrthtree<Vector, d>;
    using ADOT = AdaptiveOrthtree<Vector, d>;
    using AFMMT = AbstractFmmTree<d, field_type>;
    using BaseNode = typename AOT::BaseNode; 
    using ME = MultipoleExpansion<d>;
    using LE = LocalExpansion<d>;

    std::vector<std::vector<FmmNode*>> nodes;  // store nodes and leaves levelwise
    std::vector<std::vector<FmmLeaf*>> leaves; 

    std::size_t max_sources_per_leaf;

public: 
    AdaptiveFmmTree(std::vector<PointSource>& sources, 
        std::size_t max_sources_per_leaf, double eps, 
        double force_smoothing_eps = 0);

    double evaluatePotential(const Vector& eval_point) const override; 
    Vector evaluateForcefield(const Vector& eval_point) const override; 

    void toFile() override;

    ~AdaptiveFmmTree() { delete this->root; }

private:
    static std::array<std::vector<PointSource>, AOT::n_children> refineOrthant(
        const Vector& center, const std::vector<PointSource>& sources);
    void buildTree(FmmNode* node, std::vector<PointSource>& sources, 
            std::size_t max_sources_per_leaf);
    void computeNodeNeighbourhood(FmmNode* node); 

    void sourceToLocal(FmmNode& node);
    void constructLocalExpansion(FmmNode& node); 

    FmmLeaf* getContainingLeaf(const Vector& point) const; 

};

template<std::size_t d, bool field_type>
struct AdaptiveFmmTree<d, field_type>::FmmNode: BaseNode {

    MultipoleExpansion<d> multipole_expansion; 
    LocalExpansion<d> local_expansion;

    std::vector<FmmNode*> near_neighbours;  // list U for leaves, all NNs for nodes
                                            // (includes neighbours at all levels)
    std::vector<FmmNode*> interaction_list; // list V 
    std::vector<FmmLeaf*> X_list;           // list X 

    FmmNode() {} // Empty default constructor
    FmmNode(Vector center, double box_length, std::size_t depth, 
        FmmNode * parent, std::size_t expansion_order): 
        BaseNode(center, box_length, depth, parent), 
        local_expansion(center, expansion_order) {};

    virtual bool isLeaf() { return false; }

    virtual ~FmmNode() { for(auto child : this->children) { delete child; } }
};

template<std::size_t d, bool field_type>
struct AdaptiveFmmTree<d, field_type>::FmmLeaf: 
        AdaptiveFmmTree<d, field_type>::FmmNode {

    using Super = AdaptiveFmmTree<d, field_type>::FmmNode;

    std::vector<FmmNode*> W_list; // list W 
    std::vector<PointSource> sources;

    FmmLeaf() {};
    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent, std::size_t expansion_order): Super(center, 
        box_length, depth, parent, expansion_order), sources() {};

    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent, std::size_t expansion_order, 
        std::vector<PointSource> sources): Super(center, box_length, depth, 
        parent, expansion_order), sources(sources) {};

    virtual bool isLeaf() override { return true; }

    virtual ~FmmLeaf() {}
};

template<std::size_t d, bool field_type>
AdaptiveFmmTree<d, field_type>::AdaptiveFmmTree(std::vector<PointSource>& sources, 
        std::size_t max_sources_per_leaf, double eps, double force_smoothing_eps): 
        AdaptiveOrthtree<Vector, d>(), 
        AbstractFmmTree<d, field_type>(sources, eps, force_smoothing_eps) { 

    // Determine bounding box lenghts and center 
    auto [lower_bounds, upper_bounds] 
        = AbstractFmmTree<d, field_type>::getDataRange(); 
    auto extents = (upper_bounds - lower_bounds).data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    if(sources.size() <= max_sources_per_leaf) { // Short circuit construction

        FmmLeaf * root = new FmmLeaf(center, box_length, 0, nullptr,
            this->order, sources);
        this->computeNodeNeighbourhood(root);
        this->leaves.push_back(std::vector<FmmLeaf*>{root});

        this->root = root; 
        this->height = 0;

        return;
    }

    FmmNode * root = new FmmNode(center, box_length, 0, nullptr, this->order);;
    this->root = root; 

    this->nodes.push_back(std::vector<FmmNode*>{root});
    this->leaves.resize(1); 

    // TODO: consider discarding empty leaves entirely (replace with nullptr)
    // for optimization 
    buildTree(this->nodes[0][0], sources, max_sources_per_leaf); 

    // Compute neighbourhoods of all nodes
    this->template traverseBFSCore<FmmNode>( 
        [this](FmmNode* node) -> void { this->computeNodeNeighbourhood(node); }
    ); 

    // Upward pass: Form expansions at leaves and shift to parents

    // TODO: all leaves can be processed in parallel here, not just those on 
    // the same level. It may be advantageous to concatenate all the leaf levels 
    // into one big iterator (boost has something for that) s.t. we have only 
    // one loop to parallelize.
    
    for(std::vector<FmmLeaf*>& level : this->leaves) { 
        #pragma omp parallel for 
        for(unsigned j = 0; j < level.size(); ++j) {
            FmmLeaf* leaf = level[j];  
            leaf->multipole_expansion = ME(leaf->center, this->order, leaf->sources);
        }
    }

    assert(this->height == this->nodes.size()); 
    assert(this->height+1 == this->leaves.size()); 

    for(unsigned i = this->height-1; i > 1; --i) { 

        std::vector<FmmNode*>& node_level = this->nodes[i]; 

        #pragma omp parallel for 
        for(unsigned j = 0; j < node_level.size(); ++j) {

            FmmNode* node = node_level[j]; 
            std::vector<const ME*> children_expansions;

            for(BaseNode* child : node->children) {
                children_expansions.push_back(
                    &(static_cast<FmmNode*>(child)->multipole_expansion)
                ); 
            }

            node->multipole_expansion = ME(node->center, children_expansions); 
        }
    }
    
    // Downward pass: Convert multipole expansions to local expansions, shift 
    // parent local expansions to children, construct local expansions from
    // sources contained in leaves contained in the X list
    for(std::size_t depth = 2; depth < this->height; ++depth) {

        std::vector<FmmNode*>& node_level = this->nodes[depth];
        std::vector<FmmLeaf*>& leaf_level = this->leaves[depth];

        #pragma omp parallel 
        {
            #pragma omp for nowait schedule(dynamic)
            for(std::size_t i = 0; i < node_level.size(); i++) {

                FmmNode& node = *node_level[i]; 
                constructLocalExpansion(node); 

            }
            #pragma omp for schedule(dynamic)
            for(std::size_t i = 0; i < leaf_level.size(); i++) {

                FmmNode& node = *leaf_level[i]; 
                constructLocalExpansion(node); 
            }
        }
    }

    std::vector<FmmLeaf*>& final_level = this->leaves.at(this->height); 
    #pragma omp parallel for schedule(dynamic)
    for(std::size_t i = 0; i < final_level.size(); ++i) {

        FmmNode& node = *final_level[i]; 
        constructLocalExpansion(node); 
    }
}

template<std::size_t d, bool field_type>
double AdaptiveFmmTree<d, field_type>::evaluatePotential(const Vector& eval_point) 
        const {

    // use the grav. field here and optionally convert to Coulomb at return
    const bool grav = true; 
    const bool safe = true;

    FmmLeaf* containing_leaf = getContainingLeaf(eval_point);  
    if(containing_leaf == nullptr) { // outside of tree => must compute directly
        return fields::potential<d, field_type, !safe>(this->sources, 
            eval_point);  
    }
    
    // Contributions from local expansion of containing leaf
    double pot = containing_leaf->local_expansion.evaluatePotential(eval_point); 

    // Contributions from multipole expansions of nodes & leaves contained in
    // containing leaf's W list

    for(FmmNode* node : containing_leaf->W_list) {
        pot += node->multipole_expansion.evaluatePotential(eval_point);  
    }

    // Contributions from sources in near neighbour leaves (eval. directly).
    // The list of NNs includes the containing_leaf; this one is handled
    // separately with a potential which checks whether eval_point coincides
    // with a source location, in which case that source is ignored 

    for(auto leaf : containing_leaf->near_neighbours) {
        if(leaf != containing_leaf) {
            pot += fields::potential<d, grav, !safe>(
                static_cast<FmmLeaf*>(leaf)->sources, 
                eval_point);
        }
        else {
            pot += fields::potential<d, grav, safe>(
                static_cast<FmmLeaf*>(leaf)->sources, 
                eval_point);
        }
    }
    
    if constexpr(field_type == 0) { // Coulomb field => switch signs
        pot = -pot; 
    }

    return pot; 
}

template<std::size_t d, bool field_type>
Vector_<d> AdaptiveFmmTree<d, field_type>::evaluateForcefield(
        const Vector_<d>& eval_point) const {

    // we use the grav. field here and optionally convert to Coulomb at return
    const bool grav = 1; 
    const bool safe = 1;

    FmmLeaf* containing_leaf = getContainingLeaf(eval_point);  
    if(containing_leaf == nullptr) { // outside of tree => must compute directly
        return fields::forcefield<d, field_type, !safe>(this->sources, 
            eval_point, this->force_smoothing_eps);  
    }

    // Contributions from local expansion of containing leaf
    Vector force_vec = containing_leaf->local_expansion.evaluateForcefield(eval_point); 

    // Contributions from multipole expansions of nodes & leaves contained in
    // containing leaf's W list
    for(FmmNode* node : containing_leaf->W_list) {
        force_vec += node->multipole_expansion.evaluateForcefield(eval_point);  
    }

    // Contributions from sources in near neighbour leaves (eval. directly).
    // The list of NNs includes the containing_leaf; this one is handled
    // separately with a force which checks whether eval_point coincides
    // with a source location, in which case that source is ignored 


    for(auto leaf : containing_leaf->near_neighbours) {
        if(leaf != containing_leaf) {
            force_vec += fields::forcefield<d, grav, !safe>(
                static_cast<FmmLeaf*>(leaf)->sources, 
                eval_point, this->force_smoothing_eps);
        }
        else {
            force_vec += fields::forcefield<d, grav, safe>(
                static_cast<FmmLeaf*>(leaf)->sources, 
                eval_point, this->force_smoothing_eps);
        }
    }

    if constexpr(field_type == 0) { // Coulomb field => switch signs
        force_vec = -force_vec; 
    }

    return force_vec;
}

template<std::size_t d, bool field_type> 
void AdaptiveFmmTree<d, field_type>::toFile() {

    namespace fs =  std::filesystem;
    std::string logs_dir = "./logs"; 
    fs::create_directory(fs::path(logs_dir)); 
    
    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";
    std::string neighbours_filename = "neighbour_lists.dat";
    std::string interaction_filename = "interaction_lists.dat";
    std::string W_list_filename = "W_lists.dat";
    std::string X_list_filename = "X_lists.dat";

    std::ofstream geometry_file, data_file, neighbours_file, 
        interaction_file, W_file, X_file; 

    geometry_file.open(logs_dir + "/" + geometry_filename);
    data_file.open(logs_dir + "/" + data_filename);
    neighbours_file.open(logs_dir + "/" + neighbours_filename);
    interaction_file.open(logs_dir + "/" + interaction_filename);
    W_file.open(logs_dir + "/" + W_list_filename);  
    X_file.open(logs_dir + "/" + X_list_filename);  

    std::size_t n_node = 0;

    this->template traverseBFSCore<FmmNode>([&](FmmNode* node) {

        Vector& center = node->center;
        std::size_t depth = node->depth;
        double box_length = node->box_length;

        std::stringstream head; 
        head << n_node << ", " << node->isLeaf() << ", " << depth 
            << ", " << box_length;

        // 1. log information about the structure of the tree
        geometry_file << head.str();
        for(auto coord : center.data()) { geometry_file << ", " << coord; }
        geometry_file << std::endl;

        // 2., 3. if node is a leaf, write to file the sources contained in node
        // log the associated W list
        if(node->isLeaf()) { 

            FmmLeaf *leaf = static_cast<FmmLeaf*>(node);
            std::vector<PointSource>& sources = leaf->sources; 

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

            W_file << head.str(); 
            for(auto partner : leaf->W_list) {
                for(auto coord : partner->center.data()) {
                    W_file << ", " << coord; 
                }
            }
            W_file << std::endl;
        }

        // 4., 5., 6. log near neighbour, interaction and X lists,
        for(auto coord : center.data()) { head << ", " << coord; }

        neighbours_file << head.str(); 
        for(auto neighbour : node->near_neighbours) {
            for(auto coord : neighbour->center.data()) {
                neighbours_file << ", " << coord; 
            }
        }
        neighbours_file << std::endl;
        
        interaction_file << head.str(); 
        for(auto partner : node->interaction_list) {
            for(auto coord : partner->center.data()) {
                interaction_file << ", " << coord; 
            }
        }
        interaction_file << std::endl;

        X_file << head.str(); 
        for(auto partner : node->X_list) {
            for(auto coord : partner->center.data()) {
                X_file << ", " << coord; 
            }
        }
        X_file << std::endl;

        ++n_node;

    }); 
}

template<std::size_t d, bool field_type>
std::array<std::vector<PointSource_<d>>, AbstractOrthtree<Vector_<d>, d>::
        n_children> AdaptiveFmmTree<d, field_type>::refineOrthant(const Vector_<d>& 
        center, const std::vector<PointSource_<d>>& sources) {

    std::array<std::vector<PointSource>, AOT::n_children> child_sources; 

    for(auto& source : sources) {
        unsigned index = ADOT::getOrthant(center, source.position); 
        child_sources[index].push_back(source);   
    }
    
    return child_sources; 
}

template<std::size_t d, bool field_type>
void AdaptiveFmmTree<d, field_type>::buildTree(FmmNode* node, 
        std::vector<PointSource_<d>>& sources, std::size_t max_items_per_leaf) {

    unsigned depth = 0; 
    bool refine = true; 

    //sources of nodes at depth which still need to be assigned to 
    //leaves further down in the tree. Unfortunately we incur quite a few copies
    //of large vectors here, but storage requirements are limited to 2N + overhead, 
    //and we are cpu bound anyways 
    std::vector<std::vector<PointSource>> node_sources{sources}; 
    // temporary storage for sources assigned to the nodes on the next level
    std::vector<std::vector<PointSource>> temp_node_sources{}; 

    while(refine) { // construct the tree level by level

        ++depth;

        std::vector<FmmNode*> next_node_level; 
        std::vector<FmmLeaf*> next_leaf_level; 

        for(unsigned i = 0; i < nodes[depth-1].size(); ++i) {

            FmmNode* node = nodes[depth-1][i];

            Vector parent_center = node->center; 
            double child_box_length = node->box_length/2; 
            unsigned child_depth = node->depth + 1;

            assert(child_depth == depth);  

            auto child_sources = refineOrthant(node->center, node_sources[i]); 

            for(std::size_t i = 0; i < AOT::n_children; ++i) {

                Vector child_center = parent_center + 
                    child_box_length/2 * AOT::child_center_directions[i];

                if(child_sources[i].size() > max_items_per_leaf) { // refine

                    FmmNode* child = new FmmNode(child_center, child_box_length, 
                                child_depth, node, this->order);
                    node->children[i] = child;

                    next_node_level.push_back(child);
                    temp_node_sources.push_back(child_sources[i]);

                }
                else {
                    FmmLeaf* child = new FmmLeaf(child_center, child_box_length, 
                        child_depth, node, this->order, child_sources[i]);
                    node->children[i] = child;

                    next_leaf_level.push_back(child);

                }
            }
        }

        node_sources.swap(temp_node_sources); 
        temp_node_sources.clear(); 
        
        refine = !next_node_level.empty(); 
        if(refine) { nodes.push_back(next_node_level); }
        leaves.push_back(next_leaf_level); 
    }

    this->height = depth; 
}

// Computes near neighbour list and interaction list of a node, assuming that 
// these lists have already been computed for the parent node.
template<std::size_t d, bool field_type>
void AdaptiveFmmTree<d, field_type>::computeNodeNeighbourhood(FmmNode* node) {

    FmmNode* parent = static_cast<FmmNode*>(node->parent);
    Vector& center = node->center;
    
    // if node is root, node is node's only NN: 
    if(parent == nullptr) { 
        node->near_neighbours.push_back(node); 
        return; 
    }

    for(auto parent_neighbour : parent->near_neighbours) {

        if(parent_neighbour->isLeaf()) {
            // parent neighbour is leaf: if adjacent, add this leaf to node's NN
            // list, but also add node to that leaf's NN list. 
            if(node->adjacent(parent_neighbour)) {
                node->near_neighbours.push_back(parent_neighbour); 
                if(node->isLeaf()) { parent_neighbour->near_neighbours.push_back(node);
                }
            } 
            // non adjacent leaf => node is descendant of a neighbour of that
            // leaf but not itself adjacent to the leaf and hence part of this leaf's 
            // W list; likewise the leaf itself is part of node's X list
            else {  
                FmmLeaf* leaf = static_cast<FmmLeaf*>(parent_neighbour);
                node->X_list.push_back(leaf);
                leaf->W_list.push_back(node); 
            }
        }
        else { // the current neighbour is a node, so we look at its children
            for(auto child : parent_neighbour->children) {

                FmmNode* child_ptr = static_cast<FmmNode*>(child);

                if(node->adjacent(child)) { // child is near neighbour
                    // If node is a leaf, child gets added to its NN List only
                    // if child itself is a leaf (else some of the leaf descendants of
                    // child will get added to node's NN list later on). On the
                    // other hand if node is not a leaf, all adjacent nodes
                    // get added to the NN list in this step.
                    if(!node->isLeaf() || child_ptr->isLeaf()) {
                        node->near_neighbours.push_back(child_ptr);  
                    }
                }
                else {
                    node->interaction_list.push_back(child_ptr); 

                    double distance = (center - child->center).norm();  
                    if(distance < 1.99 * node->box_length) {
                        throw std::logic_error("Node too close to be part of the "
                            "interaction list.");
                    }
                }
            }
        }
    }
} 

template<std::size_t d, bool field_type>
void AdaptiveFmmTree<d, field_type>::sourceToLocal(FmmNode& node) {

    for(FmmLeaf* leaf : node.X_list) {
        node.local_expansion += LE(node.center, this->order, leaf->sources);          
    }
}

template<std::size_t d, bool field_type>
void AdaptiveFmmTree<d, field_type>::constructLocalExpansion(FmmNode& node) {

    this->localToLocal(node);
    this->multipoleToLocal(node);
    this->sourceToLocal(node); 
}

template<std::size_t d, bool field_type>
typename AdaptiveFmmTree<d, field_type>::FmmLeaf* AdaptiveFmmTree<d, field_type>::
        getContainingLeaf(const Vector& point) const {

    FmmNode* current = static_cast<FmmNode*>(this->root); 

    while(!current->isLeaf()) {
        unsigned child_index = this->getOrthant(current->center, point); 
        current = static_cast<FmmNode*>(current->children[child_index]);   
    }

    if(current->contains(point)) { return static_cast<FmmLeaf*>(current); }
    else { return nullptr; }
}


} // namespace fmm

#endif
