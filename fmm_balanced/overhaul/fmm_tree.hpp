#ifndef FMM_TREE_H
#define FMM_TREE_H

#include "balanced_orthtree.hpp"
#include "fmm_node.hpp"

namespace fmm {

template<typename Vector, typename Source, std::size_t d>
class BalancedFmmTree: public BalancedOrthtree<Vector, d> {

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = BalancedOrthtree<Vector, d>;

    using BaseNode = typename AOT::Node; 
    using FmmNode = typename node::FmmNode<Vector, Source, d>; 
    using FmmLeaf = typename node::FmmLeaf<Vector, Source, d>; 

public: 

    static constexpr double max_neighbour_distance = 1.1 * sqrt(d); 
    //in units of box_length + padding to avoid numerical issues

    BalancedFmmTree(std::vector<Source> sources, std::size_t items_per_cell, 
            double eps);

    void traverseBFSCore(const std::function <void(FmmNode *)>& processNode);
    void toFile() override;

    static void computeNodeNeighbourhood(FmmNode* node); 
    static std::tuple<Vector, Vector> getDataRange(const std::vector<Source> & sources);
};

template<typename Vector, typename Source, std::size_t d>
BalancedFmmTree<Vector, Source, d>::BalancedFmmTree(std::vector<Source> sources, 
        std::size_t items_per_cell, double eps): BalancedOrthtree<Vector, d>() { 

    if(d > 3lu) { 
        throw std::runtime_error("Determining near neighbours from "
            "euclidean center separation is only valid for d <= 3 "
            "dimensions, modify implementation."); 
    }

    // Determine tree height, bounding box lenghts and center 
    this->height = ceil(log((double)sources.size()/items_per_cell) / 
            log(AOT::n_children));

    Vector lower_bounds, upper_bounds;  
    std::tie(lower_bounds, upper_bounds) = getDataRange(sources); 
    auto extents = (upper_bounds - lower_bounds).data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());
    Vector center = 0.5 * (lower_bounds + upper_bounds); 

    // Build tree: 
    this->root = new FmmNode(center, box_length, 0, nullptr); 

//  std::queue<FmmNode*> node_queue({static_cast<FmmNode*>(this->root)}); 
    std::queue<FmmNode*> leaf_queue; 

    traverseBFSCore([this, &leaf_queue](FmmNode * node) {

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
    }); 

/*  while(child_depth < this->height) {

        std::size_t n_nodes = node_queue.size(); 
        assert(n_nodes == pow(AOT::n_children, child_depth));
        
        ++child_depth;

        for(std::size_t i = 0; i < n_nodes; i++) {

            FmmNode* parent = node_queue.front();
            node_queue.pop(); 

            Vector parent_center = parent->center;
            double child_box_length = parent->box_length/2;

            for(int j = 0; j < AOT::n_children; j++) {

                FmmNode* child;
                Vector child_center = parent_center + 
                    child_box_length/2 * AOT::child_center_directions[j];

                if(child_depth < this->height) {
                    child = new FmmNode(child_center, child_box_length, 
                                        child_depth, parent); 
                }
                else { // Treat leaves seperately
                    child = new FmmLeaf(child_center, child_box_length, 
                                        child_depth, parent); 
                }

                node_queue.push(child);
                parent->children[j] = child; 
            }
        }
    }
*/
    // Distribute sources to leaves
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



/* 
template<typename Vector, typename Source, std::size_t d>
BalancedFmmTree<Vector, Source, d>::BalancedFmmTree(std::vector<Vector> points, 
        std::size_t s, double eps): Super(points, s) { 

    // Theoretically we could skip the call to Super() here and construct
    // everything in one. TODO: consider this NOW

    if(d > 3lu) { 
        throw std::runtime_error("Determining near neighbours from "
            "euclidean center seperation is only valid for d <= 3 "
            "dimensions, modify implementation."); 
    }

    // Maximal distance to near neighbours, in units of box_length and with 
    // padding to avoid numerical issues
    const double max_neighbour_distance = 1.1 * sqrt(d); 

    Super::traverseBFSCore([max_neighbour_distance](typename Super::Node * node) {

        FmmData* fmm_data = new FmmData();
        node->data = fmm_data; 

        auto& interaction_list = fmm_data->interaction_list;
        auto& near_neighbours = fmm_data->near_neighbours;

        Node* parent = node->parent;
        Vector& center = node->center;

        double box_max_neighbour_distance = max_neighbour_distance 
            * node->box_length;

        //Node is not root => sort children of parent's NNs (including
        //parent itself) into NNs of node and interaction list of node
        if(parent) { 
            auto& parent_near_neighbours = 
                static_cast<FmmData*>(parent->data)->near_neighbours; 

            for(auto parent_nn : parent_near_neighbours) {
                for(auto child : parent_nn->children) {
                    double distance = (center - child->center).norm();  
                    if(distance < box_max_neighbour_distance) { // near neighbour
                        near_neighbours.push_back(child);  
                    }
                    else {
                        interaction_list.push_back(child); 
                        assert(distance > 1.99 * node->box_length); 
                    }
                }
            }
        } 
        else { near_neighbours.push_back(node); }// Root is root's only NN 

    }); 
}
*/

template<typename Vector, typename Source, std::size_t d>
void BalancedFmmTree<Vector, Source, d>::traverseBFSCore(
        const std::function <void(FmmNode *)>& processNode) {

    const std::function <void(BaseNode *)>& processNodeAdaptor = 
        [&processNode](BaseNode* node) { 
            processNode(static_cast<FmmNode*>(node)); 
        }; 

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

    // TODO: log source, geometric information
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
    double node_max_neighbour_distance = max_neighbour_distance 
        * node->box_length;

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
    // on boundaries. This simplifies the index calculations in getLeafBoxIndices().

    double paddingFactor = 1E-5;
    lower_bounds = lower_bounds - paddingFactor * lower_bounds;
    upper_bounds = upper_bounds + paddingFactor * upper_bounds;

    return make_tuple(lower_bounds, upper_bounds);
}

} // namespace fmm

#endif
