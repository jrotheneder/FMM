#ifndef POINT_ORTHTREE_H
#define POINT_ORTHTREE_H

#include "balanced_orthtree.hpp"

template<typename Vector, std::size_t d>
class PointOrthtree: public BalancedOrthtree<Vector, d> {

    using AOT = AbstractOrthtree<Vector, d>;
    using Super = BalancedOrthtree<Vector, d>;
    using Node = typename AOT::Node; 
    using AbstractOrthtree<Vector, d>::height;

public:
    struct Leaf; 

    PointOrthtree(std::vector<Vector> points, std::size_t s);

    void toFile() override; 

    ~PointOrthtree() { delete this->root; }

};

template<typename Vector, std::size_t d>
struct PointOrthtree<Vector, d>::Leaf: AbstractOrthtree<Vector, d>::Node {

    std::vector<Vector> * points;
    Leaf(Vector center, double box_length, std::size_t depth, 
        Node * parent): Node(center, box_length, depth, parent), points() {}
    virtual ~Leaf() { delete points; }
};

template<typename Vector, std::size_t d>
PointOrthtree<Vector, d>::PointOrthtree(std::vector<Vector> points, 
        std::size_t items_per_leaf): BalancedOrthtree<Vector, d>() {

    // Determine tree height, bounding box lenghts and center as well as
    // the directions in which the child centers lie relative to the center
    // of a box
    this->height = ceil(log((double)points.size()/items_per_leaf) / 
            log(AOT::n_children));

    Vector lower_bounds, upper_bounds;  
    std::tie(lower_bounds, upper_bounds) = AOT::getDataRange(points);
    auto extents = Vector{upper_bounds - lower_bounds}.data(); 
    double box_length = *std::max_element(extents.begin(), extents.end());

    Vector center = 0.5 * (lower_bounds + upper_bounds); 

//  std::cout << "Orthtree height is " << height << endl;
//  std::cout << "Orthtree extents are " << lower_bounds << std::endl << 
//      upper_bounds << std::endl;

    // Build tree: 
    std::size_t child_depth = 0; 

    this->root = new Node(center, box_length, child_depth, nullptr); 
    std::queue<Node*> node_queue({this->root}); 

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
        std::array<size_t, d> indices = Super::getLeafBoxIndices(points[k]);
        leaf_vectors[Super::getFlatIndex(indices)]->push_back(points[k]);
    }

    for(std::size_t k = 0; k < n_leaves; k++) {

        Leaf * leaf = static_cast<Leaf*>(node_queue.front());
        node_queue.pop(); 

        std::array<size_t, d> indices = Super::getLeafBoxIndices(leaf->center);
        leaf->points = leaf_vectors[Super::getFlatIndex(indices)];
    }
     
    delete[] leaf_vectors; 
}

template<typename Vector, std::size_t d>
void PointOrthtree<Vector, d>::toFile() {

    std::string geometry_filename = "geometry.dat";
    std::string data_filename =  "points.dat";

    ofstream geometry_file, data_file; 
    geometry_file.open(geometry_filename);
    data_file.open(data_filename);

    std::size_t n_node = 0;

    AOT::traverseBFSCore([&](Node * current) {

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
