#include <vector> 
#include <cmath> 
#include <algorithm> 
#include <array> 
#include <queue> 
#include <stack> 
#include <fstream> 
#include <bitset> 
#include <cassert> 

#include "debugging.hpp" 

template<typename Vector, typename NodeData, int d>
class Quadtree {

    static const int n_children = pow(2,d); 

    struct Node {

        Vector center;         // This nodes's center
        double box_length;     // This nodes's bounding box length
        std::size_t height;    // This nodes's height 

        Node * parent; 
        std::array<Node*, n_children> children; 
        NodeData * data;

        Node(Vector center, double box_length, std::size_t height, 
            Node * parent = nullptr): center(center), box_length(box_length), 
            height(height), parent(parent), children(), data() {}

        virtual ~Node() { 
            for(int i = 0; i < n_children; i++) { delete children[i]; }
            delete data; 
        }
    };

    struct Leaf: Node {

        std::vector<Vector> * points;
        Leaf(Vector center, double box_length, std::size_t height, 
            Node * parent): Node(center, box_length, height, parent), points() {}
        virtual ~Leaf() { delete points; }
    };

public:

    Node * root; 
    std::array<Vector, n_children> child_center_directions;

    Quadtree(std::vector<Vector> points, std::size_t s) {

        // Determine tree height, bounding box lenghts and center as well as
        // the directions in which the child centers lie relative to the center
        // of a box
        double height = ceil(log((double)points.size()/s)/log(n_children));

        Vector lower_bounds, upper_bounds;  
        std::tie(lower_bounds, upper_bounds) = getDataRange(points);
        auto extents = Vector{upper_bounds - lower_bounds}.data(); 
        double box_length = *std::max_element(extents.begin(), extents.end());

        Vector center = 0.5 * (lower_bounds + upper_bounds); 

        this->child_center_directions = getChildCenterDirections();

        std::cout << "Quadtree height is " << height << endl;
        std::cout << "Quadtree extents are " << lower_bounds << std::endl << 
            upper_bounds << std::endl;

        // Build tree: 
        this->root = new Node(center, box_length, height, nullptr); 

        std::size_t child_depth = 0; 
        std::queue<Node*> node_queue({this->root}); 

        while(child_depth < height) {

            std::size_t n_nodes = node_queue.size(); 
            assert(n_nodes == pow(n_children, child_depth));
            
            ++child_depth;

            for(std::size_t i = 0; i < n_nodes; i++) {

                Node * parent = node_queue.front();
                node_queue.pop(); 

                Vector parent_center = parent->center;
                double child_box_length = parent->box_length/2;

//              std::cout << "parent_center = " << parent_center << std::endl; 
//              std::cout << "parent_box_length = " << parent->box_length << std::endl; 

                for(int j = 0; j < n_children; j++) {

                    Node * child;
                    Vector child_center = parent_center + 
                        child_box_length/2 * child_center_directions[j];

//                  std::cout << ", child_center[" << j << "] = " <<  child_center
//                      << std::endl;

                    if(child_depth < height) {
                        child = new Node(child_center, child_box_length, 
                            height-child_depth, parent); 
                    }
                    else { // Treat leaves seperately
                        child = new Leaf(child_center, child_box_length, 
                            height-child_depth, parent); 
                    }

                    node_queue.push(child);
                    parent->children[j] = child; 
                }
            }
        }

        // Distribute points to leaves
        std::size_t n_leaves = node_queue.size(); 
        assert(n_leaves == pow(n_children, child_depth));

        std::vector<Vector> ** leaf_vectors = new std::vector<Vector>*[n_leaves];
        for(std::size_t i = 0; i < n_leaves; i++) {
            leaf_vectors[i] = new std::vector<Vector>;
        }

        for(std::size_t k = 0; k < points.size(); k++) {
            std::array<size_t, d> indices = getLeafBoxIndices(points[k]);
            leaf_vectors[getFlatIndex(indices)]->push_back(points[k]);
        }

        for(std::size_t k = 0; k < n_leaves; k++) {

            Leaf * leaf = static_cast<Leaf*>(node_queue.front());
            node_queue.pop(); 

            std::array<size_t, d> indices = getLeafBoxIndices(leaf->center);
            leaf->points = leaf_vectors[getFlatIndex(indices)];
        }
         
        delete[] leaf_vectors; 
    }

    ~Quadtree() { delete this->root; }

    static std::array<Vector, n_children> getChildCenterDirections() {

        std::array<Vector, n_children> child_center_directions;

        //Use a length d bit vector to represent any of the n_children = 2^d 
        //possible directions in which the centers of a given node's children 
        //may lie. A 0 in the d-i-1-st position of the bit vector corresponds to 
        //a direction vector whose i-th component is -1, a 1 in the same
        //position to a direction vector whose i-th component is +1.

        for(std::size_t i = 0; i < n_children; ++i) {

            std::bitset<d> direction_bits = bitset<d>(i); 

            Vector v;  
            for(std::size_t j = 0; j < d; ++j) {
                v[j] = direction_bits[d-j-1] ? 1 : -1;   
            }

            child_center_directions[i] = v;  
        }

        // Outputs: 
        // 2D: {-1, -1}, {-1, 1}, {1, -1}, {1, 1}
        // 3D: {-1, -1, -1}, {-1, -1, 1 }, {-1, 1, -1}, {-1, 1, 1}, {1, -1, -1}, 
        //     {1, -1, 1}, {1, 1, -1}, {1, 1, 1}

        return child_center_directions; 
    }

    // Return tuple of (x_min, x_max, y_min, y_max, ...) of set of input vectors
    // (x, y, ...)
    std::tuple<Vector, Vector> getDataRange(const std::vector<Vector> & points) 
        const {
    
        Vector lower_bounds, upper_bounds; 
        lower_bounds.data().fill(HUGE_VAL);
        upper_bounds.data().fill(-HUGE_VAL);

        std::for_each(points.begin(), points.end(), [&](Vector p) { 
            for(int i = 0; i < d; i++) {
                lower_bounds[i] = p[i] < lower_bounds[i] ? p[i] : lower_bounds[i];  
                upper_bounds[i] = p[i] > upper_bounds[i] ? p[i] : upper_bounds[i];  
            }
        }); 

        // Expand bounding box s.t. all points lie completely within it (and not
        // on boundaries. This simplifies the index calculations in getLeafBoxIndices().

        double paddingFactor = 1E-5;
        lower_bounds = lower_bounds - paddingFactor * lower_bounds;
        upper_bounds = upper_bounds + paddingFactor * upper_bounds;

        return make_tuple(lower_bounds, upper_bounds);
    }

    // Given a point, determines indices (i_1,i_2,...) of leaf that this point belongs
    // to, where leaves are indexed according to the part of they domain they
    // own (smallest x coordinate -> i_1 = 0, smallest y coordinate -> i_2 = 0 etc...) 
    std::array<std::size_t, d> getLeafBoxIndices(Vector p) const {

        Node * root = this->root;
        double box_length = root->box_length; 

        Vector diagonal;
        diagonal.data().fill(1);

        Vector lower_left_box_corner = root->center - box_length/2 * diagonal; 
        Vector upper_right_box_corner = root->center + box_length/2 * diagonal; 

        for(int i = 0; i < d; i++) {
            assert(p[i] > lower_left_box_corner[i] 
                && p[i] < upper_right_box_corner[i]); 
        }
        
        Vector ratios = (p - lower_left_box_corner) * (1 / box_length);

        const std::size_t n_boxes_per_dim = pow(2, root->height); 

        std::array<std::size_t, d> indices; 
        for(std::size_t j = 0; j < d; j++) {
            indices[j] = floor(ratios[j] * n_boxes_per_dim);
        }
        
        return indices;
    }

    void toFile(std::string geometry_filename, std::string data_filename) const {

        ofstream geometry_file, data_file; 
        geometry_file.open(geometry_filename);
        data_file.open(data_filename);

        std::queue<Node*> node_queue({this->root}); 
        std::size_t n_node = 0;
        double tree_height = this->root->height;

        Vector diagonal;
        diagonal.data().fill(1);

        while(!node_queue.empty()) {

            Node * current = node_queue.front();
            node_queue.pop(); 

            for(int i = 0; i < n_children; i++) {
                if(current->children[i] != nullptr)  {
                    node_queue.push(current->children[i]);
                }
            }

            Vector center = current->center;
            Vector lower_left_box_corner = center - current->box_length/2 * diagonal; 

            geometry_file << n_node << ", " << tree_height - current->height;

            for(auto coord : center.data()) { geometry_file << ", " << coord; }
            for(auto coord : lower_left_box_corner.data()) { 
                geometry_file << ", " << coord; 
            }
            geometry_file << std::endl;

            if(current->children[0] == nullptr) { // Leaf node => write leaf data to file

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
        }
        
        geometry_file.close();
        data_file.close();
    }

    std::size_t getFlatIndex(const std::array<std::size_t, d>& indices) const {
        
        const std::size_t n_boxes_per_dim = pow(2, this->root->height);   

        std::size_t flat_index = 0; 
        std::size_t offset_factor = 1; 

        for(std::size_t i = 0; i < d; ++i) {
            flat_index += offset_factor * indices[d-i-1];  
            offset_factor *= n_boxes_per_dim; 
        }

        return flat_index;
    }

};
