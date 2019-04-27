#include <vector> 
#include <complex> 
#include <cmath> 
#include <algorithm> 
#include <queue> 
#include <stack> 
#include <fstream> 
#include <cassert> 

#include "debugging.hpp" 

struct QuadtreeData {
};

struct NodeData: QuadtreeData {
};

struct LeafData: QuadtreeData {

    std::vector<std::complex<double>> * points;
    LeafData(std::vector<std::complex<double>> * points): points(points) {}

    ~LeafData() {
        delete points;
    }

};

class Quadtree {
public:

    std::complex<double> center;  // This tree's center
    double box_length;            // This tree's bounding box length
    std::size_t height;           // This tree's height 
     
    Quadtree * parent; 
    Quadtree * children[4];       // Children 1 - 4 for NE, NW, SW, SE
    QuadtreeData * data;

    const std::complex<double> child_center_directions[4] = {
        {1,1}, {-1, 1}, {-1, -1}, {1, -1}}; //Directions NE, NW, SW, SE

    Quadtree(std::complex<double> center, double box_length, std::size_t height, 
        Quadtree * parent = nullptr): center(center), box_length(box_length), 
        height(height), parent(parent) {}
    
    Quadtree(std::vector<std::complex<double>> points, std::size_t s): 
        parent(nullptr) {

        // Determine tree height, bounding box lenghts and center
        this->height = ceil(log((double)points.size()/s)/log(4));
//      cout << "Quadtree height is " << this->height << endl;

        double x_min, x_max, y_min, y_max;
        std::tie(x_min, x_max, y_min, y_max) = Quadtree::getDataRange(points);

        this->box_length = std::max(x_max-x_min, y_max-y_min);
        this->center = std::complex<double>(0.5*(x_max + x_min), 0.5*(y_max + y_min));

        print_vec<double>({x_min, x_max, y_min, y_max}); 

        // Build tree: 
        std::size_t child_depth = 0; 
        std::queue<Quadtree*> tree_queue({this}); 

        while(child_depth < this->height) {

            std::size_t current_stack_size = tree_queue.size(); 
            assert(current_stack_size = pow(4, child_depth));
            
            ++child_depth;

            for(std::size_t i = 0; i < current_stack_size; i++) {

                Quadtree * parent = tree_queue.front();
                tree_queue.pop(); 

                std::complex<double> parent_center = parent->center;
                double child_box_length = parent->box_length/2;

//              std::cout << "parent_center = " << parent_center << std::endl; 
//              std::cout << "parent_box_length = " << parent->box_length << std::endl; 

                for(std::size_t j = 0; j < 4; j++) {

                    std::complex<double> child_center = parent_center + 
                        child_box_length/2 * child_center_directions[j];

//                  std::cout << ", child_center[" << j << "] = " <<  child_center
//                      << std::endl;
                    Quadtree * child = new Quadtree(child_center, child_box_length, 
                        this->height- child_depth, parent); 
                    tree_queue.push(child);

                    parent->children[j] = child; 
                }
            }
        }

        // Distribute points to leaves
        std::size_t current_stack_size = tree_queue.size(); 
        assert(current_stack_size = pow(4, child_depth));

//      std::cout << "Tree stack size is " << current_stack_size << std::endl;

        std::vector<std::complex<double>> ** leaf_vectors = 
            new std::vector<std::complex<double>>*[tree_queue.size()];

        for(std::size_t i = 0; i < tree_queue.size(); i++) {
            leaf_vectors[i] = new std::vector<std::complex<double>>;
        }

        std::size_t n_boxes_per_dim =  pow(2, this->height); 
        for(std::size_t k = 0; k < points.size(); k++) {

            std::size_t i, j;
            std::tie(i, j) = this->getLeafBoxIndices(points[k]);

            assert(i < n_boxes_per_dim && j < n_boxes_per_dim);

            leaf_vectors[n_boxes_per_dim * i + j]->push_back(points[k]);
        }

        for(std::size_t k = 0; k < current_stack_size; k++) {

            Quadtree * leaf = tree_queue.front();
            tree_queue.pop(); 


            std::size_t i, j;
            std::tie(i, j) = getLeafBoxIndices(leaf->center);

//          std::cout << "Leaf " << k << " center = " << leaf->center 
//              << "indices (" << i << "," << j << ")" << std::endl; 

            leaf->data = new LeafData(leaf_vectors[n_boxes_per_dim * i + j]);

            leaf->children[0] = nullptr;
            leaf->children[1] = nullptr;
            leaf->children[2] = nullptr;
            leaf->children[3] = nullptr;
        }
         
        delete[] leaf_vectors; 
    }

    ~Quadtree() {

        if(children[0] == nullptr) {
            LeafData *ld = static_cast<LeafData*>(this->data);
            delete ld;
        }
        else {
            for(std::size_t i = 0; i < 4; i++) {
                delete children[i];  
            }
        }
    }

    // Return tuple of (x_min, x_max, y_min, y_max) of set of complex points {z = x+iy}
    static std::tuple<double,double,double,double> getDataRange(
            const std::vector<std::complex<double>> & points) {
    
        double x_min = HUGE_VAL;
        double x_max = -HUGE_VAL;
        double y_min = HUGE_VAL;
        double y_max = -HUGE_VAL;

        std::for_each(points.begin(), points.end(), [&](std::complex<double> c) { 

            double re = c.real(); 
            double im = c.imag();

            x_min = re < x_min ? re : x_min;
            x_max = re > x_max ? re : x_max;
            y_min = im < y_min ? im : y_min; 
            y_max = im > y_max ? im : y_max; 

        }); 

        // Expand bounding box s.t. all points lie completely within it (and not
        // on boundaries. This simplifies the index calculations in getLeafBoxIndices().
        //
        double growthFactor = 1E-10;
        x_min -= x_min * growthFactor;
        x_max += x_max * growthFactor;
        y_min -= y_min * growthFactor;
        y_max += y_max * growthFactor;

        return make_tuple(x_min, x_max, y_min, y_max);
    }

    // Given a point, determines indices (i,j) of leaf that this point belongs
    // to, where leaves are indexed according to the part of they domain they
    // own (smallest x -> i = 0, smallest y -> j = 0 etc...) 
    // This method should be called on the root object!
    std::pair<std::size_t, std::size_t> getLeafBoxIndices(std::complex<double> p) {

        double box_x_min = this->center.real() - this->box_length/2;
        double box_y_min = this->center.imag() - this->box_length/2;

        assert(p.real() > box_x_min && p.real() < box_x_min + this->box_length); 
        assert(p.imag() > box_y_min && p.imag() < box_y_min + this->box_length); 

        double x_ratio = (p.real() - box_x_min) / box_length; 
        double y_ratio = (p.imag() - box_y_min) / box_length; 

        std::size_t n_boxes_per_dim =  pow(2, this->height); 

        std::size_t i = floor(x_ratio * n_boxes_per_dim);
        std::size_t j = floor(y_ratio * n_boxes_per_dim);

        return std::make_pair(i,j);
    }

    void toFile(std::string geometry_filename, std::string data_filename) {

        ofstream geometry_file, data_file; 
        geometry_file.open(geometry_filename);
        data_file.open(data_filename);

        std::queue<Quadtree*> tree_queue({this}); 
        std::size_t n_node = 0;

        while(!tree_queue.empty()) {

            Quadtree * tree = tree_queue.front();
            tree_queue.pop(); 

            for(std::size_t i = 0; i < 4; i++) {
                if(tree->children[i] != nullptr)  {
                    tree_queue.push(tree->children[i]);
                }
            }
            
            double center_x = tree->center.real();
            double center_y = tree->center.imag();

            double box_x_min = center_x - 0.5 * tree->box_length; 
            double box_y_min = center_y - 0.5 * tree->box_length; 

            geometry_file << n_node << ", " 
                << this->height - tree->height << ", "
                << center_x << ", "
                << center_y << ", "
                << box_x_min << ", "
                << box_y_min << ", "
                << std::endl;

            if(tree->children[0] == nullptr) { // Leaf node => write data to file
//              std::cout << "Node " << n_node << " (leaf) writing data" << std::endl;

                data_file << n_node; 
                LeafData *ld = static_cast<LeafData*>(tree->data);
                std::vector<std::complex<double>> * points = ld->points; 

                for(std::size_t i = 0; i < points->size(); i++) {
                    data_file << ", " << (*points)[i].real() << ", " 
                        << (*points)[i].imag(); 
                }
                data_file << std::endl;
            }

            ++n_node;
        }
        
        geometry_file.close();
        data_file.close();
    }

};
