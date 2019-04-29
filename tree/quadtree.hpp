#include <vector> 
#include <complex> 
#include <cmath> 
#include <algorithm> 
#include <queue> 
#include <stack> 
#include <fstream> 
#include <cassert> 

#include "debugging.hpp" 

/*
template<typename Point>
struct NodeData {
    NodeData(); 
    virtual ~NodeData();
};

template<typename Point>
struct LeafData: NodeData<Point> {

    std::vector<Point> * points;

    LeafData(std::vector<Point> * points): points(points) {}
    virtual ~LeafData() { delete points; }

};
*/

template<typename Point, typename NodeData, int d>
class Quadtree {

    static const int n_children = pow(2,d); 

    struct Node {

        Point center;          // This nodes's center
        double box_length;     // This nodes's bounding box length
        std::size_t height;    // This nodes's height 

        Node * parent; 
        Node * children[n_children];    // Children 1 - 4 for NE, NW, SW, SE
        NodeData * data;

        Node(Point center, double box_length, std::size_t height, 
            Node * parent = nullptr): center(center), box_length(box_length), 
            height(height), parent(parent), data(nullptr) {}

        virtual ~Node() { 
            delete data; 
            for(int i = 0; i < n_children; i++) { delete children[i]; }
        }
    };

    struct Leaf: Node {

        std::vector<Point> * points;
        Leaf(Point center, double box_length, std::size_t height, 
            Node * parent = nullptr): Node(center, box_length, height, parent), 
            points(nullptr) {}
        virtual ~Leaf() { delete points; }
    };

public:

    Node * root; 

    const Point child_center_directions[4] = {
        {1,1}, {-1, 1}, {-1, -1}, {1, -1}}; //Directions NE, NW, SW, SE

    Quadtree(std::vector<Point> points, std::size_t s) {

        // Determine tree height, bounding box lenghts and center
        double height = ceil(log((double)points.size()/s)/log(n_children));
        //std::cout << "Quadtree height is " << height << endl;

        double x_min, x_max, y_min, y_max;
        std::tie(x_min, x_max, y_min, y_max) = Quadtree::getDataRange(points);

        double box_length = std::max(x_max-x_min, y_max-y_min);
        Point center = Point(0.5*(x_max + x_min), 0.5*(y_max + y_min));

        print_vec<double>({x_min, x_max, y_min, y_max}); 

        this->root = new Node(center, box_length, height, nullptr); 

        // Build tree: 
        std::size_t child_depth = 0; 
        std::queue<Node*> node_queue({this->root}); 

        while(child_depth < height) {

            std::size_t n_nodes = node_queue.size(); 
            assert(n_nodes == pow(n_children, child_depth));
            
            ++child_depth;

            for(std::size_t i = 0; i < n_nodes; i++) {

                Node * parent = node_queue.front();
                node_queue.pop(); 

                Point parent_center = parent->center;
                double child_box_length = parent->box_length/2;

//              std::cout << "parent_center = " << parent_center << std::endl; 
//              std::cout << "parent_box_length = " << parent->box_length << std::endl; 

                for(int j = 0; j < n_children; j++) {

                    Point child_center = parent_center + 
                        child_box_length/2 * child_center_directions[j];

//                  std::cout << ", child_center[" << j << "] = " <<  child_center
//                      << std::endl;
                    Node * child;

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

//      std::cout << "Tree stack size is " << n_leaves << std::endl;

        std::vector<Point> ** leaf_vectors = 
            new std::vector<Point>*[node_queue.size()];

        for(std::size_t i = 0; i < node_queue.size(); i++) {
            leaf_vectors[i] = new std::vector<Point>;
        }

        std::size_t n_boxes_per_dim =  pow(2, height); 
        for(std::size_t k = 0; k < points.size(); k++) {

            std::size_t i, j;
            std::tie(i, j) = this->getLeafBoxIndices(points[k]);

            assert(i < n_boxes_per_dim && j < n_boxes_per_dim);

            leaf_vectors[n_boxes_per_dim * i + j]->push_back(points[k]);
        }

        for(std::size_t k = 0; k < n_leaves; k++) {

            Leaf * leaf = static_cast<Leaf*>(node_queue.front());
            node_queue.pop(); 

            std::size_t i, j;
            std::tie(i, j) = getLeafBoxIndices(leaf->center);

//          std::cout << "Leaf " << k << " center = " << leaf->center 
//              << "indices (" << i << "," << j << ")" << std::endl; 

            leaf->points = leaf_vectors[n_boxes_per_dim * i + j];

            // TODO 
            leaf->children[0] = nullptr;
            leaf->children[1] = nullptr;
            leaf->children[2] = nullptr;
            leaf->children[3] = nullptr;
        }
         
        delete[] leaf_vectors; 
    }

    ~Quadtree() {

        // Alternative to the iterative implementation: 
        delete this->root;

        // I don't believe in recursion
//      std::stack<std::pair<Node*, bool>> node_stack({ {this->root, false} });

//      while(!node_stack.empty()) {

//          Node * current;
//          bool visited; 
//          std::tie(current, visited) = node_stack.top(); 

//          // Node has not been examined yet and has children => push children 
//          // to stack, flag as visited, process children
//          if(!visited && current->children[0] != nullptr) { 

//              node_stack.top().second = true; //flag node as visited

//              for(std::size_t i = 0; i < n_children; i++) {
//                  node_stack.push({current->children[i], false});  
//              }
//          }
//          // Node has been examined => all its children have been deleted
//          // => delete node, pop off stack. 
//          else { 
//              delete current; 
//              node_stack.pop();
//          }
//      }
    }

    // Return tuple of (x_min, x_max, y_min, y_max) of set of complex points {z = x+iy}
    static std::tuple<double,double,double,double> getDataRange(
            const std::vector<Point> & points) {
    
        double x_min = HUGE_VAL;
        double x_max = -HUGE_VAL;
        double y_min = HUGE_VAL;
        double y_max = -HUGE_VAL;

        std::for_each(points.begin(), points.end(), [&](Point p) { 

            double re = p.real(); 
            double im = p.imag();

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
    std::pair<std::size_t, std::size_t> getLeafBoxIndices(Point p) {

        Node * root = this->root;

        double box_x_min = root->center.real() - root->box_length/2;
        double box_y_min = root->center.imag() - root->box_length/2;

        assert(p.real() > box_x_min && p.real() < box_x_min + root->box_length); 
        assert(p.imag() > box_y_min && p.imag() < box_y_min + root->box_length); 

        double x_ratio = (p.real() - box_x_min) / root->box_length; 
        double y_ratio = (p.imag() - box_y_min) / root->box_length; 

        std::size_t n_boxes_per_dim =  pow(2, root->height); 

        std::size_t i = floor(x_ratio * n_boxes_per_dim);
        std::size_t j = floor(y_ratio * n_boxes_per_dim);

        return std::make_pair(i,j);
    }

    void toFile(std::string geometry_filename, std::string data_filename) {

        ofstream geometry_file, data_file; 

        geometry_file.open(geometry_filename);
        data_file.open(data_filename);

        std::queue<Node*> node_queue({this->root}); 
        std::size_t n_node = 0;

        while(!node_queue.empty()) {

            Node * current = node_queue.front();
            node_queue.pop(); 

            for(int i = 0; i < n_children; i++) {
                if(current->children[i] != nullptr)  {
                    node_queue.push(current->children[i]);
                }
            }

            double center_x = current->center.real();
            double center_y = current->center.imag();

            double box_x_min = center_x - 0.5 * current->box_length; 
            double box_y_min = center_y - 0.5 * current->box_length; 

            geometry_file << n_node << ", " 
                << this->root->height - current->height << ", "
                << center_x << ", "
                << center_y << ", "
                << box_x_min << ", "
                << box_y_min << ", "
                << std::endl;

            if(current->children[0] == nullptr) { // Leaf node => write data to file

                data_file << n_node; 
                Leaf * leaf = static_cast<Leaf*>(current);
                std::vector<Point> * points = leaf->points; 

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
