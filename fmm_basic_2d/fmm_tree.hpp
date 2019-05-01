#ifndef FMM_TREE_H
#define FMM_TREE_H

#include "orthtree.hpp"

template<typename Vector, std::size_t d>
class FmmTree: public Orthtree<Vector, d> {

public:
    using Super = Orthtree<Vector, d>;
    using Node = typename Super::Node; 

    struct FmmData: Super::NodeData {

//      struct MultipoleExpansion;
//      struct LocalExpansion;

//      MultipoleExpansion multipole_expansion; 
//      LocalExpansion local_expansion;

        std::vector<Node*> interaction_list;
        std::vector<Node*> near_neighbours;
        virtual ~FmmData() {};
    };

public:

    FmmTree(std::vector<Vector> points, std::size_t s): Super(points, s) {

        if(d > 3lu) { 
            throw std::runtime_error("Determining near neighbours from "
                "euclidean center seperation is only valid for d <= 3 "
                "dimensions, switch implementation."); 
        }

        // Maximal distance to near neighbours, in units of box_length and with 
        // padding to avoid numerical issues
        const double max_neighbour_distance = 1.1 * sqrt(d); 
    
        Super::traverseBFSCore([max_neighbour_distance](Node * node) {

            FmmData * fmm_data = new FmmData();
            node->data = fmm_data; 

            auto& interaction_list = fmm_data->interaction_list;
            auto& near_neighbours = fmm_data->near_neighbours;

            Node* parent = node->parent;
            Vector& center = node->center;

            double box_max_neighbour_distance = max_neighbour_distance 
                * node->box_length;

            std::cout << box_max_neighbour_distance << std::endl;

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

//  void toFile(std::string geometry_filename, std::string data_filename, 
//      std::string neighbours_filename, std::string interaction_filename) {
    virtual void toFile() override {
        std::string neighbours_filename = "neighbours.dat";
        std::string interaction_filename = "interactions.dat";

        ofstream neighbours_file, interaction_file; 
        neighbours_file.open(neighbours_filename);
        interaction_file.open(interaction_filename);

        std::size_t n_node = 0;

        Super::toFile(); 
        Super::traverseBFSCore([&](Node* node) {

            auto node_data = static_cast<FmmData*>(node->data); 
            auto& interaction_list = node_data->interaction_list;
            auto& near_neighbours = node_data->near_neighbours;

            Vector& center = node->center;
            std::size_t depth = this->getHeight() - node->height;
            double box_length = node->box_length;

            std::stringstream head; 

            head << n_node++ << ", " << depth << ", " << box_length;
            for(auto coord : center.data()) { head << ", " << coord; }

            neighbours_file << head.str(); 
            interaction_file << head.str(); 

            for(auto node : near_neighbours) {
                for(auto coord : node->center.data()) {
                    neighbours_file << ", " << coord; 
                }
            }
            neighbours_file << std::endl;
            
            for(auto node : interaction_list) {
                for(auto coord : node->center.data()) {
                    interaction_file << ", " << coord; 
                }
            }
            interaction_file << std::endl;

        }); 
        
        neighbours_file.close();
        interaction_file.close(); 
    }
};

#endif
