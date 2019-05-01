#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 

#include "debugging.hpp" 
#include "vector.hpp" 

#include "fmm_tree.hpp" 

//using namespace std;

int main(int argc, char *argv[]) {
     
    std::size_t N = 5000;
    const std::size_t d = 2;
    std::vector<Vector<d>> pts;

    for(std::size_t i = 0; i < N; i++) {
        Vector<d> v;      
        for(std::size_t j = 0; j < d; ++j) {
            v[j] =  64 * ((double) rand() / (RAND_MAX));
        }
        pts.push_back(Vector<d>(v)); 
    }

//  print_vec(pts); 
    
    FmmTree<Vector<d>, d>q(pts, 100);
    std::cout << "Orthtree height is " << q.getHeight() << ", centered at " <<
        q.getCenter() << std::endl;
    q.toFile();

    q.traverseBFSCore([&q](const FmmTree<Vector<d>, d>::Node * node) {
            std::cout << q.getHeight() - node->height << ", " 
            << node->center << std::endl; }); 


//  auto dirs =  Orthtree<Vector<3>, bool, 3>::getChildCenterDirections();
//  for(auto vec : dirs) {
//      std::cout << vec << "\n";
//  }
    
    
    

    return 0;

}

