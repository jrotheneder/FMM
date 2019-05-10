#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 

#include "debugging.hpp" 
#include "vector.hpp" 

//#include "point_orthtree.hpp" 
#include "fmm_tree.hpp" 

//using namespace std;
using namespace fmm; 

int main(int argc, char *argv[]) {
     
    size_t N = 5000;
    const size_t d = 2;
    const size_t order = 10; 
    const size_t seed = 42; // Reference data is hardcoded and requires this seed
    srand(seed); 

    vector<PointCharge<d>> sources;
    Vector<d> center{}; // Origin 
    Vector<d> ones; ones.fill(1);
    double extent = 32;

    for(size_t i = 0; i < N; i++) {
        Vector<d> v;      
        for(size_t j = 0; j < d; ++j) {
            v[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }
        double q = (double) rand() / (RAND_MAX) * (i % 2 ? 1 : -1);
        PointCharge<d> src {v, q};

        sources.push_back(src); 
    }

//  PointOrthtree<Vector<d>, d>q(pts, 100);
    BalancedFmmTree<Vector<d>, PointCharge<d>, d>q(sources, 100, 1E-5);
    std::cout << "Orthtree height is " << q.getHeight() << ", centered at " <<
        q.getCenter() << std::endl;
    q.toFile();

//  q.traverseBFSCore([&q](const AbstractOrthtree<Vector<d>, d>::Node * node) {
//          std::cout << q.getHeight() - node->height << ", " 
//          << node->center << std::endl; }); 


//  auto dirs =  Orthtree<Vector<3>, bool, 3>::getChildCenterDirections();
//  for(auto vec : dirs) {
//      std::cout << vec << "\n";
//  }
    
    return 0;

}
