#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 

#include "debugging.hpp" 
#include "vector.hpp" 

#include "balanced_point_orthtree.hpp" 
#include "adaptive_point_orthtree.hpp" 

//using namespace std;

int main(int argc, char *argv[]) {
     
    size_t N = 2000;
    const int s = 20;
    const size_t d = 2;
    const size_t seed = 42; // Reference data is hardcoded and requires this seed
    srand(seed); 

    using Vec = Vector<d>;

    double extent = 32;
    auto centers = AbstractOrthtree<Vec, d>::getChildCenterDirections(); 
    for(std::size_t i = 0; i < (1 << d) ; ++i) {
        centers[i] = centers[i] * (((rand() % 2) ? 1 : 1.5) * extent); 
    }
     
    Vec ones; ones.fill(1);

    std::vector<Vec> sources;
    for(size_t i = 0; i < N; i++) {
        Vec v;      
        for(size_t j = 0; j < d; ++j) {
            v[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }
        v += centers[rand() % AbstractOrthtree<Vec, d>::n_children];

        sources.push_back(v); 
    }

    //BalancedPointOrthtree<Vec, d>q(sources, s);
    AdaptivePointOrthtree<Vec, d>aq(sources, s);

    std::cout << "Orthtree height is " << aq.getHeight() << ", centered at " <<
        aq.getCenter() << std::endl;
    aq.toFile();

//  q.traverseBFSCore([&q](const AbstractOrthtree<Vec, d>::Node * node) {
//          std::cout << q.getHeight() - node->height << ", " 
//          << node->center << std::endl; }); 


//  auto dirs =  Orthtree<Vector<3>, bool, 3>::getChildCenterDirections();
//  for(auto vec : dirs) {
//      std::cout << vec << "\n";
//  }
    
    return 0;

}
