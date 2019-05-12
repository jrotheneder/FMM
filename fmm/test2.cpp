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

using namespace std;
using namespace fmm; 

int main(int argc, char *argv[]) {
     
    size_t N = 5000;
    const size_t d = 2;
    const double eps = 1E-5; 
    const size_t order = ceil(log(1/eps) / log(2)); 
    const size_t seed = 42; 
    srand(seed); 

    vector<PointCharge<d>> sources;
    Vector<d> center{}; // Origin 
    double extent = 32;

    Vector<d> eval_point; eval_point.fill(96); 
    Vector<d> shift; shift.fill(8);
    const Vector<d> shifted_center = center + shift;

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
    cout << "Orthtree height is " << q.getHeight() << ", centered at " <<
        q.getCenter() << " with parent box size = " << q.getBoxLength() <<  endl;
    q.toFile();

    fmm::MultipoleExpansion<Vector<d>, PointCharge<d>, d> me(center, order, sources); 
    std::vector<fmm::MultipoleExpansion<Vector<d>, PointCharge<d>, d>*> vme{&me};
    fmm::MultipoleExpansion<Vector<d>, PointCharge<d>, d> se(shifted_center, vme); 

    /*
    cout << "Order is " << order << endl;
    cout << me.evaluatePotential(eval_point) << endl; 
    cout <<  me.evaluateForcefield(eval_point) << endl;

    cout << se.evaluatePotential(eval_point) << endl; 
    cout << se.evaluateForcefield(eval_point) << endl;
    */

    /*
    for(unsigned i = 0; i < 8; ++i) {
        for(unsigned j = 0; j < 8; ++j) {
            std::cout << i << "\t" << j << "\t" << q.getMortonIndex({i,j}) << endl; 
        }
    }
    */
    

    return 0;

}

