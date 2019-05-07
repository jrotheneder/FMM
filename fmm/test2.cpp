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
    vector<PointCharge<d>> sources;

    for(size_t i = 0; i < N; i++) {
        Vector<d> v;      
        for(size_t j = 0; j < d; ++j) {
            v[j] =  64 * ((double) rand() / (RAND_MAX));
        }
        double q = (double) rand() / (RAND_MAX) * (i % 2 ? 1 : -1);
        PointCharge<d> src {v, q};

        sources.push_back(src); 
    }

//  PointOrthtree<Vector<d>, d>q(pts, 100);
    BalancedFmmTree<Vector<d>, PointCharge<d>, d>q(sources, 100, 1E-5);
    cout << "Orthtree height is " << q.getHeight() << ", centered at " <<
        q.getCenter() << endl;
    q.toFile();

    Vector<2> center = q.getCenter();
    fmm::MultipoleExpansion<Vector<2>, PointCharge<2>, 2> me(sources, center, 6); 

    Vector<2> eval_point{{96, 96}};
    cout << me.evaluatePotential(eval_point) << endl; 
    
    Vector<2> f = me.evaluateForcefield(eval_point);
    cout << f << endl;

    return 0;

}

