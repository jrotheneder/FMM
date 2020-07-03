#include "../adaptive_fmm_tree.hpp"

using namespace fmm;

int main(int argc, char *argv[]) {

    const size_t N = 10000;
    const size_t items_per_leaf = 50; 
    const size_t d = 2;
    const double eps = 1E1; 
    const double extent = 32;

    using Vec = Vector_<d>;
    using Src = PointSource_<d>;

    const bool field_type = 1; // gravitational

    std::vector<Src> sources;
    for(size_t i = 0; i < N; i++) {
        Vec v;      
        for(size_t j = 0; j < d; ++j) {
            v[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }
        double q = (double) rand() / (RAND_MAX) * (i % 2 ? 1 : -1);
        sources.push_back({v,q}); 
    }

    AdaptiveFmmTree<d, field_type> t(sources, items_per_leaf, eps);  
    Vec eval_point = sources[0].position; 

    std::cout << "Potential at " << eval_point << " evaluated via the FMM is " << 
       t.evaluatePotential(eval_point) << " vs " 
       << fields::potential<d, field_type>(sources, eval_point) 
       << " when evaluated directly. \n"; 
      



}

