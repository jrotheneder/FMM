#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <fstream> 
#include <string> 
#include <sstream> 
#include <cstdlib> 
#include <cmath> 
#include <chrono> 
#include <filesystem> 

#include <valgrind/callgrind.h> 

#include "../debugging.hpp" 
#include "../balanced_fmm_tree.hpp" 
#include "../adaptive_fmm_tree.hpp" 

using namespace fmm; 
using namespace fmm::fields; 

int main(int argc, char *argv[]) {

    size_t N = 1000;
    const size_t items_per_leaf = 90;
    const size_t d = 2;
    const double eps = 1E-4; 
    const size_t order = ceil(log(1/eps) / log(2));
    const double extent = 1; 
    const double force_smoothing_eps = 0;
    
    #define adaptive true
    const bool uniform = false; 
    const bool from_file = (d == 2) && true; 
    const bool accuracy_check = true; 
    const bool tree_to_file = true; 
    const bool field_type = false;  // true for gravitational

    const size_t seed = 42; 
    srand(seed); 

    using Vec = Vector_<d>;
    using Src = PointSource_<d>;

    std::vector<Src> sources;

    if(from_file) {
        if(!std::filesystem::exists("input/sources.dat")) {
            throw std::runtime_error("File not found.");  
        }
        sources = readSourcesFromFile<d>("input/sources.dat");  
        N = sources.size(); 
    }
    else {
        Vec center{}; // Origin 
        Vec shift(8);
        auto centers = AbstractOrthtree<Vec, d>::getChildCenterDirections(); 
        for(std::size_t i = 0; i < (1 << d) ; ++i) {
            centers[i] = centers[i] * (((rand() % 2) ? 1 : 1.5) * extent); 
        }

        for(size_t i = 0; i < N; i++) {
            Vec v;      
            for(size_t j = 0; j < d; ++j) {
                v[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
            }

            if(!uniform) {
                v += centers[rand() % AbstractOrthtree<Vec, d>::n_children];
            }

            double q = (double) rand() / (RAND_MAX) * (i % 2 ? 1 : -1);
            Src src {v, q};

            sources.push_back(src); 
        }
    }

    
    auto t1 = std::chrono::high_resolution_clock::now();

    CALLGRIND_TOGGLE_COLLECT;
    #if adaptive
        AdaptiveFmmTree<d, field_type> fmm_tree(sources, items_per_leaf, eps, 
            force_smoothing_eps);
    #else 
        BalancedFmmTree<d, field_type> fmm_tree(sources, items_per_leaf, eps, 
            force_smoothing_eps);
    #endif
    CALLGRIND_TOGGLE_COLLECT;

    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "N = " << N << ", " << (uniform ? "uniform" : "non-uniform") 
        << " charge distribution. Orthtree is " 
        << (adaptive ? "adaptive" : "non-adaptive") << ", height is " 
        << fmm_tree.getHeight() << ", centered at " << fmm_tree.getCenter() 
        << " with parent box size = " << fmm_tree.getBoxLength() <<  std::endl;

    std::cout << "Tree building took " << chrono_duration(t2-t1) << "s, Order is " 
        << order << std::endl;
    if(tree_to_file) {
        fmm_tree.toFile();
    }

    Vec eval_point = sources[1].position;  
    std::cout << fmm_tree.evaluatePotential(eval_point) << std::endl; 
    std::cout << potential<d, field_type>(sources, eval_point, 
            force_smoothing_eps) << std::endl;

    std::cout << fmm_tree.evaluateForcefield(eval_point) << std::endl;
    std::cout << forcefield<d, field_type>(sources, eval_point, 
            force_smoothing_eps) << std::endl;

    if(accuracy_check) {
//      t1 = std::chrono::high_resolution_clock::now();
//      auto potentials = fmm_tree.evaluateParticlePotentialEnergies(); 
//      t2 = std::chrono::high_resolution_clock::now();
//      std::cout << "fmm potential computation: " 
//          << chrono_duration(t2-t1) << "s" << std::endl; 


        t1 = std::chrono::high_resolution_clock::now();
        CALLGRIND_TOGGLE_COLLECT;
        auto forces = fmm_tree.evaluateParticleForces(); 
        CALLGRIND_TOGGLE_COLLECT;
        CALLGRIND_STOP_INSTRUMENTATION;    
        t2 = std::chrono::high_resolution_clock::now();

        std::cout << "fmm force computation: " 
            << chrono_duration(t2-t1) << "s" << std::endl; 

//      std::vector<Vec> ref_forces(sources.size()); 
//      std::vector<double> ref_potentials(sources.size()); 
//      t1 = std::chrono::high_resolution_clock::now();
//      #pragma omp parallel for
//      for(std::size_t i = 0; i < sources.size(); ++i) {
//          ref_potentials[i] =  sources[i].sourceStrength() 
//              * potential<d, field_type>(sources, sources[i].position,
//                      force_smoothing_eps);
//      }
//      t2 = std::chrono::high_resolution_clock::now();
//      double t_pot_dir = chrono_duration(t2-t1); 

//      std::vector<double> diffs(N); 
//      std::transform (
//              potentials.begin(), potentials.end(), ref_potentials.begin(), 
//              diffs.begin(), [](double a, double b) -> double 
//              { return std::abs(((a-b)/b)); } 
//              );

//      std::vector<double> force_diffs(N);  
//      t1 = std::chrono::high_resolution_clock::now();
//      #pragma omp parallel for
//      for(std::size_t i = 0; i < sources.size(); ++i) {
//          ref_forces[i] =  sources[i].sourceStrength() 
//              * forcefield<d, field_type>(sources, 
//              sources[i].position, force_smoothing_eps);
//      }
//      t2 = std::chrono::high_resolution_clock::now();
//      double t_frc_dir = chrono_duration(t2-t1); 

//      std::transform (
//              forces.begin(), forces.end(), ref_forces.begin(), 
//              force_diffs.begin(), [](Vec a, Vec b) -> double 
//              { return (a-b).norm()/b.norm(); } 
//              );

//      std::cout << "direct potential computation: " 
//          << t_pot_dir << "s" << std::endl; 
//      std::cout << "direct force computation: " 
//          << t_frc_dir << "s" << std::endl; 
//      std::cout << "Mean relative potential deviation: " << 
//          std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() 
//          << std::endl;
//      std::cout << "Mean relative force deviation: " << 
//          std::accumulate(force_diffs.begin(), force_diffs.end(), 0.0)
//          / force_diffs.size() << std::endl;

    }
  

    return 0;

}

