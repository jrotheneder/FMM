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

#include <omp.h> 
#include <valgrind/callgrind.h> 

#include "../debugging.hpp" 
#include "../balanced_fmm_tree.hpp" 
#include "../adaptive_fmm_tree.hpp" 

using namespace fmm; 
using namespace fmm::fields; 

int main(int argc, char *argv[]) {

    omp_set_num_threads(1); 

    size_t N = 7500;
    const size_t items_per_leaf = 400;
    const size_t d = 3;
    const double eps = 1E-4; 
    const double extent = 5; 
    const double force_smoothing_eps = 0;
    
    #define adaptive true
    const bool uniform = true; 
    const bool from_file = false; 
    const bool accuracy_check = true; 
    const bool tree_to_file = false; 
    const bool field_type = true;  // true for gravitational

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

    sourceLocationsToFile(sources, "logs/source_locs2.dat"); 
    
    auto t1 = std::chrono::high_resolution_clock::now();

    //CALLGRIND_TOGGLE_COLLECT;
    #if adaptive
        AdaptiveFmmTree<d, field_type> fmm_tree(sources, items_per_leaf, eps, 
            force_smoothing_eps);
    #else 
        std::cout << "this" << "\n";
        BalancedFmmTree<d, field_type> fmm_tree(sources, items_per_leaf, eps, 
            force_smoothing_eps);
    #endif
    //CALLGRIND_TOGGLE_COLLECT;

    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "N = " << N << "\n" 
        << "eps = " << eps << "\n"
        << "Order is " << fmm_tree.getOrder() << "\n"
        << (uniform ? "uniform" : "non-uniform") << " charge distribution. " 
        << "\nOrthtree is " << (adaptive ? "adaptive" : "non-adaptive") 
        << "\nheight is " << fmm_tree.getHeight() 
        << "\ncentered at " << fmm_tree.getCenter() 
        << "\nparent box size = " << fmm_tree.getBoxLength() << std::endl;

    if(tree_to_file) {
        fmm_tree.toFile();
    }

    Vec eval_point = sources[1].position;  
    std::cout << fmm_tree.evaluatePotential(eval_point) << std::endl; 
    std::cout << potential<d, field_type>(sources, eval_point) << std::endl;

    std::cout << fmm_tree.evaluateForcefield(eval_point) << std::endl;
    std::cout << forcefield<d, field_type>(sources, eval_point, 
            force_smoothing_eps) << std::endl;

    std::cout << "\n\nTree building took " << chrono_duration(t2-t1) << "s.\n";
    if(accuracy_check) {
        t1 = std::chrono::high_resolution_clock::now();
        auto potentials = fmm_tree.evaluateParticlePotentialEnergies(); 
        t2 = std::chrono::high_resolution_clock::now();
        std::cout << "fmm potential computation: " 
            << chrono_duration(t2-t1) << "s" << std::endl; 


        t1 = std::chrono::high_resolution_clock::now();
        CALLGRIND_TOGGLE_COLLECT;
        std::vector<Vector_<d>> forces = fmm_tree.evaluateParticleForces(); 
        CALLGRIND_TOGGLE_COLLECT;
        CALLGRIND_STOP_INSTRUMENTATION;    
        t2 = std::chrono::high_resolution_clock::now();

        std::cout << "fmm force computation: " 
            << chrono_duration(t2-t1) << "s" << std::endl; 

        std::vector<double> diffs(N); 
        std::vector<double> force_diffs(N);  

        t1 = std::chrono::high_resolution_clock::now();
        std::vector<double> ref_potentials
            = particlePotentialEnergies<d, field_type>(sources); 
        t2 = std::chrono::high_resolution_clock::now();
        double t_pot_dir = chrono_duration(t2-t1); 

        t1 = std::chrono::high_resolution_clock::now();
        std::vector<Vec> ref_forces 
            = particleForces<d, field_type>(sources, force_smoothing_eps); 
        t2 = std::chrono::high_resolution_clock::now();
        double t_frc_dir = chrono_duration(t2-t1); 

        std::transform (
                potentials.begin(), potentials.end(), ref_potentials.begin(), 
                diffs.begin(), [](double a, double b) -> double 
                { return std::abs(((a-b)/b)); } 
                );


        std::transform (
                forces.begin(), forces.end(), ref_forces.begin(), 
                force_diffs.begin(), [](Vec a, Vec b) -> double 
                { return (a-b).norm()/b.norm(); } 
                );

        std::cout << "direct potential computation: " 
            << t_pot_dir << "s" << std::endl; 
        std::cout << "direct force computation: " 
            << t_frc_dir << "s" << std::endl; 
        std::cout << "Mean relative potential deviation: " << 
            std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() 
            << std::endl;
        std::cout << "Max relative potential deviation: " << 
            *std::max_element(diffs.begin(), diffs.end()) << std::endl;
        std::cout << "Mean relative force deviation: " << 
            std::accumulate(force_diffs.begin(), force_diffs.end(), 0.0)
            / force_diffs.size() << std::endl;
        std::cout << "Max relative force deviation: " << 
            *std::max_element(force_diffs.begin(), force_diffs.end()) << std::endl;



//      t1 = std::chrono::high_resolution_clock::now();
//      std::vector<Vec> ref_forces2 
//          = particleForces<d, field_type>(sources, force_smoothing_eps); 
//      t2 = std::chrono::high_resolution_clock::now();
//      t_frc_dir = chrono_duration(t2-t1); 

//      t1 = std::chrono::high_resolution_clock::now();
//      std::vector<double> ref_potentials2 
//          = particlePotentialEnergies<d, field_type>(sources); 
//      t2 = std::chrono::high_resolution_clock::now();
//      t_pot_dir = chrono_duration(t2-t1); 

//      std::vector<double> diffs2(N); 
//      std::vector<double> force_diffs2(N);  

//      std::transform (
//              ref_potentials2.begin(), ref_potentials2.end(), ref_potentials.begin(), 
//              diffs2.begin(), [](double a, double b) -> double 
//              { return std::abs(((a-b)/b)); } 
//              );


//      std::transform (
//              ref_forces2.begin(), ref_forces2.end(), ref_forces.begin(), 
//              force_diffs2.begin(), [](Vec a, Vec b) -> double 
//              { return (a-b).norm()/b.norm(); } 
//              );

//      std::cout << "direct potential computation (par): " 
//          << t_pot_dir << "s" << std::endl; 
//      std::cout << "direct force computation (par): " 
//          << t_frc_dir << "s" << std::endl; 
//      std::cout << "Mean relative potential deviation: " << 
//          std::accumulate(diffs2.begin(), diffs2.end(), 0.0)/diffs2.size() 
//          << std::endl;
//      std::cout << "Mean relative force deviation: " << 
//          std::accumulate(force_diffs2.begin(), force_diffs2.end(), 0.0)
//          / force_diffs2.size() << std::endl;

    }
  

    return 0;

}

