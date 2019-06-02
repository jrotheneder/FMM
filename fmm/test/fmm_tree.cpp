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

#include <valgrind/callgrind.h> 

#include "../debugging.hpp" 
#include "../balanced_fmm_tree.hpp" 
#include "../adaptive_fmm_tree.hpp" 

using namespace fmm; 

int main(int argc, char *argv[]) {
    
    #define adaptive true
    const bool uniform = false; 
    const bool from_file = true; 

    size_t N = 10000;
    const size_t items_per_leaf = 50; 
    const size_t d = 2;
    const double eps = 1E-4; 
    const size_t order = ceil(log(1/eps) / log(2)); 
    const double extent = 1E-6;

    const size_t seed = 42; 
    srand(seed); 

    using Vec = Vector_<d>;
    using Src = PointSource_<d>;

    constexpr auto EPot = fields::electrostaticPotential_s<Vec, Src, d>;
    constexpr auto EFrc = fields::electrostaticForce_s<Vec, Src, d>;
    constexpr auto evalScalarInteraction = evaluateInteraction<Vec, Src, double>;
    constexpr auto evalVectorInteraction = evaluateInteraction<Vec, Src, Vec>;

    std::vector<Src> sources;

    if(from_file) {
        sources = getSourcesFromFile<d>("sources.dat");  
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
        AdaptiveFmmTree<d> fmm_tree(sources, items_per_leaf, eps);
    #else 
        BalancedFmmTree<d> fmm_tree(sources, items_per_leaf, eps);
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
    fmm_tree.toFile();

    Vec eval_point = sources[1].position;  
    std::cout << fmm_tree.evaluatePotential(eval_point) << std::endl; 
    std::cout << evalScalarInteraction(sources, eval_point, EPot) << std::endl;

    std::cout << fmm_tree.evaluateForcefield(eval_point) << std::endl;
    std::cout << evalVectorInteraction(sources, eval_point, EFrc) << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    CALLGRIND_TOGGLE_COLLECT;
    auto potentials = fmm_tree.evaluateParticlePotentialEnergies(); 
    CALLGRIND_TOGGLE_COLLECT;
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "fmm potential computation: " 
        << chrono_duration(t2-t1) << "s" << std::endl; 

    t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> ref_potentials(sources.size()); 

    #pragma omp parallel for
    for(std::size_t i = 0; i < sources.size(); ++i) {
        ref_potentials[i] =  sources[i].sourceStrength() 
            * evalScalarInteraction(sources, sources[i].position, EPot);
    }
    t2 = std::chrono::high_resolution_clock::now();

    std::vector<double> diffs(potentials.size()); 
    std::transform (
        potentials.begin(), potentials.end(), ref_potentials.begin(), 
        diffs.begin(), [](double a, double b) -> double { return std::abs((a-b)); } 
    );
    std::cout << "direct potential computation: " 
        << chrono_duration(t2-t1) << "s" << std::endl; 

    std::cout << "Mean relative potential deviation: " << 
        std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() << std::endl;

  
    t1 = std::chrono::high_resolution_clock::now();
    auto forces = fmm_tree.evaluateParticleForces(); 
    t2 = std::chrono::high_resolution_clock::now();
    std::vector<Vec> ref_forces(sources.size()); 
    std::cout << "fmm force computation: " 
        << chrono_duration(t2-t1) << "s" << std::endl; 


    t1 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(std::size_t i = 0; i < sources.size(); ++i) {
        ref_forces[i] =  sources[i].sourceStrength() * evalVectorInteraction(sources, 
            sources[i].position, EFrc);
    }

    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "direct force computation: " 
        << chrono_duration(t2-t1) << "s" << std::endl; 

    std::transform (
        forces.begin(), forces.end(), ref_forces.begin(), 
        diffs.begin(), [](Vec a, Vec b) -> double { return (a-b).norm()/b.norm(); } 
    );


    std::cout << "Mean relative force deviation: " << 
        std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() << std::endl;
  

    return 0;

}

