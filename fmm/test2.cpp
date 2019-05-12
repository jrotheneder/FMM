#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <chrono> 

#include "debugging.hpp" 
#include "vector.hpp" 

//#include "point_orthtree.hpp" 
#include "fmm_tree.hpp" 

using namespace fmm; 

int main(int argc, char *argv[]) {
    

    const size_t N = 1000;
    const size_t items_per_leaf = 50; 
    const size_t d = 2;
    const double eps = 1E-5; 
    const size_t order = ceil(log(1/eps) / log(2)); 
    const size_t seed = 42; 
    srand(seed); 

    using Vec = Vector<d>;
    using Src = PointCharge<d>;

    constexpr auto EPot = fields::electrostaticPotential_s<Vec, Src, d>;
    constexpr auto EFrc = fields::electrostaticForce_s<Vec, Src, d>;
    constexpr auto evalScalarInteraction = evaluateInteraction<Vec, Src, double>;
    constexpr auto evalVectorInteraction = evaluateInteraction<Vec, Src, Vec>;

    std::vector<Src> sources;
    Vec center{}; // Origin 
    double extent = 32;

    Vec shift(8);

    for(size_t i = 0; i < N; i++) {
        Vec v;      
        for(size_t j = 0; j < d; ++j) {
            v[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }
        double q = (double) rand() / (RAND_MAX) * (i % 2 ? 1 : -1);
        Src src {v, q};

        sources.push_back(src); 
    }

//  PointOrthtree<Vec, d>q(pts, 100);
    BalancedFmmTree<Vec, Src, d>q(sources, items_per_leaf, eps);
    std::cout << "Orthtree height is " << q.getHeight() << ", centered at " <<
        q.getCenter() << " with parent box size = " << q.getBoxLength() <<  std::endl;
    std::cout << "Order is " << order << std::endl;
    q.toFile();

    Vec eval_point = sources[0].position;  
    std::cout << q.evaluatePotential(eval_point) << std::endl; 
    std::cout << evalScalarInteraction(sources, eval_point, EPot) << std::endl;

    std::cout << q.evaluateForcefield(eval_point) << std::endl;
    std::cout << evalVectorInteraction(sources, eval_point, EFrc) << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto potentials = q.evaluateParticlePotentials(); 
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "took " << chrono_duration(t2-t1) << std::endl; 

    std::vector<double> ref_potentials(sources.size()); 
    for(std::size_t i = 0; i < sources.size(); ++i) {
        ref_potentials[i] =  evalScalarInteraction(sources, 
            sources[i].position, EPot);
    }

    std::vector<double> diffs(potentials.size()); 
    std::transform (
        potentials.begin(), potentials.end(), ref_potentials.begin(), 
        diffs.begin(), [](double a, double b) -> double { return std::abs(a-b); } 
    );

    std::cout << "Mean abs potential deviation: " << 
        std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() << std::endl;

    auto forces = q.evaluateParticleForcefields(); 
    std::vector<Vec> ref_forces(sources.size()); 

    for(std::size_t i = 0; i < sources.size(); ++i) {
        ref_forces[i] =  evalVectorInteraction(sources, 
            sources[i].position, EFrc);
    }

    std::transform (
        forces.begin(), forces.end(), ref_forces.begin(), 
        diffs.begin(), [](Vec a, Vec b) -> double { return (a-b).norm(); } 
    );

    std::cout << "Mean abs force deviation: " << 
        std::accumulate(diffs.begin(), diffs.end(), 0.0)/diffs.size() << std::endl;

    return 0;

}

