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

using namespace std; 
using namespace fmm; 
using namespace fmm::fields; 

const size_t d = 3;
const bool field_type = true;  // true for gravitational

template<typename T>
double vec_avg(vector<T> v) {
    return accumulate(v.begin(), v.end(), 0.0) / v.size(); 
}

template<typename T>
double vec_max(vector<T> v) { return *max_element(v.begin(), v.end()); }

template<unsigned d>
vector<PointSource_<d>> generateUniformDistribution(Vector_<d> center, double extent, 
        int N) {

    vector<PointSource_<d>> sources; 

    for(int i = 0; i < N; i++) {

        Vector_<d> v = center;      
        for(unsigned j = 0; j < d; ++j) {
            v[j] +=  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }

        int sign = i % 2 ? 1 : (!field_type ? -1 : 1); 
        double q = (double) rand() / (RAND_MAX) * sign;

        sources.push_back({v, q}); 
    }

    return sources; 
}

template<unsigned d>
vector<PointSource_<d>> generateNonUniformDistribution(double extent, int N) {

    double shift_dist = 8 * extent; 
    
    vector<Vector_<d>> centers; 
    if constexpr(d == 2) { // Equilateral triangle
        centers = {
            {{ 0,1 }},
            {{ sqrt(3)/2.,	-0.5 }},
            {{ -sqrt(3)/2.,	-0.5 }}
        };
    }
    else if(d == 3) {
        centers = {
            {{0,	1,	0}},
            {{sqrt(3)/2.,	-0.5, 0}},
            {{-sqrt(3)/2.,	-0.5, 0}},
            {{0,	0,	sqrt(3)}}
        };
    }

    int N_per_center = ceil(N/(double)centers.size());
    vector<PointSource_<d>> sources; 

    for(auto center : centers) {
        auto new_sources = generateUniformDistribution<d>(shift_dist * center, extent, 
            N_per_center);
        sources.insert(sources.end(), new_sources.begin(), new_sources.end()); 
    }

    return sources; 
}


// abs max, abs avg, rel max, rel avg errors
array<double, 4> abs_rel_pot_error(vector<double> fmm_pots, 
        vector<double> ref_pots) {

    vector<double> diffs(fmm_pots.size());

    std::transform (
        fmm_pots.begin(), fmm_pots.end(), ref_pots.begin(), 
        diffs.begin(), [](double a, double b) -> double 
        { return std::abs(((a-b))); } 
    );

    double max_abs_error = vec_max(diffs);
    double avg_abs_error = vec_avg(diffs);

    std::transform (
        diffs.begin(), diffs.end(), ref_pots.begin(), 
        diffs.begin(), [](double a, double b) -> double 
        { return a/std::abs(b); } 
    );

    double max_rel_error = vec_max(diffs);
    double avg_rel_error = vec_avg(diffs);

    return {max_abs_error, avg_abs_error, max_rel_error, avg_rel_error};
}

array<double, 4> abs_rel_frc_error(vector<Vector_<d>> fmm_forces, 
        vector<Vector_<d>> ref_forces) {

    vector<double> diffs(fmm_forces.size());

    std::transform (
        fmm_forces.begin(), fmm_forces.end(), ref_forces.begin(), 
        diffs.begin(), [](Vector_<d> a, Vector_<d> b) -> double 
        { return (a-b).norm(); } 
    );

    double max_abs_error = vec_max(diffs);
    double avg_abs_error = vec_avg(diffs);

    std::transform (
        diffs.begin(), diffs.end(), ref_forces.begin(), 
        diffs.begin(), [](double a, Vector_<d> b) -> double 
        { return a/b.norm(); } 
    );

    double max_rel_error = vec_max(diffs);
    double avg_rel_error = vec_avg(diffs);

    return {max_abs_error, avg_abs_error, max_rel_error, avg_rel_error};
}

std::array<vector<double>, 5> par_fmm_efficiencies(int N, double eps, size_t s, 
        bool uniform_tree, bool uniform_sources) {

    double extent = 10; 

    vector<PointSource_<d>> sources; 
    if(uniform_sources) {
        sources = generateUniformDistribution<d>(Vector_<d>(0), extent, N); 
    }
    else {
        sources = generateNonUniformDistribution<d>(extent, N); 
    }

    std::cout << "s = " << s << "\n";
    sourceLocationsToFile(sources, "./logs/source_locs.dat");

    vector<double> nprocs, tree_timings, pot_eval_timings, 
        frc_eval_timings, total_timings; 

    int min_procs = 1;
    int max_procs = 4; 
    int incr = 1; 

    for(int i = min_procs; i <= max_procs; i+=incr) {

        omp_set_num_threads(i);
        std::cout << omp_get_max_threads() << "\n";
        double t1 = omp_get_wtime();      

        AbstractFmmTree<d, field_type>* fmm_tree; 
        if(uniform_tree) {
            fmm_tree = new BalancedFmmTree<d, field_type>(sources, s, eps);  
        }
        else {
            fmm_tree = new AdaptiveFmmTree<d, field_type>(sources, s, eps);  
        }

        double t2 = omp_get_wtime();      
        auto potentials = fmm_tree->evaluateParticlePotentialEnergies();
        double t3 = omp_get_wtime();      
        auto forces = fmm_tree->evaluateParticleForces();
        double t4 = omp_get_wtime();      

//      auto ref_forces = particleForces<d, field_type>(sources, 0);
//      auto accuracies = abs_rel_frc_error(forces, ref_forces);  

//      printVec(accuracies, 4);  

        nprocs.push_back(i);
        tree_timings.push_back(t2-t1);  
        pot_eval_timings.push_back(t3-t2);  
        frc_eval_timings.push_back(t4-t3);  
        total_timings.push_back(t4-t1);

        delete fmm_tree; 
    }

    return {nprocs, tree_timings, pot_eval_timings, frc_eval_timings, total_timings}; 
}

void logEfficiencesToFile(array<vector<double>, 5> efficiencies, string head, string filename) {

    auto [nprocs, tree_timings, pot_eval_timings, frc_eval_timings, total_timings] 
        = efficiencies;


    std::ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific;
    ofile << head; 
    for(unsigned i = 0; i < tree_timings.size(); ++i) {
        ofile << nprocs[i] << "\t"
            << tree_timings[i] << "\t"
            << pot_eval_timings[i] << "\t" 
            << frc_eval_timings[i] << "\t" 
            << total_timings[i] << "\n";
    }
    ofile.close();

    printIterable(nprocs);
    printIterable(tree_timings);
    printIterable(pot_eval_timings);
    printIterable(frc_eval_timings);
    printIterable(total_timings);
}

int main(int argc, char *argv[]) {

    const size_t seed = 42; 
    srand(seed); 

    size_t N_exp = 5; 
    size_t N = pow(10, N_exp);
    const size_t items_per_leaf = 200;
    const double eps = 1E-4; 
    
    bool uniform_sources = true; 
    bool uniform_tree = true; 

    string head = "nprocs\ttreebuild.\tpotential\tforce\ttotal\n";
    string filename = string("logs/Eff_") + (field_type ? "gr" : "cl")
        + "_d" + to_string(d) + "_N" + to_string(N_exp) + "_"
        + (uniform_sources ? "us" : "ns") + "_" + (uniform_tree ? "bt" : "at"); 

    auto efficiencies = par_fmm_efficiencies(N, eps, items_per_leaf, 
            uniform_tree, uniform_sources); 
    logEfficiencesToFile(efficiencies, head, filename);

    uniform_sources = true; 
    uniform_tree = true; 

    head = "nprocs\ttreebuild.\tpotential\tforce\ttotal\n";
    filename = string("logs/Eff_") + (field_type ? "gr" : "cl")
        + "_d" + to_string(d) + "_N" + to_string(N_exp) + "_"
        + (uniform_sources ? "us" : "ns") + "_" + (uniform_tree ? "bt" : "at"); 

    efficiencies = par_fmm_efficiencies(N, eps, items_per_leaf, 
            uniform_tree, uniform_sources); 
    logEfficiencesToFile(efficiencies, head, filename);

    return 0;

}

