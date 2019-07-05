#include <filesystem> 

#include <gsl/gsl_odeiv2.h> 
#include <gsl/gsl_errno.h>

#include "../adaptive_fmm_tree.hpp"

#include <valgrind/callgrind.h> 

using namespace fmm; 
using namespace std; 

// Writes first half of vectors in positions together with charges to output
// vector sources
template<std::size_t d>
void arrayToSources(std::vector<PointSource_<d>>& sources, 
        const double positions[], const double charges[]) {

    for(std::size_t i = 0; i < sources.size(); ++i) {
        Vector_<d> v; 
        for(std::size_t j = 0; j < d; ++j) {
            v[j] = positions[i*d + j];   
        }
        sources[i] = PointSource_<d>{v, charges[i]};  
    }
}

// Writes source positions into the first half of the positions vector and
// charges into the corresponding charge vector
template<std::size_t d> 
void sourcesToArray(const std::vector<PointSource_<d>>& sources, 
        double positions[], double charges[]) {

    for(std::size_t i = 0; i < sources.size(); ++i) {
        for(std::size_t j = 0; j < d; ++j) {
            positions[i*d + j] = sources[i][j];   
        }

        // technically superflous assuming that the tree does not reorder its
        // input.
        charges[i] = sources[i].sourceStrength();   
    }
}

struct deriv_params {

    double* charges; 
    double * masses; 

    std::size_t N; 
    std::size_t items_per_leaf; 

    double fmm_accuracy_eps; 
    double force_smoothing_eps; 
};

template<std::size_t d>
int derivative_func(double t, const double y[], double f[], void *params) {

    (void)(t); // avoid unused parameter warning

    deriv_params* info = (deriv_params*) params;

    double* charges = info->charges; 
    double* masses = info->masses; 
    std::size_t N = info->N; 

    std::vector<PointSource_<d>> sources(N);
    arrayToSources(sources, y, charges);  

    AdaptiveFmmTree<d> fmm_tree(sources, info->items_per_leaf, 
            info->fmm_accuracy_eps, info->force_smoothing_eps); 

    std::vector<Vector_<d>> forces = fmm_tree.evaluateParticleForces(); 

    for(std::size_t i = 0; i < N; ++i) {
        for(std::size_t j = 0; j < d; ++j) {
            f[i*d + j] = y[N*d + i*d + j]; 
            f[N*d + i*d + j] = forces[i][j]/masses[i];  
        }
    }

    std::cout << "t = " << t << "\n";

    return GSL_SUCCESS;
}

int main(int argc, char *argv[]) {

    string input_folder = "./input/"; 
    string sources_file = argc > 1 ? argv[1] : input_folder + "sources.dat"; 
    string velocities_file = argc > 2 ? argv[2] : input_folder + "velocities.dat"; 
    string masses_file = argc > 3 ? argv[3] : input_folder + "masses.dat"; 

    if(!filesystem::exists(sources_file) 
        || !filesystem::exists(velocities_file)
        || !filesystem::exists(masses_file)) {

        throw runtime_error("File not found."); 
    }

    string trajectory_filename = "logs/trajectory.dat"; 

    ofstream ofile; 
    ofile.open(trajectory_filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) 
        << std::scientific;


    const unsigned d = 2;       // Dimension 
    using Src = PointSource_<d>; 

    double t = 0.0, t1 = 10;  // start and end times
    const unsigned n_saves = 100; 

    unsigned N; 
    const unsigned items_per_leaf = 90; 
    const double fmm_accuracy_eps = 1E-3; 
    const double odeint_eps = 1E-3; 
    const double force_smoothing_eps = 1E-6; 

    double* positions;  
    double* charges;  
    double* masses;  

    {
        vector<Src> sources = readSourcesFromFile<d>(sources_file);  
        N = sources.size(); 

        positions = new double[2*d*N]; // contains positions and velocities
        charges = new double[N];   
        masses = new double[N]; 

        sourcesToArray(sources, positions, charges); 
        vectorFromFile(velocities_file, positions + d*N);  
        vectorFromFile(masses_file, masses);  
    }

    deriv_params info {charges, masses, N, items_per_leaf, fmm_accuracy_eps, 
        force_smoothing_eps};

    CALLGRIND_TOGGLE_COLLECT; 
    // https://www.gnu.org/software/gsl/doc/html/ode-initval.html
    gsl_odeiv2_system sys = {derivative_func<d>, NULL, 2*d*N, &info};
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, 
        gsl_odeiv2_step_rk4, 1E-6, odeint_eps, odeint_eps);

    std::cout << "\n";
    for (unsigned i = 1; i <= n_saves; i++) {

        double ti = i * t1 / n_saves; 
        int status = gsl_odeiv2_driver_apply (driver, &t, ti, positions);

        if (status != GSL_SUCCESS) {
            cout << "error, return value= " <<  status << "\n";
            break;
        }

        ofile << ti << "\n"; 
        for(unsigned j = 0; j < 2*d*N; ++j) { ofile << positions[j] << "\t";  }
        ofile << "\n"; 

        std::cout << "Iteration " << i << ", t = " << ti << "\n";

    }

    gsl_odeiv2_driver_free(driver);
    delete[] positions; 
    delete[] charges; 
    delete[] masses; 

    return 0;
}

