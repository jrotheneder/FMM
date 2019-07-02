#include <vector> 
#include <stdexcept> 
#include <algorithm> 
#include <cmath> 

#include "WolframLibrary.h"
#include "../adaptive_fmm_tree.hpp" 

#ifndef FIELD_TYPE // True => Grav. field, false => Electrostatic
#error "Field type undefined."
#endif

#ifndef DIM
#error "Dimension undefined."
#endif

using namespace fmm;

EXTERN_C DLLEXPORT int FMMParticleForces(WolframLibraryData libData, mint Argc, 
        MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mint items_per_leaf = MArgument_getInteger(Args[2]); 
    const mreal accuracy_eps = MArgument_getReal(Args[3]); 
    const mreal force_smoothing_eps = MArgument_getReal(Args[4]); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_sources * DIM};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        Vec v; 
        for(int j = 0; j < DIM; ++j) {
            v[j] = positions_data[DIM*i + j];    
        }
        double q = charge_data[i];  
        sources[i] = Src{v,q}; 
    }

    AdaptiveFmmTree<DIM, FIELD_TYPE> fmm_tree(sources, items_per_leaf, 
            accuracy_eps, force_smoothing_eps);  

    std::vector<Vec> forces = fmm_tree.evaluateParticleForces();  

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        for(int j = 0; j < DIM; ++j) {
            out_data[DIM*i + j] = forces[i][j];     
        }
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int FMMParticlePotentialEnergies(WolframLibraryData libData, 
        mint Argc, MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mint items_per_leaf = MArgument_getInteger(Args[2]); 
    const mreal accuracy_eps = MArgument_getReal(Args[3]); 
    const mreal force_smoothing_eps = MArgument_getReal(Args[4]); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_sources};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        Vec v; 
        for(int j = 0; j < DIM; ++j) {
            v[j] = positions_data[DIM*i + j];    
        }
        double q = charge_data[i];  
        sources[i] = Src{v,q}; 
    }

    AdaptiveFmmTree<DIM, FIELD_TYPE> fmm_tree(sources, items_per_leaf, 
            accuracy_eps, force_smoothing_eps);  

    std::vector<double> potential_energies 
        = fmm_tree.evaluateParticlePotentialEnergies();  

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        out_data[i] = potential_energies[i];     
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}


EXTERN_C DLLEXPORT int DirectParticleForces(WolframLibraryData libData, mint Argc, 
        MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mreal force_smoothing_eps = MArgument_getReal(Args[2]); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_sources * DIM};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        Vec v; 
        for(int j = 0; j < DIM; ++j) {
            v[j] = positions_data[DIM*i + j];    
        }
        double q = charge_data[i];  
        sources[i] = Src{v,q}; 
    }

    std::vector<Vec> forces 
        = fields::particleForces<DIM, FIELD_TYPE>(sources, force_smoothing_eps);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        for(int j = 0; j < DIM; ++j) {
            out_data[DIM*i + j] = forces[i][j];     
        }
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int DirectParticlePotentialEnergies(WolframLibraryData libData, 
        mint Argc, MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mreal force_smoothing_eps = MArgument_getReal(Args[2]); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_sources};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        Vec v; 
        for(int j = 0; j < DIM; ++j) {
            v[j] = positions_data[DIM*i + j];    
        }
        double q = charge_data[i];  
        sources[i] = Src{v,q}; 
    }

    std::vector<double> potential_energies 
        = fields::particlePotentialEnergies<DIM, FIELD_TYPE>(sources, 
            force_smoothing_eps);  

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < n_sources; ++i) {
        out_data[i] = potential_energies[i];     
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}


EXTERN_C DLLEXPORT int FMMEvaluatePotentials(WolframLibraryData libData, 
        mint Argc, MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mint items_per_leaf = MArgument_getInteger(Args[2]); 
    const mreal accuracy_eps = MArgument_getReal(Args[3]); 
    const mreal force_smoothing_eps = MArgument_getReal(Args[4]); 

    MTensor flat_eval_points = MArgument_getMTensor(Args[5]); 
    mbool tree_to_file = MArgument_getBoolean(Args[6]); 

    mint n_eval_points = libData->MTensor_getDimensions(flat_eval_points)[0]/DIM;
    mreal* eval_points_data = libData->MTensor_getRealData(flat_eval_points); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_eval_points};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);
    std::vector<Vec> eval_points(n_eval_points);

    #pragma omp parallel 
    {
        #pragma omp for schedule(dynamic)
        for(int i = 0; i < n_sources; ++i) {
            Vec v; 
            for(int j = 0; j < DIM; ++j) {
                v[j] = positions_data[DIM*i + j];    
            }
            double q = charge_data[i];  
            sources[i] = Src{v,q}; 
        }

        #pragma omp for schedule(dynamic)
        for(int i = 0; i < n_eval_points; ++i) {
            Vec v; 
            for(int j = 0; j < DIM; ++j) {
                v[j] = eval_points_data[DIM*i + j];    
            }
            eval_points[i] = v;  
        }
    }
    
    AdaptiveFmmTree<DIM, FIELD_TYPE> fmm_tree(sources, items_per_leaf, 
            accuracy_eps, force_smoothing_eps);  

    if(tree_to_file) {
        fmm_tree.toFile(); 
    }

    
    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < n_eval_points; ++i) {
        out_data[i] = fmm_tree.evaluatePotential(eval_points[i]);     
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}


EXTERN_C DLLEXPORT int DirectEvaluatePotentials(WolframLibraryData libData, 
        mint Argc, MArgument *Args, MArgument Res) {

    using Vec = Vector_<DIM>;
    using Src = PointSource_<DIM>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/DIM;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mreal force_smoothing_eps = MArgument_getReal(Args[2]); 
    MTensor flat_eval_points = MArgument_getMTensor(Args[3]); 

    mint n_eval_points = libData->MTensor_getDimensions(flat_eval_points)[0]/DIM;
    mreal* eval_points_data = libData->MTensor_getRealData(flat_eval_points); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_eval_points};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources(n_sources);
    std::vector<Vec> eval_points(n_eval_points);

    #pragma omp parallel 
    {
        #pragma omp for schedule(dynamic)
        for(int i = 0; i < n_sources; ++i) {
            Vec v; 
            for(int j = 0; j < DIM; ++j) {
                v[j] = positions_data[DIM*i + j];    
            }
            double q = charge_data[i];  
            sources[i] = Src{v,q}; 
        }

        #pragma omp for schedule(dynamic)
        for(int i = 0; i < n_eval_points; ++i) {
            Vec v; 
            for(int j = 0; j < DIM; ++j) {
                v[j] = eval_points_data[DIM*i + j];    
            }
            eval_points[i] = v;  
        }
        
        #pragma omp for schedule(dynamic) 
        for(long int i = 0; i < n_eval_points; ++i) {
            out_data[i] = fields::potential<DIM, FIELD_TYPE, true>(
                sources, eval_points[i], force_smoothing_eps
                );
        }
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}
