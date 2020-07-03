#include <vector> 

#include "WolframLibrary.h"
#include "../adaptive_fmm_tree.hpp" 

using namespace fmm;

EXTERN_C DLLEXPORT int FMM2D_Forces(WolframLibraryData libData, mint Argc, 
        MArgument *Args, MArgument Res) {

    const int d = 2;

    using Vec = Vector_<d>;
    using Src = PointSource_<d>;
 
    MTensor positions = MArgument_getMTensor(Args[0]); 
    MTensor charges = MArgument_getMTensor(Args[1]); 
    mint count = libData->MTensor_getDimensions(positions)[0]; 
    mint n_sources = count/d;

    mreal* positions_data = libData->MTensor_getRealData(positions); 
    mreal* charge_data = libData->MTensor_getRealData(charges); 

    const mint items_per_leaf = MArgument_getInteger(Args[2]); 
    const mreal eps = MArgument_getReal(Args[3]); 

    MTensor out; 
    mint out_rank = 1;
    mint out_dims[1] {n_sources * d};
    mint err = libData->MTensor_new(MType_Real, out_rank, out_dims, &out);

    if(err) {
        throw std::runtime_error("Error in libData->MTensor_new()");  
    }

    mreal* out_data = libData->MTensor_getRealData(out);

    std::vector<Src> sources;

    for(int i = 0; i < n_sources; ++i) {
        Vec v; 
        for(int j = 0; j < d; ++j) {
            v[j] = positions_data[d*i + j];    
        }
        double q = charge_data[i];  
        sources.push_back({v,q}); 
    }

    AdaptiveFmmTree<d> fmm_tree(sources, items_per_leaf, eps);  
    std::vector<Vec> forces = fmm_tree.evaluateParticleForces();  

    for(int i = 0; i < n_sources; ++i) {
        for(int j = 0; j < d; ++j) {
            out_data[d*i + j] = forces[i][j];     
        }
    }
    
    MArgument_setMTensor(Res,out);
    return LIBRARY_NO_ERROR;
}
