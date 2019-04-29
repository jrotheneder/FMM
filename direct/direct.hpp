#include <vector> 
#include <functional> 

// Generic class for computing pairwise interactions between N sources,
// e.g. potentials or forces in a system of charges in electrostatics
template<typename SourceType, typename InteractionType>
class Interaction {

public:
    const std::size_t N; 
    const std::vector<SourceType>& sources; 
    std::vector<InteractionType> interactions; 
    const std::function <InteractionType (const SourceType&, 
        const SourceType&)> interactionFunction;
    
    Interaction(const std::vector<SourceType>& sources, const std::function 
        <InteractionType (const SourceType&, const SourceType&)>& 
        interactionFunction): N(sources.size()), sources(sources), 
        interactions(std::vector<InteractionType>(N)),
        interactionFunction(interactionFunction) {
    }  

    // Computes interactions with an arbitrary (i.e. nonsymmetric) kernel,
    // ignores self-interaction. Relies on InteractionType being default 
    // initialized to 'zero'
    void Compute() {

        for(std::size_t i = 0; i < N; i++) {

            SourceType src = sources[i];
            InteractionType ia = InteractionType(); // Default must be 'zero' el.

            for(std::size_t j = 0; j < i; j++) {
                ia += interactionFunction(src, sources[j]);
            }
            for(std::size_t j = i+1; j < N; j++) {
                ia += interactionFunction(src, sources[j]);
            }

            interactions[i] = ia;
        }
    }

    // Computes interactions with a symmetric kernel, i.e. if 
    // interactionFunction(a,b) == interactionFunction(b,a) for all a,b. 
    // Ignores self-interaction, relies on InteractionType being default 
    // initialized to 'zero'
    void ComputeSymmetric() {

        for(std::size_t i = 0; i < N; i++) {

            SourceType src = sources[i];
            InteractionType ia = InteractionType(); // Default must be 'zero' el.

            for(std::size_t j = 0; j < i; j++) {
                InteractionType new_ia = interactionFunction(src, sources[j]);
                ia += new_ia;
                interactions[j] += new_ia;
            }

            interactions[i] = ia;
        }
    }

    // Computes interactions with an antisymmetric kernel, i.e. if 
    // interactionFunction(a,b) == -interactionFunction(b,a) for all a,b. 
    // Ignores self-interaction, relies on InteractionType being default 
    // initialized to 'zero'
    void ComputeAntisymmetric() {

        for(std::size_t i = 0; i < N; i++) {

            SourceType src = sources[i];
            InteractionType ia = InteractionType(); // Default must be 'zero' el.

            for(std::size_t j = 0; j < i; j++) {
                InteractionType new_ia = interactionFunction(src, sources[j]);
                ia += new_ia;
                interactions[j] -= new_ia;
            }

            interactions[i] = ia;
        }
    }

};
