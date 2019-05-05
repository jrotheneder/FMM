#include "abstract_orthtree.hpp"

namespace fmm {
namespace node {

template<typename Vector, typename Source, std::size_t d>
struct FmmNode: AbstractOrthtree<Vector, d>::Node {

protected:
    using Super = typename AbstractOrthtree<Vector, d>::Node; 

public:
    struct MultipoleExpansion;
    struct LocalExpansion;

    MultipoleExpansion multipole_expansion; 
    LocalExpansion local_expansion;

    std::vector<FmmNode*> interaction_list;
    std::vector<FmmNode*> near_neighbours;

    FmmNode(Vector center, double box_length, std::size_t depth, 
        FmmNode * parent = nullptr): Super(center, box_length, depth, parent) {};

    virtual ~FmmNode() {};
};

template<typename Vector, typename Source, std::size_t d>
struct FmmLeaf: FmmNode<Vector, Source, d> {

    using Super = FmmNode<Vector, Source, d>;
    using BaseNode = typename Super::Node; 

    std::vector<Source> * sources;

    FmmLeaf(Vector center, double box_length, std::size_t depth, 
        Super * parent = nullptr): Super(center, box_length, depth, parent),
        sources() {};

    virtual ~FmmLeaf() { delete sources; }
};

template<typename Vector, typename Source, std::size_t d>
struct FmmNode<Vector, Source, d>::FmmNode::MultipoleExpansion {

    std::vector<double> coefficients; 

    MultipoleExpansion(std::vector<Source> sources = {}) {} //TODO

    std::vector<double> shift(Vector v) { return {}; }
};

template<typename Vector, typename Source, std::size_t d>
struct FmmNode<Vector, Source, d>::FmmNode::LocalExpansion {

    using MultipoleExpansion = 
        FmmNode<Vector, Source, d>::MultipoleExpansion;

    std::vector<double> coefficients; 

    LocalExpansion(MultipoleExpansion me = {}) {} // TODO 
    std::vector<double> shift(Vector v);
};

} // namespace node
} // namespace fmm
