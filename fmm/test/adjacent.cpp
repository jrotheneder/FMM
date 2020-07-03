#define protected public // this is risky to say the least, but it works in this case
#include "../adaptive_fmm_tree.hpp"

using namespace fmm; 

int main(int argc, char *argv[]) {

    const unsigned d = 3; 

    using Node = AdaptiveFmmTree<d>::FmmNode;
    using Vec = Vector_<d>;

    const Vec center{};
    const Vec ones(1); 
    const Vec shift_y{{0, 1 + 5 * std::numeric_limits<double>::epsilon()}};

    Node root = Node(center, 1, 0, nullptr, 5);  
    Node child1 = Node(center + ones, 1, 0, &root, 5);  
    Node child2 = Node(center - ones, 1, 0, &root, 5);  
    Node child3 = Node(center - shift_y, 1, 0, &root, 5);  

    std::cout << "root.adjacent(root) = "     << root.adjacent(&root) << "\n";
    std::cout << "root.adjacent(child1) = "   << root.adjacent(&child1) << "\n";
    std::cout << "root.adjacent(child2) = "   << root.adjacent(&child2) << "\n";
    std::cout << "child1.adjacent(root) = "   << child1.adjacent(&root) << "\n";
    std::cout << "child1.adjacent(child1) = " << child1.adjacent(&child1) << "\n";
    std::cout << "child1.adjacent(child2) = " << child1.adjacent(&child2) << "\n";
    std::cout << "child1.adjacent(child3) = " << child1.adjacent(&child3) << "\n";
    std::cout << "child3.adjacent(root) = "   << child3.adjacent(&root) << "\n";
    std::cout << "child3.adjacent(child1) = " << child3.adjacent(&child1) << "\n";
    std::cout << "child3.adjacent(child2) = " << child3.adjacent(&child2) << "\n";
    std::cout << "child3.adjacent(child3) = " << child3.adjacent(&child3) << "\n";

}
