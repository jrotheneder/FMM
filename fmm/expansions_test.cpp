#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 

#include "debugging.hpp" 
#include "vector.hpp" 

//#include "point_orthtree.hpp" 
#include "multipole_expansion.hpp" 
#include "local_expansion.hpp" 
#include "../direct/direct.hpp" 

using namespace std;
using namespace fmm; 

int main(int argc, char *argv[]) {

    const size_t d = 2;

    using Vec = Vector<d>;
    using Src = PointCharge<d>;
    using ME = fmm::MultipoleExpansion<Vec, Src, d>;
    using LE = fmm::LocalExpansion<Vec, Src, d>;
     
    size_t N = 5000;
    const double extent = 16;

    const double eps = 1E-10; 
    const size_t order = ceil(log(1/eps) / log(2)); 
    const size_t seed = 42; 
    srand(seed); 

    Vec unit_1{}; unit_1[0] = 1;
    Vec unit_2{}; unit_2[1] = 1;
    const Vec ones(1); 

    const Vec center1 = - 2 * extent * unit_1; // Multipole expansion 1
    const Vec center2 = 2 * extent * unit_1;  // Multipole expansion 2
    const Vec se_center = 2 * extent * unit_2;  // center of shifted multipole expansion
    const Vec le_center{}; // Origin, center of local expansion
    const Vec sle_center = -0.5 * extent * ones; // Center of shifted local expansion

    const Vec me_eval_point(6 * extent); 
    const Vec le_eval_point = extent/2 * ones; 

    vector<Src> sources1, sources2;

    for(size_t i = 0; i < N; i++) {
        Vec v1, v2;      

        for(size_t j = 0; j < d; ++j) {
            v1[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
            v2[j] =  extent * ((double) rand() / (RAND_MAX)) - extent/2;
        }

        double q1 = (double) rand() / (RAND_MAX) * (rand() % 2 ? 1 : -1);
        double q2 = (double) rand() / (RAND_MAX) * (rand() % 2 ? 1 : -1);

        v1 += center1; 
        v2 += center2; 

        Src src1 {v1, q1};
        Src src2 {v2, q2};

        sources1.push_back(src1); 
        sources2.push_back(src2); 
    }

    // Compute reference values directly: 

    const double me_pot_ref = direct::evaluateInteraction<Vec, Src, 
        double>(sources1, me_eval_point, direct::Electrostatic_Potential<Vec, 
        Src, d>) 
        + direct::evaluateInteraction<Vec, Src, double>(sources2, 
        me_eval_point, direct::Electrostatic_Potential<Vec, Src, d>);

    const double le_pot_ref = direct::evaluateInteraction<Vec, Src, 
        double>(sources1, le_eval_point, direct::Electrostatic_Potential<Vec, 
        Src, d>) 
        + direct::evaluateInteraction<Vec, Src, double>(sources2, 
        le_eval_point, direct::Electrostatic_Potential<Vec, Src, d>);

    const Vec me_force_ref = direct::evaluateInteraction<Vec, 
        Src, Vec>(sources1, me_eval_point, 
        direct::Electrostatic_Force<Vec, Src, d>)
        + direct::evaluateInteraction<Vec, Src, 
        Vec>(sources2, me_eval_point, direct::Electrostatic_Force<Vec, 
        Src, d>);

    const Vec le_force_ref = direct::evaluateInteraction<Vec, 
        Src, Vec>(sources1, le_eval_point, direct::Electrostatic_Force<Vec, Src, d>)
        + direct::evaluateInteraction<Vec, Src, Vec>(sources2, le_eval_point, 
        direct::Electrostatic_Force<Vec, Src, d>);

    // Set up multipole expansion from sources:
    ME me1(center1, order, sources1); 
    ME me2(center2, order, sources2); 

    // ... and multipole expansion by shift: 
    std::vector<ME*> vme{&me1, &me2};
    ME se(se_center, vme); 

    // Set up local expansions from multipole expansions: 
    LE le(le_center, vme); 
    // And a shifted local expansion: 
    std::vector<LE*> vle{&le};
    LE sle(sle_center, vle); 

    double me_pot_error = abs(me1.evaluatePotential(me_eval_point) 
        + me2.evaluatePotential(me_eval_point) - me_pot_ref); 
    double se_pot_error = abs(se.evaluatePotential(me_eval_point) - me_pot_ref); 
    double le_pot_error = abs(le.evaluatePotential(le_eval_point) - le_pot_ref); 
    double sle_pot_error = abs(sle.evaluatePotential(le_eval_point) - le_pot_ref); 

    double me_force_error = (me1.evaluateForcefield(me_eval_point) 
        + me2.evaluateForcefield(me_eval_point) - me_force_ref).norm(); 
    double se_force_error = (se.evaluateForcefield(me_eval_point) - me_force_ref).norm(); 
    double le_force_error = (le.evaluateForcefield(le_eval_point) - le_force_ref).norm(); 
    double sle_force_error = (sle.evaluateForcefield(le_eval_point) 
            - le_force_ref).norm(); 

    cout << "Order is " << order << endl;
    cout << "Potential errors:" << endl 
        << "Multipole: " << me_pot_error << endl 
        << "Multipole (shifted): " << se_pot_error << endl 
        << "Local:" << le_pot_error << endl
        << "Local (shifted):" << sle_pot_error << endl << endl
        << "Force errors: " << endl
        << "Multipole: " << me_force_error << endl 
        << "Multipole (shifted): " << se_force_error << endl
        << "Local: " << le_force_error << endl
        << "Local (shifted): " << sle_force_error << endl; 

//  cout << "diff local vs local (shifted) = " << 
//      abs(le.evaluatePotential(le_eval_point) 
//      - sle.evaluatePotential(le_eval_point)) 
//      << endl;

//  cout << "Coeff le: "; 
//  for(auto c : le.coefficients) { 
//      cout << " " << c.real() << " + " << c.imag() << "i\n"; 
//  }
//  cout << endl << "Coeff sle: "; 
//  for(auto c : sle.coefficients) { 
//      cout << " " << c.real() << " + " << c.imag() << "i\n"; 
//  }

    bool failure = 
        me_pot_error >= eps 
        || se_pot_error >= eps 
        || le_pot_error >= eps 
        || sle_pot_error >= eps 
        || me_force_error >= eps 
        || se_force_error >= eps
        || le_force_error >= eps
        || sle_force_error >= eps; 

    if(failure) { std::cout << "Test failed in expansions_test.cpp" << std::endl; }

    return failure;
}

