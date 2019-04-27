#include <vector> 
#include <complex> 
#include <iostream> 
#include <cstdlib> 
#include <cmath> 

#include "quadtree.hpp" 

using namespace std;

int main(int argc, char *argv[]) {
     
    size_t N = 1000;
    vector<complex<double>> points;
    for(std::size_t i = 0; i < N; i++) {

        double x = 64 * ((double) rand() / (RAND_MAX));
        double y = 64 * ((double) rand() / (RAND_MAX));
        points.push_back(std::complex<double>(x,y)); 
    }


    Quadtree q(points, 10);
    cout << "Quadtree height is " << q.height << endl;
    q.toFile("geometry.dat", "points.dat");

//  cout << "Diff to next before " << x_min << " is " << 
//      std::nextafter(x_min, -HUGE_VAL) - x_min << ", diff to next after " 
//      << x_max << " is " << std::nextafter(x_max, HUGE_VAL) - x_max << endl;

//  cout << "Comparisons: std::nextafter(x_min, -HUGE_VAL) < x_min: "  << 
//      (std::nextafter(x_min, -HUGE_VAL) < x_min) << 
//      ", std::nextafter(x_max, HUGE_VAL) > x_max: "
//      << (std::nextafter(x_max, HUGE_VAL) > x_max) << endl;

}

