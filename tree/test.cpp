#include <vector> 
#include <complex> 
#include <iostream> 
#include <cstdlib> 
#include <cmath> 

#include "quadtree.hpp" 

typedef std::complex<double> Point;

using namespace std;

int main(int argc, char *argv[]) {
     
    size_t N = 1000;
    vector<complex<double>> points;
    for(std::size_t i = 0; i < N; i++) {

        double x = 64 * ((double) rand() / (RAND_MAX));
        double y = 64 * ((double) rand() / (RAND_MAX));
        points.push_back(std::complex<double>(x,y)); 
    }

    Quadtree<Point, bool, 2>q(points, 1);
    cout << "Quadtree height is " << q.root->height << endl;
    q.toFile("geometry.dat", "points.dat");

    return 0;

}

