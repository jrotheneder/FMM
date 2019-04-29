#include <vector> 
#include <array> 
#include <complex> 
#include <iostream> 
#include <cstdlib> 
#include <cmath> 

#include "debugging.hpp" 
#include "quadtree.hpp" 

template<int d>
struct Vector {

    std::array<double, d> coords;

    std::array<double, d>& data() { return coords; }

    double& operator[](std::size_t index) { return coords[index]; }
    const double& operator[](std::size_t index) const { return coords[index]; }

    Vector operator+(const Vector& rhs) const {
        std::array<double, d> summed_coords;
        for(std::size_t i = 0; i < rhs.coords.size(); ++i) {
            summed_coords[i] = (*this)[i] + rhs[i];  
        }
        return Vector{summed_coords};
    }
    Vector operator*(const double s) const {
        std::array<double, d> scaled_coords = this->coords;
        for(double& c : scaled_coords) { c *= s; }
        return Vector{scaled_coords};
    }
    Vector operator-() const { return -1 * *this; }
    Vector operator-(const Vector& rhs) const { return *this + (-rhs); }
//  Vector& operator+=(const Vector& rhs) { return *this = *this + rhs; }
//  Vector& operator-=(const Vector& rhs) { return *this = *this - rhs; }

    friend Vector operator*(const double s, const Vector& rhs) { return rhs * s; }
    friend std::ostream& operator<<(std::ostream& o, const Vector& p) {
       for(std::size_t i = 0; i < d - 1; i++) { o << p[i] << ", "; }
       o << p[d-1];  
       return o; 
    }
};

using namespace std;

int main(int argc, char *argv[]) {
     
    size_t N = 5000;
    const size_t d = 3;
    vector<Vector<d>> pts;

    for(std::size_t i = 0; i < N; i++) {
        Vector<d> v;      
        for(std::size_t j = 0; j < d; ++j) {
            v[j] =  64 * ((double) rand() / (RAND_MAX));
        }
        pts.push_back(Vector<d>(v)); 
    }

//  print_vec(pts); 

    Quadtree<Vector<d>, bool, d>q(pts, 100);
    cout << "Quadtree height is " << q.root->height << endl;
    q.toFile("geometry.dat", "points.dat");

//  auto dirs =  Quadtree<Vector<3>, bool, 3>::getChildCenterDirections();
//  for(auto vec : dirs) {
//      cout << vec << "\n";
//  }
    
    
    
    

    return 0;

}

