#include <iostream> 
#include <utility> 
#include <stack> 
#include <tuple> 
#include <array> 
#include <vector> 

using namespace std;

template<int d>
struct Point {

    double coords[d];
    double m;

    double& operator[](std::size_t index) { return coords[idx]; }
    const double& operator[](std::size_t index) const { return coords[idx]; }
};


int main(int argc, char *argv[]) {

//  Point<6> p = {{1,2,3,4}, 6}; 
//  std::cout << p.coords[4] << " " << p.m << std::endl;

//  array<Point<2>*, 10> arr {0}; 
//  for(auto pp : arr) {
//      cout << (pp == nullptr) << std::endl;
//  }
    
    std::array<double, 5> arr;
    arr.fill(1); 
    for(auto x : arr) {
        cout << x++ << " ";
    }
    
    for(auto x : arr) {
        cout << x++ << " ";
    }

}

