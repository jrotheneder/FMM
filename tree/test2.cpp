#include <iostream> 
#include <utility> 
#include <stack> 
#include <tuple> 

template<int d>
struct Point {
    double coords[d];
    double m;
//  Point(double coords[d]) {
//      for(int i = 0; i < d; i++) {
//          this->coords[i] = coords[i];
//      }
//  }
};


int main(int argc, char *argv[]) {

//  Point<6> p = {{1,2,3,4}, 6}; 
//  std::cout << p.coords[4] << " " << p.m << std::endl;

    std::stack<std::pair<int, bool>> st({{3, false}});

    st.top().second = true;
    std::cout << st.top().second;
}

