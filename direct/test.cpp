#include <vector> 
#include <complex> 
#include <iostream> 
#include <cstdlib> 
#include <cmath> 

#include "direct.hpp" 
#include "debugging.hpp" 

using namespace std;

struct Vector2D {
    double x;
    double y;

    Vector2D(): x(0), y(0) {}
    Vector2D(double x, double y): x(x), y(y) {}

    Vector2D& operator+=(const Vector2D& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }
    Vector2D& operator-=(const Vector2D& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
    }
    Vector2D operator-() { return Vector2D(-this->x, -this->y); }
    friend Vector2D operator+(Vector2D lhs, const Vector2D &rhs) {
        return lhs += rhs;
    };
    friend ostream& operator<<(ostream& os, const Vector2D& v) {
        os << v.x << " " << v.y << " "; 
        return os;
    }
};

struct PointMass2D: Vector2D {

    double m; 
    PointMass2D(double x, double y, double m): Vector2D(x,y), m(m) {}

    friend ostream& operator<<(ostream& os, const PointMass2D& p) {
        os << static_cast<const Vector2D&>(p) << p.m; 
        return os;
    } 
};


int main(int argc, char *argv[]) {
     
    size_t N = 1000;
    vector<PointMass2D> points;
    for(std::size_t i = 0; i < N; i++) {

        double x = 64 * ((double) rand() / (RAND_MAX));
        double y = 64 * ((double) rand() / (RAND_MAX));
        double m =      ((double) rand() / (RAND_MAX));
        points.push_back(PointMass2D(x,y,m)); 
    }

    Interaction<PointMass2D, Vector2D> forceInteraction(points, 
        [] (const PointMass2D& s1, const PointMass2D& s2) -> Vector2D {

            double dx = s1.x - s2.x; 
            double dy = s1.y - s2.y; 
            double rsq = dx*dx + dy*dy;

            double mult = s1.m * s2.m; 
            double fx = mult * dx / rsq; 
            double fy = mult * dy / rsq; 

            return Vector2D(fx, fy);

        });


    Interaction<PointMass2D, double> potentialInteraction(points, 
        [] (const PointMass2D &s1, const PointMass2D& s2) -> double {

            double dx = s1.x - s2.x; 
            double dy = s1.y - s2.y; 
            double rsq = dx*dx + dy*dy;

            double phi = -0.5 * s1.m * s2.m * log(rsq); 
            return phi;

        });

    forceInteraction.ComputeAntisymmetric();
    potentialInteraction.ComputeSymmetric();

    VectorToFile(points, "points.dat"); 
    VectorToFile(forceInteraction.interactions, "forces.dat"); 
    VectorToFile(potentialInteraction.interactions, "potentials.dat"); 


}

