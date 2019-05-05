#ifndef VECTOR_H
#define VECTOR_H

#include <array> 
#include <iostream> 
#include <numeric> 
#include <cassert> 

template<std::size_t d>
struct Vector {
    std::array<double, d> coords;

    Vector(std::array<double, d> coords = {}): coords(coords) {}
    Vector(std::initializer_list<double> il) {
        assert(il.size() <= d);
        coords = std::array<double, d>{il};
    }

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

    std::array<double, d>& data() { return coords; }
    const std::array<double, d>& data() const { return coords; }
    void fill(double val) { coords.fill(val); }

    double dot(const Vector& rhs) const {
        return std::inner_product(this->data().begin(), this->data().end(), 
            rhs.data().begin(), 0.); 
    }
    double norm() const { return sqrt(this->dot(*this)); }
};

template<std::size_t d>
struct PointCharge {

    Vector<d> position; 
    double q; 

    PointCharge(Vector<d> position = {}, double q = 0): position(position), q(q) {}
    PointCharge(std::array<double, d> position, double q): position(position), q(q) {}

    double sourceStrength() const { return q; }

    friend std::ostream& operator<<(std::ostream& o, const PointCharge& p) {
        o << static_cast<const Vector<d>&>(p) << ", q = " << p.q; 
        return o; 
    }
};

#endif
