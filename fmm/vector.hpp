#ifndef VECTOR_H
#define VECTOR_H

#include <array> 
#include <iostream> 
#include <numeric> 
#include <cassert> 

template<std::size_t d>
struct Vector {
    
    using NumType = double; 

    std::array<NumType, d> coords;

    Vector(std::array<NumType, d> coords = {}): coords(coords) {}
    Vector(NumType fill_value): coords() { this->fill(fill_value); }

    NumType& operator[](std::size_t index) { return coords[index]; }
    const NumType& operator[](std::size_t index) const { return coords[index]; }

    Vector operator+(const Vector& rhs) const {
        std::array<NumType, d> summed_coords;
        for(std::size_t i = 0; i < rhs.coords.size(); ++i) {
            summed_coords[i] = (*this)[i] + rhs[i];  
        }
        return Vector{summed_coords};
    }

    Vector operator*(const NumType s) const {
        std::array<NumType, d> scaled_coords = this->coords;
        for(NumType& c : scaled_coords) { c *= s; }
        return Vector{scaled_coords};
    }

    Vector operator-() const { return -1 * *this; }
    Vector operator-(const Vector& rhs) const { return *this + (-rhs); }
    Vector& operator+=(const Vector& rhs) { return *this = *this + rhs; }
    Vector& operator-=(const Vector& rhs) { return *this = *this - rhs; }

    friend Vector operator*(const NumType s, const Vector& rhs) { return rhs * s; }
    friend std::ostream& operator<<(std::ostream& o, const Vector& p) {
       for(std::size_t i = 0; i < d - 1; i++) { o << p[i] << ", "; }
       o << p[d-1];  
       return o; 
    }

    bool operator==(const Vector& rhs) const { return this->coords == rhs.coords; }

    std::array<NumType, d>& data() { return coords; }
    const std::array<NumType, d>& data() const { return coords; }
    void fill(NumType val) { coords.fill(val); }

    NumType dot(const Vector& rhs) const {
        return std::inner_product(coords.begin(), coords.end(), rhs.data().begin(), 0.); 
    }

    NumType norm() const { //TODO take care with extended precision...
        if constexpr(d == 1) {
            return std::abs(coords[0]); 
        }
        else if(d==2) {
            return std::hypot(coords[0], coords[1]);
        }
        else if(d == 3)  {
            return std::hypot(coords[0], coords[1], coords[2]);
        }
        else {
            return sqrt(this->dot(*this)); 
        }
    } 

    // Generalized spherical coordinates: returns vector of (r, φ) if d == 2, 
    // (r, θ, φ) if d == 3 and is not implemented for other cases.
    Vector toSpherical() const {

        if constexpr (d == 2) {

            NumType r = this->norm();
            NumType phi = std::atan2(coords[1], coords[0]); 

            return Vector{{r, phi}};
        }
        else if(d == 3) {

            NumType r = this->norm(); 
            NumType theta = std::atan2(std::hypot(coords[0], coords[1]), coords[2]); 
            NumType phi = std::atan2(coords[1], coords[0]); 

            return Vector{{r, theta, phi}};
        }
        else {
            throw std::logic_error("toSpherical() implemented only "
                "for 2 & 3 dimensions.");
        }
    }

    // Returns vector of (x, y) from vector (r, φ) if d == 2, vector of (x, y, z) 
    // from vector (r, θ, φ) if d == 3 and is not implemented for other cases.
    Vector toCartesian() const {

        if constexpr (d == 2) {

            NumType x = coords[0] * std::cos(coords[1]); 
            NumType y = coords[0] * std::sin(coords[1]); 

            return Vector{{x, y}};
        }
        else if(d == 3) {

            NumType cos_theta = std::cos(coords[1]); 
            NumType sin_theta = std::sin(coords[1]); 
            NumType cos_phi = std::cos(coords[2]); 
            NumType sin_phi = std::sin(coords[2]); 

            NumType x = coords[0] * sin_theta * cos_phi; 
            NumType y = coords[0] * sin_theta * sin_phi; 
            NumType z = coords[0] * cos_theta; 

            return Vector{{x, y, z}};
        }
        else {
            throw std::logic_error("toCartesian() implemented only "
                "for 2 & 3 dimensions.");
        }
    }
};

template<std::size_t d>
struct PointCharge {

    using NumType = double; 

    Vector<d> position; 
    NumType q; 

    PointCharge(Vector<d> position = {}, NumType q = 0): position(position), q(q) {}
    PointCharge(std::array<NumType, d> position, NumType q): position(position), q(q) {}

    bool operator==(const PointCharge& rhs) const { 
        return this->q == rhs.q && this->position == rhs.position;
    }

    NumType sourceStrength() const { return q; }

    friend std::ostream& operator<<(std::ostream& o, const PointCharge& p) {
        o << p.position << ", q = " << p.q; 
        return o; 
    }
};

#endif
