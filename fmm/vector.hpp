#ifndef VECTOR_H
#define VECTOR_H

#include <array> 
#include <iostream> 
#include <numeric> 
#include <cassert> 
#include <cmath> 

namespace fmm {

template<std::size_t d>
struct Vector_ {
    
    std::array<double, d> coords;

    Vector_(std::array<double, d> coords = {}): coords(coords) {}
    Vector_(double fill_value): coords() { this->fill(fill_value); }

    double& operator[](std::size_t index) { return coords[index]; }
    const double& operator[](std::size_t index) const { return coords[index]; }

    Vector_ operator+(const Vector_& rhs) const {
        std::array<double, d> summed_coords;
        for(std::size_t i = 0; i < rhs.coords.size(); ++i) {
            summed_coords[i] = (*this)[i] + rhs[i];  
        }
        return Vector_{summed_coords};
    }

    Vector_ operator*(const double s) const {
        std::array<double, d> scaled_coords = this->coords;
        for(double& c : scaled_coords) { c *= s; }
        return Vector_{scaled_coords};
    }

    Vector_ operator-() const { return -1. * *this; }
    Vector_ operator-(const Vector_& rhs) const { return *this + (-rhs); }
    Vector_& operator*=(const double s) { return *this = *this * s; }
    Vector_& operator/=(const double s) { return *this *= (1/s); }
    Vector_& operator+=(const Vector_& rhs) { return *this = *this + rhs; }
    Vector_& operator-=(const Vector_& rhs) { return *this = *this - rhs; }

    friend Vector_ operator*(const double s, const Vector_& rhs) { return rhs * s; }
    friend std::ostream& operator<<(std::ostream& o, const Vector_& p) {
       for(std::size_t i = 0; i < d - 1; i++) { o << p[i] << ", "; }
       o << p[d-1];  
       return o; 
    }

    bool operator==(const Vector_& rhs) const { return this->coords == rhs.coords; }
    bool operator!=(const Vector_& rhs) const { return !(*this == rhs); }

    std::array<double, d>& data() { return coords; }
    const std::array<double, d>& data() const { return coords; }
    void fill(double val) { coords.fill(val); }

    double dot(const Vector_& rhs) const {
        return std::inner_product(coords.begin(), coords.end(), 
                rhs.data().begin(), 0.); 
    }

    double norm() const { 
        if constexpr(d == 1) {
            return std::abs(coords[0]); 
        }
        else if(d==2) {
            return std::hypot(coords[0], coords[1]);
        }
        else if(d == 3)  {
            return std::hypot(coords[0], coords[1], coords[2]);
        }

        return std::sqrt(this->dot(*this)); 
    } 

    double norm_sq() const {
        if constexpr(d == 1) {
            return coords[0] * coords[0]; 
        }
        else if(d==2) {
            return coords[0] * coords[0] + coords[1] * coords[1]; 
        }
        else if(d == 3)  {
            return coords[0] * coords[0] 
                + coords[1] * coords[1]
                + coords[2] * coords[3]; 
        }

        return this->dot(*this); 
    }

    // Returns vector of  (r, φ) of polar coordinates from vector (x, y) of
    // cartesian coordinates if d == 2, vector of (r, θ, φ) spherical 
    // from vector (x, y, z) of cartesian coordinates if d == 3 and is not 
    // implemented for other cases.
    Vector_ toSpherical() const {

        if constexpr (d == 2) {

            double r = this->norm();
            double phi = std::atan2(coords[1], coords[0]); 

            return Vector_{{r, phi}};
        }
        else if(d == 3) {

            double r = this->norm(); 
            double theta = std::atan2(std::hypot(coords[0], coords[1]), coords[2]); 
            double phi = std::atan2(coords[1], coords[0]); 

            return Vector_{{r, theta, phi}};
        }
        else {
            throw std::logic_error("toSpherical() implemented only "
                "for 2 & 3 dimensions.");
        }
    }

    // Returns vector of (x, y) cartesian coordinates from vector (r, φ) of 
    // polar coordinates if d == 2, vector of (x, y, z) cartesian 
    // from vector (r, θ, φ) of spherical coordinates if d == 3 and is not 
    // implemented for other cases.
    Vector_ toCartesian() const {

        if constexpr (d == 2) {

            double x = coords[0] * std::cos(coords[1]); 
            double y = coords[0] * std::sin(coords[1]); 

            return Vector_{{x, y}};
        }
        else if(d == 3) {

            double cos_theta = std::cos(coords[1]); 
            double sin_theta = std::sin(coords[1]); 
            double cos_phi = std::cos(coords[2]); 
            double sin_phi = std::sin(coords[2]); 

            double x = coords[0] * sin_theta * cos_phi; 
            double y = coords[0] * sin_theta * sin_phi; 
            double z = coords[0] * cos_theta; 

            return Vector_{{x, y, z}};
        }
        else {
            throw std::logic_error("toCartesian() implemented only "
                "for 2 & 3 dimensions.");
        }
    }

    // Transforms components (v_r, v_t, v_p) of a vector (*this) given w.r.t. 
    // the basis vectors
    // \hat r     = Cos(phi) Sin(theta)     Sin(theta) Sin(phi) 	Cos(theta)
    // \hat theta = Cos(theta) Cos(phi)	    Cos(theta) Sin(phi)	    -Sin(theta)
    // \hat phi   = -Sin(phi)               Cos(phi)	            0
    // into components (x,y,z) of the cartesian standard basis
    Vector_ toCartesianBasis(const double theta, const double phi) const {

        if constexpr(d == 3) {

            const double sin_theta = std::sin(theta);  
            const double cos_theta = std::cos(theta);  
            const double sin_phi = std::sin(phi);  
            const double cos_phi = std::cos(phi);  

            return Vector_{{ 
                sin_theta * cos_phi * coords[0] + cos_theta * cos_phi * coords[1] 
                    - sin_phi * coords[2],
                sin_theta * sin_phi * coords[0] + cos_theta * sin_phi * coords[1] 
                    + cos_phi * coords[2],
                cos_theta * coords[0] - sin_theta * coords[1] 
            }};
        }
        else {
            throw std::logic_error("toCartesianBasis() implemented only "
                "for 3 dimensions.");
        }
    }
};

template<std::size_t d>
struct PointSource_ {

    Vector_<d> position; 
    double q; 

    PointSource_(Vector_<d> position = {}, double q = 0): 
        position(position), q(q) {}
    PointSource_(std::array<double, d> position, double q): 
        position(position), q(q) {}

    double& operator[](std::size_t index) { return position[index]; }
    const double& operator[](std::size_t index) const { return position[index]; }

    bool operator==(const PointSource_& rhs) const { 
        return this->q == rhs.q && this->position == rhs.position;
    }

    double sourceStrength() const { return q; }

    friend std::ostream& operator<<(std::ostream& o, const PointSource_& p) {
        o << p.position << ", q = " << p.q; 
        return o; 
    }
};

} //namespace fmm

#endif
