#ifndef _VECTOR_3D_HPP
#define _VECTOR_3D_HPP

#include <iostream>
#include <cmath>

class vector3D {
private:
    double _crdVec[3] = {0, 0, 0};
public:
    constexpr vector3D() = default;
    constexpr vector3D(const vector3D&) = default;
    
    constexpr explicit vector3D(double x, double y, double z) : _crdVec{x, y, z} {}
    
    constexpr vector3D& operator=(const vector3D& oth) = default;
    
    constexpr vector3D& operator+=(const vector3D& oth) {
        for(uint8_t i = 0; i < 3; ++i)
            _crdVec[i] += oth._crdVec[i];
        return *this;
    }
    constexpr vector3D& operator-=(const vector3D& oth) {
        for(uint8_t i = 0; i < 3; ++i)
            _crdVec[i] -= oth._crdVec[i];
        return *this;
    }
    constexpr vector3D& operator*=(double coef) {
        for(uint8_t i = 0; i < 3; ++i)
            _crdVec[i] *= coef;
        return *this;
    }
    
    double& operator[](const size_t i) { return (i < 3) ? _crdVec[i] : _crdVec[0]; }
    constexpr const double& operator[](const size_t i) const { return (i < 3) ? _crdVec[i] : _crdVec[0]; }
    
    friend constexpr vector3D operator+(const vector3D& lhs, const vector3D& rhs) {
        return vector3D(lhs) += rhs;
    }
    friend constexpr vector3D operator-(const vector3D& lhs, const vector3D& rhs) {
        return vector3D(lhs) -= rhs;
    }
    friend constexpr vector3D operator*(const vector3D& lhs, double coef) {
        return vector3D(lhs) *= coef;
    }
    
    friend constexpr double operator*(const vector3D& lhs, const vector3D& rhs) {
        double res = 0.0;
        for(uint8_t i = 0; i < 3; ++i)
            res += lhs[i] * rhs[i];
        return res;
    }
    
    friend constexpr vector3D cross_product(const vector3D& lhs, const vector3D& rhs) {
        return vector3D(lhs[1]*rhs[2] - lhs[2]*rhs[1], lhs[2]*rhs[0] - lhs[0]*rhs[2], lhs[0]*rhs[1] - lhs[1]*rhs[0]);
    }
    
    friend bool operator==(const vector3D& lhs, const vector3D& rhs) {
        double res = true;
        for(uint8_t i = 0; i < 3 && res; ++i)
            res = ( std::abs(lhs[i] - rhs[i]) < 1e-12 );
        return res;
    }
    
    double getLength() const { return std::sqrt(*this * *this); }
    
    ~vector3D() = default;
};

constexpr vector3D operator "" _i(long double arg) { return vector3D(arg, 0.0, 0.0); }
constexpr vector3D operator "" _j(long double arg) { return vector3D(0.0, arg, 0.0); }
constexpr vector3D operator "" _k(long double arg) { return vector3D(0.0, 0.0, arg); }

#endif
