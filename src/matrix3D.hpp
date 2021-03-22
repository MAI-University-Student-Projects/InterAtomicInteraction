#ifndef _MATRIX_3D_HPP
#define _MATRIX_3D_HPP

#include "vector3D.hpp"

namespace inter_atomic {

    class matrix3D {
    private:
        vector3D _arr3x3[3] = { vector3D{}, vector3D{}, vector3D{} };
    public:
        constexpr matrix3D() = default;
        constexpr matrix3D(const matrix3D&) = default;
        
        constexpr explicit matrix3D(vector3D x, vector3D y, vector3D z) : _arr3x3{x, y, z} {}
        
        constexpr matrix3D& operator=(const matrix3D& oth) = default;
        
        vector3D& operator[](const size_t i) { return (i < 3) ? _arr3x3[i] : _arr3x3[2]; }
        constexpr const vector3D& operator[](const size_t i) const { return (i < 3) ? _arr3x3[i] : _arr3x3[2]; }
        
        friend constexpr vector3D operator*(const matrix3D& lhs, const vector3D& rhs) {
            vector3D res;
            for(uint8_t i = 0; i < 3; ++i) {
                for(uint8_t j = 0; j < 3; ++j)
                    res[i] += lhs._arr3x3[i][j] * rhs[j];
            }
            return res;
        }
        
        bool isIdentity() const {
            double res = true;
            for(uint8_t i = 0; i < 3 && res; ++i) {
                for(uint8_t j = 0; j < 3 && res; ++j)
                    res = (i == j) ? ( std::abs(_arr3x3[i][j] - 1.0) < 1e-12 ) : ( std::abs(_arr3x3[i][j] - 0.0) < 1e-12 );
            }
            return res;
        }
    };

    constexpr matrix3D operator "" _identity(long double arg) {
        return matrix3D{ 1.0_i * arg, 1.0_j * arg, 1.0_k * arg };
    }

    constexpr matrix3D bulk(double alpha) {
        return matrix3D{ 1.0_i * (alpha + 1.0), 1.0_j * (alpha + 1.0), 1.0_k * (alpha + 1.0) };
    }

    constexpr matrix3D elasticity11(double alpha) {
        return matrix3D{ 1.0_i * (alpha + 1.0), 1.0_j * (alpha + 1.0), 1.0_k };
    }

    constexpr matrix3D elasticity12(double alpha) {
        return matrix3D{ 1.0_i * (alpha + 1.0), 1.0_j * (1.0 - alpha), 1.0_k };
    }

    constexpr matrix3D elasticity44(double alpha) {
        return matrix3D{ vector3D{1.0, alpha, 0.0}, vector3D{alpha, 0.1, 0.0}, 1.0_k * (1 / (1 - alpha) / (1 - alpha)) };
    }
}

#endif
