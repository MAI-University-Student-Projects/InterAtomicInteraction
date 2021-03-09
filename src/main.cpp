#include "matrix3D.hpp"
#include "lattice.h"

int main(int argc, char* argv[]) {
    constexpr vector3D vec_2 = 1.0_i + 5.0_j + 0.890_k;
    std::cout << vec_2[2] << std::endl;
    return 0;
}

