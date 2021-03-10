#include "matrix3D.hpp"
#include "lattice.h"
#include <iostream>

int main(int argc, const char * argv[]) {
#pragma omp parallel
#pragma omp critical
    std::cout << "Greetings from thread "<< omp_get_thread_num() << std::endl;
    return 0;
}

