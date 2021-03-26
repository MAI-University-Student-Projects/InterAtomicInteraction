#include "optimizer.h"

#include <map>
#include <random>
#include <omp.h>

namespace inter_atomic {

std::valarray<double> NelderMeadOptimizer::optimize(int iteration_limit) const {
    size_t n = _inpt_ptncl_bounds.first.size();
    //use openmp to parallel each optimization start with bottleneck at the convergence condition
    std::random_device rd_init_state;
    std::mt19937 gen(rd_init_state());
    for (int launch = 0; launch < iteration_limit; ++launch) {
        TableEstimator tabl_prms_cntr;
        std::multimap<double, parameters> f_x;
        for(size_t i = 0; i < n + 1; ++i) {
            parameters x_i(n);
            std::transform(std::begin(_inpt_ptncl_bounds.first), std::end(_inpt_ptncl_bounds.first), std::begin(_inpt_ptncl_bounds.second),
                           std::begin(x_i), [&gen](double x, double y) { return std::uniform_real_distribution<>(x, y)(gen); } );
            f_x.emplace(errorFunctional(x_i, tabl_prms_cntr), std::move(x_i));
        }
        
        // ...
        
    }
}

}
