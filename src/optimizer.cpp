#include "optimizer.h"

#include <map>
#include <random>
#include <limits>
#include <fstream>

#include <omp.h>

namespace inter_atomic {

    parameters NelderMeadOptimizer::optimize(int iteration_limit) const {
        int stop = (iteration_limit > 0) ? iteration_limit : std::numeric_limits<int16_t>::max();
        const size_t n = _inpt_ptncl_bounds.first.size();
        volatile bool to_finish = false;
        parameters result(n);
        std::ofstream logger("nelder_mead_log.txt", std::ios::out);
        //use openmp to parallel each optimization start with bottleneck at the convergence condition
        #pragma omp parallel
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            TableEstimator tabl_prms_cntr;
            std::vector< std::pair<double, parameters> > f_x;
            f_x.reserve(n + 1);
            for(size_t i = 0; i < n + 1; ++i) {
                parameters x_i(n);
                std::transform(std::begin(_inpt_ptncl_bounds.first), std::end(_inpt_ptncl_bounds.first), std::begin(_inpt_ptncl_bounds.second),
                               std::begin(x_i), [&gen](double x, double y) { return std::uniform_real_distribution<>(x, y)(gen); } );
                f_x.emplace_back(errorFunctional(x_i, tabl_prms_cntr), std::move(x_i));
            }
            auto it_f_h = f_x.begin();
            auto it_f_g = std::next(it_f_h);
            auto it_f_l = std::prev(f_x.end());
            double converg = std::numeric_limits<double>::max();
            
            #pragma omp for schedule(dynamic, 1)
            for (int launch = 0; launch < stop; ++launch) {
                if(to_finish)
                    continue;
                std::sort(f_x.begin(), f_x.end(), [](auto &lhs_f_x, auto &rhs_f_x) { return lhs_f_x.first > rhs_f_x.first; });
                bool to_shrink = false;

                // works in godbolt clang 11.0.0, doesn't work on my clang 11.0.1, may be because of openmp
//                parameters x_c = std::accumulate(it_f_g, f_x.end(), parameters(f_x[0].second.size()),
//                                  [&f_x](const parameters& res, auto &pp_x_i) {
//                                                    return res + pp_x_i.second / (f_x.size() - 1);
//                                                });
                parameters x_c(n);
                for(auto it = it_f_g; it != f_x.end(); ++it)
                    x_c += (*it).second / (f_x.size() - 1);

                parameters x_r = (1 + _alpha) * x_c - _alpha * (*it_f_h).second;
                double f_r = errorFunctional(x_r, tabl_prms_cntr);
                if(f_r < (*it_f_l).first) {
                    parameters x_e = (1 - _gamma) * x_c + _gamma * x_r;
                    double f_e = errorFunctional(x_e, tabl_prms_cntr);
                    *it_f_h = (f_e < f_r) ? std::make_pair(f_e, std::move(x_e)) : std::make_pair(f_e, std::move(x_r));
                }
                else if(f_r > (*it_f_l).first and f_r < (*it_f_h).first) {
                    if(f_r > (*it_f_g).first)
                        to_shrink = true;
                    *it_f_h = std::make_pair(f_r, std::move(x_r));
                }
                else
                    to_shrink = true;

                if(to_shrink) {
                    parameters x_s = _betta * (*it_f_h).second + (1 - _betta) * x_c;
                    double f_s = errorFunctional(x_s, tabl_prms_cntr);
                    if(f_s < (*it_f_h).first)
                        *it_f_h = std::make_pair(f_s, std::move(x_s));
                    else
                        std::for_each(f_x.begin(), it_f_l, [&](auto &f_x_i) {
                            parameters x_i = (*it_f_l).second + _sigma * (f_x_i.second - (*it_f_l).second);
                            f_x_i = std::make_pair(errorFunctional(x_i, tabl_prms_cntr), std::move(x_i));
                        });
                }
                double mean = std::accumulate(f_x.cbegin(), f_x.cend(), 0.0, [&f_x](double res, auto &f) { return res + f.first / f_x.size(); });
                converg = std::sqrt(
                                    std::transform_reduce(f_x.cbegin(), f_x.cend(), f_x.cbegin(), 0.0,
                                                          std::plus<>(),
                                                          [mean](auto &lhs_f, auto &rhs_f) { return (lhs_f.first - mean) * (rhs_f.first - mean);}
                                                          )
                                    );
                #pragma omp critical(stopCriteria)
                logger << "i = " << launch << ": func_min=" << (*it_f_l).first << ", func_max=" << (*it_f_h).first << ", cnvg = " << converg << ", thr=" << omp_get_thread_num() << "\n";
                if((converg < _epsln) || (launch == iteration_limit - 1)) {
                    to_finish = true;
                    result = (*it_f_l).second;
                }
            }
        }
        return result;
    }

}
