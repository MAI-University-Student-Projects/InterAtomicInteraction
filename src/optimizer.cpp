#include "optimizer.h"

#include <tuple>
#include <random>
#include <fstream>

#include <omp.h>

namespace inter_atomic {

    std::pair<parameters, parameters> NelderMeadOptimizer::optimize(int iteration_limit) const {
        int stop = (iteration_limit > 0) ? iteration_limit : std::numeric_limits<int16_t>::max();
        const size_t n = _inpt_ptncl_bounds.first.size();
        volatile bool to_finish = false;
        std::pair<parameters, parameters> result;
        std::ofstream logger("nelder_mead_log.txt", std::ios::out);
        //use openmp to parallel each optimization start with bottleneck at the convergence condition
        #pragma omp parallel
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            TableEstimator tabl_prms_cntr;
            // vector <error, ptncl_prms, table_prms_for_ptncl_prms>
            std::vector< std::tuple<double, parameters, parameters> > f_x;
            f_x.reserve(n + 1);
            for(size_t i = 0; i < n + 1; ++i) {
                parameters x_i(n);
                std::transform(std::begin(_inpt_ptncl_bounds.first), std::end(_inpt_ptncl_bounds.first), std::begin(_inpt_ptncl_bounds.second),
                               std::begin(x_i), [&gen](double x, double y) { return std::uniform_real_distribution<>(x, y)(gen); } );
                auto err_res = errorFunctional(x_i, tabl_prms_cntr);
                f_x.emplace_back(err_res.first, std::move(x_i), std::move(err_res.second));
            }
            auto it_f_h = f_x.begin();
            auto it_f_g = std::next(it_f_h);
            auto it_f_l = std::prev(f_x.end());
            double converg = std::numeric_limits<double>::max();
            
            #pragma omp for schedule(dynamic, 1)
            for (int launch = 0; launch < stop; ++launch) {
                if(to_finish)
                    continue;
                std::sort(f_x.begin(), f_x.end(), [](auto &lhs_f_x, auto &rhs_f_x) { return std::get<0>(lhs_f_x) > std::get<0>(rhs_f_x); });
                bool to_shrink = false;

                // works in godbolt clang 11.0.0, doesn't work on my clang 11.0.1, may be because of openmp
//                parameters x_c = std::accumulate(it_f_g, f_x.end(), parameters(f_x[0].second.size()),
//                                  [&f_x](const parameters& res, auto &pp_x_i) {
//                                                    return res + pp_x_i.second / (f_x.size() - 1);
//                                                });
                parameters x_c(n);
                for(auto it = it_f_g; it != f_x.end(); ++it)
                    x_c += std::get<1>(*it) / (f_x.size() - 1);

                parameters x_r = (1 + _alpha) * x_c - _alpha * std::get<1>(*it_f_h);
                auto refl_res = errorFunctional(x_r, tabl_prms_cntr);
                if(refl_res.first < std::get<0>(*it_f_l)) {
                    parameters x_e = (1 - _gamma) * x_c + _gamma * x_r;
                    auto expand_res = errorFunctional(x_r, tabl_prms_cntr);
                    *it_f_h = (expand_res.first < refl_res.first) ? std::make_tuple(expand_res.first, std::move(x_e), std::move(expand_res.second)) : std::make_tuple(expand_res.first, std::move(x_r), std::move(refl_res.second));
                }
                else if(refl_res.first > std::get<0>(*it_f_l) and refl_res.first < std::get<0>(*it_f_h)) {
                    if(refl_res.first > std::get<0>(*it_f_g))
                        to_shrink = true;
                    *it_f_h = std::make_tuple(refl_res.first, std::move(x_r), std::move(refl_res.second));
                }
                else
                    to_shrink = true;

                if(to_shrink) {
                    parameters x_s = _betta * std::get<1>(*it_f_h) + (1 - _betta) * x_c;
                    auto shrink_res = errorFunctional(x_s, tabl_prms_cntr);
                    if(shrink_res.first < std::get<0>(*it_f_h))
                        *it_f_h = std::make_tuple(shrink_res.first, std::move(x_s), std::move(shrink_res.second));
                    else
                        std::for_each(f_x.begin(), it_f_l, [&](auto &f_x_i) {
                            parameters x_i = std::get<1>(*it_f_l) + _sigma * (std::get<1>(f_x_i) - std::get<1>(*it_f_l));
                            auto res = errorFunctional(x_i, tabl_prms_cntr);
                            f_x_i = std::make_tuple(res.first, std::move(x_i), std::move(res.second));
                        });
                }
                double mean = std::accumulate(f_x.cbegin(), f_x.cend(), 0.0, [&f_x](double res, auto &f) { return res + std::get<0>(f) / f_x.size(); });
                converg = std::sqrt(
                                    std::transform_reduce(f_x.cbegin(), f_x.cend(), f_x.cbegin(), 0.0,
                                                          std::plus<>(),
                                                          [mean](auto &lhs_f, auto &rhs_f) { return (std::get<0>(lhs_f) - mean) * (std::get<0>(rhs_f) - mean);}
                                                          )
                                    );
                #pragma omp critical(stopCriteria)
                logger << "i = " << launch << ": func_min=" << std::get<0>(*it_f_l) << ", func_max=" << std::get<0>(*it_f_h) << ", cnvg = " << converg << ", thr=" << omp_get_thread_num() << "\n";
                if((converg < _epsln) || (launch == iteration_limit - 1) || std::get<0>(*it_f_l) < _epsln) {
                    to_finish = true;
                    result = std::make_pair(std::get<1>(*it_f_l), std::get<2>(*it_f_l));
                }
            }
        }
        return result;
    }

}
