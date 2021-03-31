#include <iostream>
#include <chrono>

#include <nlohmann/json.hpp>

#include "optimizer.h"

using namespace inter_atomic;

int main(int argc, const char * argv[]) {
    parameters in_tbl_prms = { 4.085, -2.96, 1.08, 0.97, 0.51, 0.881, 0.32, -0.18, -4.10 };
    std::pair<parameters, parameters> ptncl = {
        {
            0.06, -0.1, 0.7853, 7.2853, 2.0927, 1.9257,
            0.06, -0.1, 0.7853, 7.2853, 2.0927, 1.9257,
            0.06, -0.1, 0.7853, 7.2853, 2.0927, 1.9257
        },
        {
            0.1370, 0.1, 1.57066667, 14.5706, 4.1853, 3.8514,
            0.1370, 0.1, 1.57066667, 14.5706, 4.1853, 3.8514,
            0.1370, 0.1, 1.57066667, 14.5706, 4.1853, 3.8514
        }
    };
    std::unique_ptr<AbstrOptimazer> optimizer = std::make_unique<NelderMeadOptimizer>(in_tbl_prms, ptncl, 0.001);
    Solver solvr{std::move(optimizer)};
    auto start = std::chrono::high_resolution_clock::now();
    parameters res = solvr.solve();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::chrono::seconds::period> elapsedTime = finish - start;
    std::cout << "Elapsed time " << elapsedTime.count() << 's' << std::endl;
    std::for_each(std::begin(res), std::end(res), [](double val) { std::cout << val << ' '; });
}
