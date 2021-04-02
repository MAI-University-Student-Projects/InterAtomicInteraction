#include <iostream>
#include <chrono>
#include <fstream>
#include <exception>
#include <cstring>

#include <nlohmann/json.hpp>

#include "optimizer.h"

using json = nlohmann::json;
using namespace inter_atomic;

int main(int argc, const char * argv[]) {
    const uint16_t iter_limit = (argc > 1) ? std::stoi(argv[1]) : 1000;
    
    std::ifstream inp_js("../../../src/input_data.json", std::ifstream::in);
    json jsonConfig;
    if(!inp_js.is_open()) {
        std::cerr << "Invalid filename" << std::endl;
        exit(-1);
    }
    inp_js >> jsonConfig;
    
    parameters in_tbl = jsonConfig["init_table_params"].get<parameters>();
    parameters ptncl_left(PtclPrmID::PTCL_SIZE * 3);
    parameters ptncl_right(PtclPrmID::PTCL_SIZE * 3);
    for(int i = 0; i < 3; ++i) {
        ptncl_left[std::slice(i * PtclPrmID::PTCL_SIZE, PtclPrmID::PTCL_SIZE, 1)] = jsonConfig["ptncl_params_left_bound"].get<parameters>();
        ptncl_right[std::slice(i * PtclPrmID::PTCL_SIZE, PtclPrmID::PTCL_SIZE, 1)] = jsonConfig["ptncl_params_right_bound"].get<parameters>();
    }
    jsonConfig.clear();
    
    std::unique_ptr<AbstrOptimazer> optimizer = std::make_unique<NelderMeadOptimizer>(std::move(in_tbl), std::make_pair(std::move(ptncl_left), std::move(ptncl_right)), 0.001);
    Solver solvr{std::move(optimizer)};
    
    auto start = std::chrono::high_resolution_clock::now();
    parameters res = solvr.solve(iter_limit);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::chrono::seconds::period> elapsedTime = finish - start;
    
    jsonConfig["optimized_ptncl_params"] = res;
    jsonConfig["elapsed time"] = std::to_string(elapsedTime.count()) + "sec";
    jsonConfig["nOMPThreads"] = std::getenv("OMP_NUM_THREADS");
    std::ofstream outp_js("../../../src/optimized_params.json", std::ios::out);
    if(!outp_js.is_open()) {
        std::cerr << "Invalid filename" << std::endl;
        std::cerr << "Unsaved result of optimization:" << std::endl;
        std::for_each(std::begin(res), std::end(res), [](double vl) { std::cerr << vl << ' '; });
        std::cerr << std::endl;
        exit(-1);
    }
    outp_js << jsonConfig.dump(4);
    
    return 0;
}
