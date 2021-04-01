#include <iostream>
#include <fstream>
#include <exception>

#include <nlohmann/json.hpp>
#include <omp.h>

#include "atom.h"

using json = nlohmann::json;
using namespace inter_atomic;

int main(int argc, const char * argv[]) {
    constexpr std::array<std::pair<Atom, Atom>, 3> atm_arr = {
        std::make_pair( Atom(Atom::AtomType::A, vector3D{0.0, 0.0, 0.0}),
                       Atom(Atom::AtomType::A, vector3D{0.0, 0.0, 1.0}) ) ,
        std::make_pair( Atom(Atom::AtomType::A, vector3D{0.0, 0.0, 0.0}),
                       Atom(Atom::AtomType::B, vector3D{0.0, 0.0, 1.0}) ) ,
        std::make_pair( Atom(Atom::AtomType::B, vector3D{0.0, 0.0, 0.0}),
                       Atom(Atom::AtomType::B, vector3D{0.0, 0.0, 1.0}) )
    };
    
    std::array<std::string, 3> lbl_arr = { "AA", "AB", "BB" };
    std::ifstream inp_js("optimized_params.json", std::ifstream::in);
    json jsonConfig;
    if(!inp_js.is_open()) {
        std::cerr << "File doesn't exist" << std::endl;
        exit(-1);
    }
    try {
        inp_js >> jsonConfig;
    } catch (nlohmann::detail::parse_error) {
        std::cerr << "File of unexpected type or empty, run optimization procedure first" << std::endl;
        exit(-1);
    }
    parameters ptncl_prms = jsonConfig["optimized_ptncl_params"].get<parameters>();
    
    jsonConfig.clear();
    std::ofstream outp_js("plot_data.json", std::ios::out);
    
    double a_left = 0.1;
    double a_right = 20.1;
    int limit = 1000;
    double step = (a_right - a_left) / limit;
    
    #pragma omp parallel num_threads(4)
    {
        #pragma omp single nowait
        {
            std::valarray<double> r(limit);
            std::generate(std::begin(r), std::end(r), [&, idx = 0]() mutable {
                ++idx;
                return distance(atm_arr[0].first, atm_arr[0].second) * (a_left + (idx - 1) * step);
            });
            jsonConfig["r"] = r;
        }
        
        #pragma omp for schedule(static, 3)
        for(size_t i = 0; i < atm_arr.size(); ++i) {
            std::valarray<double> u_r(limit);
            std::generate(std::begin(u_r), std::end(u_r), [&, idx = 0]() mutable {
                ++idx;
                return energy(atm_arr[i].first, atm_arr[i].second, ptncl_prms, a_left + (idx - 1) * step);
            });
            
            jsonConfig[lbl_arr[i]] = u_r;
        }
    
    }
    outp_js << jsonConfig;
    
    return 0;
}

