#ifndef _SOLVER_H
#define _SOLVER_H

#include <memory>

#include "lattice.h"

enum TableParamID : uint8_t {
    LATTICE_CONST_ID = 0,
    ECOH_ID, B_ID,
    C11_ID, C12_ID, C44_ID,
    ESOL_ID, EIN_ID, EON_ID
};

class TableParamCounter {
public:
    TableParamCounter(Lattice* latt, std::valarray<double> _table_prms) : _lttc_ptr{latt}, _inpt_table_prms{std::move(_table_prms)} { }
    double ErrorFunctional(const std::valarray<double>& ptncl_prms) {
        return std::sqrt(
                         std::pow( (_inpt_table_prms / CountTableParams(ptncl_prms)) - 1.0, 2)
                         .sum() / (_inpt_table_prms.size() - 1)
                         );
    }
private:
    std::valarray<double> CountTableParams(const std::valarray<double>& ptncl_prms);
    
    std::valarray<double> _inpt_table_prms;
    const Lattice* _lttc_ptr;
};

#endif
