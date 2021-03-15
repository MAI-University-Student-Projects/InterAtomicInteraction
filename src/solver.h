#ifndef _SOLVER_H
#define _SOLVER_H

#include "lattice.h"

enum TableParamID : uint8_t {
    LATTICE_CONST_ID = 0,
    ECOH_ID, B_ID,
    C11_ID, C12_ID, C44_ID,
    ESOL_ID, EIN_ID, EON_ID
};

class TableParamCounter {
public:
    
private:
    std::valarray<double> _inpt_table_prms;
};

#endif
