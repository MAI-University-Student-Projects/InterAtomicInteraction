#ifndef _TABLE_ESTIMATOR_H
#define _TABLE_ESTIMATOR_H

#include <memory>

#include "lattice.h"

namespace inter_atomic {

    enum TblPrmID {
        LTTC_CONST_ID = 0,
        ECOH_ID, B_ID,
        C11_ID, C12_ID, C44_ID,
        ESOL_ID, EIN_ID, EON_ID,
        ECOH_A_ID,
        TBL_SIZE
    };

    class TableEstimator {
    public:
        explicit TableEstimator(double input_lttc_cnst = 1.0) {
            _lttc_ptr = std::make_unique<Lattice>(Atom::AtomType::B, std::array<int, 3>{3, 3, 3}, input_lttc_cnst);
        }
        
        parameters estimateTblPrms(const parameters& ptncl_prms, double coh_energy_oth);
    private:
        double estimateLttcConstnt(const parameters& ptncl_prms, double a_left, double a_right, double epsln = 0.01) const;
        
        std::unique_ptr<Lattice> _lttc_ptr;
        
        static constexpr double _delta = 0.001;
        static constexpr matrix3D _ident_dfrm = 1.0_identity;

        static constexpr matrix3D _b_dfrm_pos = bulk(_delta);
        static constexpr matrix3D _b_dfrm_neg = bulk(-_delta);

        static constexpr matrix3D _c11_dfrm_pos = elasticity11(_delta);
        static constexpr matrix3D _c11_dfrm_neg = elasticity11(-_delta);

        static constexpr matrix3D _c12_dfrm_pos = elasticity12(_delta);
        static constexpr matrix3D _c12_dfrm_neg = elasticity12(-_delta);

        static constexpr matrix3D _c44_dfrm_pos = elasticity44(_delta);
        static constexpr matrix3D _c44_dfrm_neg = elasticity44(-_delta);
    };

}

#endif
