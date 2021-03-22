#ifndef _SOLVER_H
#define _SOLVER_H

#include <memory>

#include "lattice.h"

namespace inter_atomic {

    enum TblPrmID {
        LTTC_CONST_ID = 0,
        ECOH_ID, B_ID,
        C11_ID, C12_ID, C44_ID,
        ESOL_ID, EIN_ID, EON_ID,
        TBL_SIZE
    };

    /* только для решетки с периодом (3,3,3) как и должно быть, так как на ней и вычисляем примерные параметры => есть абстракция => необходим класс
    
     стратегия: ОценщикТабличныхПараметров ИнтерфейсМетодовОптимиацииСФункционаломОшибки МетодОптимизации Решатель
    */
    double ErrorFunctional(const std::valarray<double>& in_tbl_prms, const std::valarray<double>& ptncl_prms, const Lattice& lttc);

    class TableEstimator {
    public:
        explicit TableEstimator(double input_lttc_cnst, double coh_energy_oth) : _coh_energy_A{coh_energy_oth} {
            _lttc_ptr = std::make_unique<Lattice>(Atom::AtomType::B, std::array<int, 3>{3, 3, 3}, input_lttc_cnst);
        }
        
        std::valarray<double> estimateTblPrms(const std::valarray<double>& ptncl_prms);
    private:
        double estimateLttcConstnt(const std::valarray<double>& ptncl_prms) const;
        
        std::unique_ptr<Lattice> _lttc_ptr;
        
        const double _coh_energy_A;
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

//    double ErrorFunctional(const std::valarray<double>& ptncl_prms) {
//        return std::sqrt(
//                         std::pow( (_inpt_table_prms / CountTableParams(ptncl_prms)) - 1.0, 2)
//                         .sum() / (_inpt_table_prms.size() - 1)
//                         );
//    }
}

#endif
