#include "atom.h"

namespace inter_atomic {

    double energy(const Atom& lhs, const Atom& rhs, const parameters& ptncl_prms, double lttc_cnst) {
        Bond_Tp intr_tp = interact_type(lhs, rhs);
        double r_ij = distance(lhs, rhs) * lttc_cnst;
        auto tp_prm_it = std::next(std::begin(ptncl_prms), intr_tp * PtclPrmID::PTCL_SIZE);
        return (tp_prm_it[PtclPrmID::A1_ID] * (r_ij - tp_prm_it[PtclPrmID::R0_ID]) + tp_prm_it[PtclPrmID::A0_ID]) * std::exp(-tp_prm_it[PtclPrmID::P_ID] * (r_ij / tp_prm_it[PtclPrmID::R0_ID] - 1)) - std::sqrt(std::pow(tp_prm_it[PtclPrmID::KSI_ID], 2) * std::exp(-2 * tp_prm_it[PtclPrmID::Q_ID] * (r_ij / tp_prm_it[PtclPrmID::R0_ID] - 1)));
    }

}
