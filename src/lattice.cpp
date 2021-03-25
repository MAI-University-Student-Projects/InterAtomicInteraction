#include "lattice.h"

#include <numeric>

namespace inter_atomic {

    //can be parallelized with execution_policy from c++17
    double Lattice::fullEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, size_t atm_id_r) const {
        size_t end_idx = (atm_id_r == 0) ? _atoms.size() : atm_id_r;
        return std::transform_reduce(
                                     _atoms.cbegin(), _atoms.cbegin() + end_idx,
                                     0.0, std::plus<>(),
                                     [&, idx = 0](const Atom&) mutable {
                                         ++idx;
                                         return cohesiveEnergy(ptncl_prms, dfrm, idx - 1, end_idx);
                                     } );
    }

    double Lattice::cohesiveEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l, size_t atm_id_r) const {
        size_t end_idx = (atm_id_r == 0) ? _atoms.size() : atm_id_r;
        return repulsiveEnergy(ptncl_prms, dfrm, atm_id_l, end_idx) + bandEnergy(ptncl_prms, dfrm, atm_id_l, end_idx);
    }

    double Lattice::repulsiveEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l, size_t atm_id_r) const {
        return std::transform_reduce(
                                     _atoms.cbegin(), _atoms.cbegin() + atm_id_r,
                                     0.0, std::plus<>(),
                                     [&, idx = 0](const Atom& atm_oth) mutable {
                                         double res = 0.0;
                                         if(idx != atm_id_l) {
                                             Bond_Tp intr_tp = interact_type(_atoms[atm_id_l], atm_oth);
                                             double r_ij = (intr_tp == AB) ? distance(_atoms[atm_id_l], atm_oth, _period, dfrm, false) * _a : distance(_atoms[atm_id_l], atm_oth, _period, dfrm) * _a;
                                             auto tp_prm_it = std::begin(ptncl_prms) + intr_tp * PtclPrmID::PTCL_SIZE;
                                             res = (tp_prm_it[PtclPrmID::A1_ID] * (r_ij - tp_prm_it[PtclPrmID::R0_ID]) + tp_prm_it[PtclPrmID::A0_ID]) * std::exp(-tp_prm_it[PtclPrmID::P_ID] * (r_ij / tp_prm_it[PtclPrmID::R0_ID] - 1));
                                         }
                                         ++idx;
                                         return res;
                                     } );
    }

    double Lattice::bandEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l, size_t atm_id_r) const {
        return -std::sqrt(
                          std::transform_reduce(
                                                _atoms.cbegin(), _atoms.cbegin() + atm_id_r,
                                                0.0, std::plus<>(),
                                                [&, idx = 0](const Atom& atm_oth) mutable {
                                                double res = 0.0;
                                                if(idx != atm_id_l) {
                                                    Bond_Tp intr_tp = interact_type(_atoms[atm_id_l], atm_oth);
                                                    double r_ij = (intr_tp == Bond_Tp::AB) ? distance(_atoms[atm_id_l], atm_oth, _period, dfrm, false) * _a : distance(_atoms[atm_id_l], atm_oth, _period, dfrm) * _a;
                                                    auto tp_prm_it = std::begin(ptncl_prms) + intr_tp * PtclPrmID::PTCL_SIZE;
                                                    res = std::pow(tp_prm_it[PtclPrmID::KSI_ID], 2) * std::exp(-2 * tp_prm_it[PtclPrmID::Q_ID] * (r_ij / tp_prm_it[PtclPrmID::R0_ID] - 1));
                                                }
                                                ++idx;
                                                return res;
                                                } )
                          );
    }

}

