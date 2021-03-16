#include "lattice.h"

#include <numeric>

//can be parallelized with execution_policy from c++17
double Lattice::fullEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm) const {
    return std::transform_reduce(
                                 _atoms.cbegin(), _atoms.cend(),
                                 0.0, std::plus<>(),
                                 [&, it = _atoms.cbegin()](const Atom&) mutable {
                                     ++it;
                                     return cohesiveEnergy(ptncl_prms, dfrm, it);
                                 } );
}

double Lattice::repulsiveEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const {
    return 2 * std::transform_reduce(
                                     atm_it, _atoms.cend(),
                                     0.0, std::plus<>(),
                                     [&](const Atom& atm_oth) {
                                         InteractType intr_tp = interact_type(*atm_it, atm_oth);
                                         double r_ij = (intr_tp == AB) ? distance(*(atm_it - 1), atm_oth, _period, dfrm, false) * _a : distance(*(atm_it - 1), atm_oth, _period, dfrm) * _a;
                                         return (ptncl_prms[A1_ID*intr_tp] * (r_ij - ptncl_prms[R0_ID*intr_tp]) + ptncl_prms[A0_ID*intr_tp]) * std::exp(-ptncl_prms[P_ID*intr_tp] * (r_ij / ptncl_prms[R0_ID*intr_tp] - 1));
                                     } );
}

double Lattice::bandEnergy(const std::valarray<double>& ptncl_prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const {
    return 2 * -std::sqrt(
                          std::transform_reduce(
                                                atm_it, _atoms.cend(),
                                                0.0, std::plus<>(),
                                                [&](const Atom& atm_oth) mutable {
                                                    InteractType intr_tp = interact_type(*atm_it, atm_oth);
                                                    double r_ij = (intr_tp == AB) ? distance(*(atm_it - 1), atm_oth, _period, dfrm, false) * _a : distance(*(atm_it - 1), atm_oth, _period, dfrm) * _a;
                                                    return std::pow(ptncl_prms[KSI_ID*intr_tp], 2) * std::exp(-2 * ptncl_prms[Q_ID*intr_tp] * (r_ij / ptncl_prms[R0_ID*intr_tp] - 1));
                                                } )
                          );
}

