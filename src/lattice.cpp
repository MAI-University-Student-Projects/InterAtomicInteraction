#include "lattice.h"

#include <numeric>

//can be parallelized with execution_policy from c++17
double Lattice::fullEnergy(const std::valarray<double>& prms, const matrix3D& dfrm) const {
    return std::transform_reduce(
                                 _atoms.cbegin(), _atoms.cend(),
                                 0.0, std::plus<>(),
                                 [&, it = _atoms.cbegin()]() mutable {
                                     ++it;
                                     return cohesiveEnergy(prms, dfrm, it);
                                 } );
}

double Lattice::repulsiveEnergy(const std::valarray<double>& prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const {
    return 2 * std::transform_reduce(
                                     atm_it, _atoms.cend(),
                                     0.0, std::plus<>(),
                                     [&](const Atom& atm_oth) {
                                         double r_ij = distance(*(atm_it - 1), atm_oth, _period, dfrm) * _a;
                                         InteractType intr_tp = interact_type(*atm_it, atm_oth);
                                         return (prms[A1_ID*intr_tp] * (r_ij - prms[R0_ID*intr_tp]) + prms[A0_ID*intr_tp]) * std::exp(-prms[P_ID*intr_tp] * (r_ij / prms[R0_ID*intr_tp] - 1));
                                     } );
}

double Lattice::bandEnergy(const std::valarray<double>& prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const {
    return 2 * -std::sqrt(
                          std::transform_reduce(
                                                atm_it, _atoms.cend(),
                                                0.0, std::plus<>(),
                                                [&](const Atom& atm_oth) mutable {
                                                    double r_ij = distance(*(atm_it - 1), atm_oth, _period, dfrm) * _a;
                                                    InteractType intr_tp = interact_type(*atm_it, atm_oth);
                                                    return std::pow(prms[KSI_ID*intr_tp], 2) * std::exp(-2 * prms[Q_ID*intr_tp] * (r_ij / prms[R0_ID*intr_tp] - 1));
                                                } )
                         );
}

