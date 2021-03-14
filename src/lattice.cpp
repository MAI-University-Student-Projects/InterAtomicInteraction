#include "lattice.h"

#include <numeric>

//can be parallelized with execution_policy from c++17
double Lattice::fullEnergy(const Parameters& prms, const matrix3D& dfrm) const {
    return std::transform_reduce(_atoms.cbegin(), _atoms.cend(),
                                 0.0, std::plus<>(),
                                 [&, idx = 0](const Atom& atm) mutable {
        idx++;
        return repulsiveEnergy(atm, prms, dfrm, idx - 1) + bandEnergy(atm, prms, dfrm, idx - 1);
    });
}

double Lattice::repulsiveEnergy(const Atom &atm, const Parameters &prms, const matrix3D &dfrm, size_t atm_i = 0) const {
    return std::transform_reduce(_atoms.cbegin(), _atoms.cend(),
                                 0.0, std::plus<>(),
                                 [&, atm_j = 0](const Atom& atm_oth) mutable {
        double r_ij = (atm_i == atm_j) ? 0.0 : distance(atm, atm_oth, _period, dfrm);
        InteractType intr_tp = interact_type(atm, atm_oth);
        atm_j++;
        return (atm_i == atm_j) ? 0.0 : (prms[A1_ID*intr_tp] * (r_ij - prms[R0_ID*intr_tp]) + prms[A0_ID*intr_tp]) * std::exp(-prms[P_ID*intr_tp] * (r_ij / prms[R0_ID*intr_tp] - 1));
    });
}

double Lattice::bandEnergy(const Atom &atm, const Parameters &prms, const matrix3D &dfrm, size_t atm_i = 0) const {
    return std::sqrt(std::transform_reduce(_atoms.cbegin(), _atoms.cend(),
                                           0.0, std::plus<>(),
                                           [&, atm_j = 0](const Atom& atm_oth) mutable {
        double r_ij = (atm_i == atm_j) ? 0.0 : distance(atm, atm_oth, _period, dfrm);
        InteractType intr_tp = interact_type(atm, atm_oth);
        return (atm_i == atm_j) ? 0.0 : std::pow(prms[KSI_ID*intr_tp], 2) * std::exp(-2 * prms[Q_ID*intr_tp] * (r_ij / prms[R0_ID*intr_tp] - 1));
    }) );
}

