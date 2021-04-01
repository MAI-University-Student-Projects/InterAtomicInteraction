#ifndef _ATOM_H
#define _ATOM_H

#include <array>
#include <valarray>

#include "matrix3D.hpp"

namespace inter_atomic {

    using parameters = std::valarray<double>;

    enum PtclPrmID { A0_ID = 0, A1_ID, KSI_ID, P_ID, Q_ID, R0_ID, PTCL_SIZE };

    enum Bond_Tp { AA = 0, AB, BB };

    struct Atom {
        enum AtomType { A, B };
        
        constexpr explicit Atom(AtomType tp, const vector3D& vec) : _pos{vec}, _type{tp} {}
        
        friend double distance(const Atom& lhs, const Atom& rhs,
                               std::array<int, 3> period = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()},
                               const matrix3D& m = 1.0_identity, bool z_periodic = true) {
            vector3D diff = lhs._pos - rhs._pos;
            const uint8_t end_i = z_periodic ? 3 : 2;
            for(uint8_t i = 0; i < end_i; ++i) {
                if(diff[i] > 0.5 * period[i])
                    diff[i] -= period[i];
                else if (diff[i] < -0.5 * period[i])
                    diff[i] += period[i];
            }
            return m.isIdentity() ? diff.getLength() : (m * diff).getLength();
        }
        
        friend double energy(const Atom& lhs, const Atom& rhs, const parameters& ptncl_prms, double lttc_cnst);
        
        friend constexpr Bond_Tp interact_type(const Atom& lhs, const Atom& rhs) {
            return static_cast<Bond_Tp>(lhs._type + rhs._type);
        }
        
        AtomType _type;
        vector3D _pos;
    };

}

#endif
