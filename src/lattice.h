#ifndef _LATTICE_H
#define _LATTICE_H

#include <vector>

#include "atom.h"
#include <valarray>

enum InitParamID : uint8_t { A0_ID = 0, A1_ID, KSI_ID, P_ID, Q_ID, R0_ID };

class Lattice {
public:
    Lattice() = default;
    
    template<size_t N>
    explicit Lattice(double lat_const, const Atom (&basis)[N], const vector3D& prd) : _a{lat_const}, _period{prd} {
        _atoms.reserve(N * static_cast<size_t>(_period[0] * _period[1] * _period[2]));
        for(uint8_t i = 0; i < prd[0]; ++i) {
            for(uint8_t j = 0; j < prd[1]; ++j) {
                for(uint8_t k = 0; k < prd[2]; ++k) {
                    for(size_t b = 0; b < N; ++b)
                        _atoms.emplace_back(Atom(basis[b]._type, basis[b]._pos + vector3D(i, j, k)));
                }
            }
        }
    }
    
    Lattice(const Lattice&) = default;
    Lattice(Lattice&&) = default;
    
    double fullEnergy(const std::valarray<double>& prms, const matrix3D& dfrm) const;
    double cohesiveEnergy(const std::valarray<double>& prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const {
        return repulsiveEnergy(prms, dfrm, atm_it) + bandEnergy(prms, dfrm, atm_it);
    }
    
    const Atom& operator[](size_t i) const { return ( i < _atoms.size() ) ? _atoms[i] : _atoms[0]; };
    Atom& operator[](size_t i) { return ( i < _atoms.size() ) ? _atoms[i] : _atoms[0]; };
    
    ~Lattice() = default;
private:
    double repulsiveEnergy(const std::valarray<double>& prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const;
    double bandEnergy(const std::valarray<double>& prms, const matrix3D& dfrm, std::vector<Atom>::const_iterator atm_it) const;
    
    std::vector<Atom> _atoms;
    vector3D _period;
    double _a;
};

#endif
