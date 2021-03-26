#ifndef _LATTICE_H
#define _LATTICE_H

#include <vector>
#include <valarray>

#include "atom.hpp"

namespace inter_atomic {

    using parameters = std::valarray<double>;

    enum PtclPrmID { A0_ID = 0, A1_ID, KSI_ID, P_ID, Q_ID, R0_ID, PTCL_SIZE };

    class Lattice {
    public:
        Lattice() = default;
        
        explicit Lattice(Atom::AtomType tp, const std::array<int, 3>& prd, double lat_const) : _period{prd}, _a{lat_const} {
            _atoms.reserve(_unit_vtcs.size() * (_period[0]+1) * (_period[1]+1) * (_period[2]+1));
            for(int k = 0; k < prd[2]; ++k) {
                for(int j = 0; j < prd[1]; ++j) {
                    for(int i = 0; i < prd[0]; ++i) {
                        for(size_t b = 0; b < _unit_vtcs.size(); ++b)
                            _atoms.emplace_back(Atom{tp, _unit_vtcs[b] + vector3D(i, j, k)});
                    }
                }
            }
        }
        
        Lattice(const Lattice&) = delete;
        Lattice(Lattice&&) = default;
        
        Lattice& operator=(const Lattice&) = delete;
        Lattice& operator=(Lattice&&) = default;

    //    fullEnergy(parameters, deformation, end())
        double fullEnergy(const parameters& ptncl_prms, const matrix3D& dfrm, size_t atm_id_r = 0) const;
    //    cohesiveEnergy(parameters, deformation, begin(), end())
        double cohesiveEnergy(const parameters& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l = 0, size_t atm_id_r = 0) const;
        
        // atom period id-cies i = ; j = ; k = ; unit_cell_idx =
        // index taken from the way lattice constructed - watch nested loops in explicit constructor
        const Atom& at(size_t i, size_t j, size_t k, size_t unit_id = 0) const { return _atoms[id_at(i, j, k, unit_id)]; }
        Atom& at(size_t i, size_t j, size_t k, size_t unit_id = 0) { return _atoms[id_at(i, j, k, unit_id)]; }
        size_t id_at(size_t i, size_t j, size_t k, size_t unit_id = 0) const {
            return k*_period[1]*_period[0]*_unit_vtcs.size() + j*_period[0]*_unit_vtcs.size() + i*_unit_vtcs.size() + unit_id;
        }
        
        const Atom& operator[](size_t i) const { return ( i < _atoms.size() ) ? _atoms[i] : _atoms[_atoms.size() - 1]; };
        Atom& operator[](size_t i) { return ( i < _atoms.size() ) ? _atoms[i] : _atoms[_atoms.size() - 1]; };
        
        double volume() const { return std::pow(_a, 3) * _period[0] * _period[1] * _period[2]; }
        size_t size() const { return _atoms.size(); }
        
        void set_constant(double val) { _a = val; }
        double get_constant() const { return _a; }
        
        ~Lattice() = default;
    private:
        
        double repulsiveEnergy(const parameters& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l, size_t atm_id_r) const;
        double bandEnergy(const parameters& ptncl_prms, const matrix3D& dfrm, size_t atm_id_l, size_t atm_id_r) const;
        
        std::vector<Atom> _atoms;
        std::array<int, 3> _period;
        double _a;
        
        static constexpr std::array<vector3D, 4> _unit_vtcs = {
            vector3D{0.0, 0.0, 0.0}, vector3D{0.5, 0.5, 0.0},
            vector3D{0.0, 0.5, 0.5}, vector3D{0.5, 0.0, 0.5}
        };
    };

}

#endif
