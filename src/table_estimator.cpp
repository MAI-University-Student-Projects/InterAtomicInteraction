#include "table_estimator.h"

namespace inter_atomic {

    //    possible parallel, but many lattice setters - better not here
    parameters TableEstimator::estimateTblPrms(const parameters& ptncl_prms, double coh_energy_oth) {
        parameters table_res(TblPrmID::TBL_SIZE);
        table_res[TblPrmID::ECOH_A_ID] = coh_energy_oth;
        
        // lttc NON-CONST here
        double a = estimateLttcConstnt(ptncl_prms, 1, 8);
        _lttc_ptr->set_constant(a);
        table_res[TblPrmID::LTTC_CONST_ID] = a;
        
        double lttc_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm);
        double coh_energy = lttc_energy / _lttc_ptr->size();
        table_res[TblPrmID::ECOH_ID] = coh_energy;
        double V0 = _lttc_ptr->volume();
        // indecies taken from the way lattice constructed - watch nested loops in lattice constructor
        // able to do that, because table parameters are estimated on a certain lattice, with certain period (3, 3, 3) and atom type B
        // id_at(period_x_id, period_y_id, period_z_id, current_unit_cell_id)
        double surface_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2));
        
        //when openmp, use default shared to avoid copies of variables, global var-s implicitly shared
        //lambda should be private in parallel block because of ref grab
        auto second_deriv = [&](const matrix3D& neg_dfrm, const matrix3D& pos_dfrm) -> double {
            return (_lttc_ptr->fullEnergy(ptncl_prms, neg_dfrm) - 2 * lttc_energy + _lttc_ptr->fullEnergy(ptncl_prms, pos_dfrm)) / (_delta * _delta);
        };
        table_res[TblPrmID::B_ID] = 1.602 * second_deriv(_b_dfrm_neg, _b_dfrm_pos) / (9.0 * V0);
        double deriv_C11 = second_deriv(_c11_dfrm_neg, _c11_dfrm_pos);
        double deriv_C12 = second_deriv(_c12_dfrm_neg, _c12_dfrm_pos);
        table_res[TblPrmID::C11_ID] = 1.602 * (deriv_C11 + deriv_C12) / (4.0 * V0);
        table_res[TblPrmID::C12_ID] = 1.602 * (deriv_C11 - deriv_C12) / (4.0 * V0);
        table_res[TblPrmID::C44_ID] = 1.602 * second_deriv(_c44_dfrm_neg, _c44_dfrm_pos) / (4.0 * V0);
        
        // lttc NON-CONST here
        _lttc_ptr->operator[](0)._type = Atom::AtomType::A; //change type of any atom
        double mixed_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm);
        table_res[TblPrmID::ESOL_ID] = mixed_energy - lttc_energy - coh_energy_oth + coh_energy;
        _lttc_ptr->operator[](0)._type = Atom::AtomType::B; //change atom type back
        
        // lttc NON-CONST here
        _lttc_ptr->at(0, 0, 1, 3)._type = Atom::AtomType::A; // change type of atom in surface layer
        double in_adatom_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2));
        _lttc_ptr->at(1, 1, 1, 2)._type = Atom::AtomType::A; // making dimer of A atoms in surface layer
        double in_dim_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2));
        table_res[TblPrmID::EIN_ID] = in_dim_surf_energy - surface_energy - 2 * (in_adatom_surf_energy - surface_energy);
        _lttc_ptr->at(0, 0, 1, 3)._type = Atom::AtomType::B; // change atom type back
        _lttc_ptr->at(1, 1, 1, 2)._type = Atom::AtomType::B; // change atom type back
        
        // lttc NON-CONST here
        _lttc_ptr->at(0, 0, 2)._type = Atom::AtomType::A; // change type of atom on surface layer
        double on_adatom_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2, 1));
        _lttc_ptr->at(0, 0, 2, 1)._type = Atom::AtomType::A; // making dimer of A atoms on surface layer
        double on_dim_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2, 2));
        table_res[TblPrmID::EON_ID] = on_dim_surf_energy - surface_energy - 2 * (on_adatom_surf_energy - surface_energy);
        _lttc_ptr->at(0, 0, 2)._type = Atom::AtomType::B; // change atom type back
        _lttc_ptr->at(0, 0, 2, 1)._type = Atom::AtomType::B; // change atom type back
        return table_res;
    }

    double TableEstimator::estimateLttcConstnt(const parameters& ptncl_prms, double a_left, double a_right, double epsln) const {
        const uint8_t knots_size = 5;
        double cnt_epsln = std::numeric_limits<double>::max();
        std::pair cnt_min_pair = std::make_pair(0.0, std::numeric_limits<double>::max());
        std::pair prev_min_pair = std::make_pair(0.0, std::numeric_limits<double>::max());
        while(cnt_epsln > epsln) {
            double step = (a_right - a_left) / knots_size;
            for(uint8_t i = 0; i < knots_size; ++i) {
                double a_knot = a_left + i * step;
                _lttc_ptr->set_constant(a_knot);
                double knot_energy = _lttc_ptr->cohesiveEnergy(ptncl_prms, _ident_dfrm);
                if(knot_energy < cnt_min_pair.second) {
                    cnt_min_pair.first = a_knot;
                    cnt_min_pair.second = knot_energy;
                }
            }
            a_left = cnt_min_pair.first - step / 2;
            a_right = cnt_min_pair.first - step / 2;
            cnt_epsln = std::abs(prev_min_pair.second - cnt_min_pair.second);
            prev_min_pair = cnt_min_pair;
        }
        return cnt_min_pair.first;
    }

}
