#include "solver.h"

namespace inter_atomic {

    //    possible parallel, but many lattice setters - better not here
    std::valarray<double> TableEstimator::estimateTblPrms(const std::valarray<double>& ptncl_prms) {
        std::valarray<double> table_res(TblPrmID::TBL_SIZE);
        //count lattice constant here estimateLttcConstnt(ptncl_prms); set_constant =
        //table_res[TblPrmID::LTTC_CONST_ID] =
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
        table_res[TblPrmID::ESOL_ID] = mixed_energy - lttc_energy - _coh_energy_A + coh_energy;
        _lttc_ptr->operator[](0)._type = Atom::AtomType::B; //change atom type back
        
        // lttc NON-CONST here
        _lttc_ptr->at(0, 0, 1, 3)._type = Atom::AtomType::A; // change type of atom in surface layer
        double in_adatom_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2));
        _lttc_ptr->at(0, 0, 1, 2)._type = Atom::AtomType::A; // making dimer of A atoms in surface layer
        double in_dim_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2));
        table_res[TblPrmID::EIN_ID] = in_dim_surf_energy - surface_energy - 2 * (in_adatom_surf_energy - surface_energy);
        _lttc_ptr->at(0, 0, 1, 3)._type = Atom::AtomType::B; // change atom type back
        _lttc_ptr->at(0, 0, 1, 2)._type = Atom::AtomType::B; // change atom type back
        
        // lttc NON-CONST here
        _lttc_ptr->at(0, 0, 2)._type = Atom::AtomType::A; // change type of atom on surface layer
        double on_adatom_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2, 1));
        _lttc_ptr->at(0, 0, 2, 1)._type = Atom::AtomType::A; // making dimer of A atoms on surface layer
        double on_dim_surf_energy = _lttc_ptr->fullEnergy(ptncl_prms, _ident_dfrm, _lttc_ptr->id_at(0, 0, 2, 2));
        table_res[TblPrmID::EON_ID] = on_dim_surf_energy - surface_energy - 2 * (on_adatom_surf_energy - surface_energy);
        _lttc_ptr->at(0, 0, 2)._type = Atom::AtomType::B; // change atom type back
        _lttc_ptr->at(0, 0, 2, 1)._type = Atom::AtomType::B; // change atom type back
    }

    double TableEstimator::estimateLttcConstnt(const std::valarray<double>& ptncl_prms) const {
        // use gradient descent
    }

//double ErrorFunctional(const std::valarray<double>& in_tbl_prms, const std::valarray<double>& ptncl_prms, Lattice& lttc) {
//    return std::sqrt(
//                     std::pow( (in_tbl_prms / CountTblPrms(ptncl_prms, lttc)) - 1.0, 2)
//                     .sum() / in_tbl_prms.size()
//                     );
//}

}
