#ifndef _OPTIMIZER_H
#define _OPTIMIZER_H

#include "table_estimator.h"

// strategy through pimpl (Pointer to implementation)

namespace inter_atomic {

    class AbstrOptimazer {
    public:
        explicit AbstrOptimazer(parameters in_tbl_prms,
                                std::pair<parameters, parameters> ptncl_prms_bounds) : _inpt_table_prms{std::move(in_tbl_prms)}, _inpt_ptncl_bounds{std::move(ptncl_prms_bounds)} {}
        
        AbstrOptimazer(const AbstrOptimazer&) = delete;
        AbstrOptimazer(AbstrOptimazer&&) = default;
        
        AbstrOptimazer& operator=(const AbstrOptimazer&) = delete;
        AbstrOptimazer& operator=(AbstrOptimazer&&) = default;
        
        std::pair<double, parameters> errorFunctional(const parameters& ptncl_prms, TableEstimator& tabl_esmtr) const {
            parameters table_tmp_prms = tabl_esmtr.estimateTblPrms(ptncl_prms, _inpt_table_prms[TblPrmID::ECOH_A_ID]);
            return { std::pow( (_inpt_table_prms / table_tmp_prms) - 1.0, 2).sum() / (_inpt_table_prms.size() - 1), table_tmp_prms } ;
        }
        virtual std::pair<parameters, parameters> optimize(int iteration_limit) const = 0;
        virtual ~AbstrOptimazer() {}
    protected:
        parameters _inpt_table_prms;
        std::pair<parameters, parameters> _inpt_ptncl_bounds;
    };

    class NelderMeadOptimizer : public AbstrOptimazer {
    public:
        explicit NelderMeadOptimizer(parameters in_tbl_prms,
                                     std::pair<parameters, parameters> ptncl_prms_bounds, double epsln = 1e-4,
                                     double alpha = 1.0, double betta = 0.5, double gamma = 2, double sigma = 0.5) : AbstrOptimazer{std::move(in_tbl_prms), std::move(ptncl_prms_bounds)}, _epsln{epsln}, _alpha{alpha}, _betta{betta}, _gamma{gamma}, _sigma{sigma} {}
        std::pair<parameters, parameters> optimize(int iteration_limit) const override;
    private:
        const double _epsln, _alpha, _betta, _gamma, _sigma;
    };

    class Solver {
    public:
        Solver(std::unique_ptr<AbstrOptimazer> optimizer) : _optimizer{std::move(optimizer)} {}
        std::pair<parameters, parameters> solve(int iteration_limit = -1) {
            return _optimizer->optimize(iteration_limit);
        }
    private:
        std::unique_ptr<AbstrOptimazer> _optimizer;
    };

}

#endif
