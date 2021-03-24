#ifndef _OPTIMIZER_H
#define _OPTIMIZER_H

#include "table_estimator.h"

// strategy through pimpl (Pointer to implementation)

namespace inter_atomic {

    class AbstrOptimazer {
    public:
        explicit AbstrOptimazer(std::valarray<double> in_tbl_prms,
                                std::pair<std::valarray<double>, std::valarray<double>> ptncl_prms_bounds) : _inpt_table_prms{std::move(in_tbl_prms)}, _inpt_ptncl_bounds{std::move(ptncl_prms_bounds)} {}
        
        AbstrOptimazer(const AbstrOptimazer&) = delete;
        AbstrOptimazer(AbstrOptimazer&&) = default;
        
        AbstrOptimazer& operator=(const AbstrOptimazer&) = delete;
        AbstrOptimazer& operator=(AbstrOptimazer&&) = default;
        
        virtual double errorFunctional(const std::valarray<double>& ptncl_prms, TableEstimator& tabl_esmtr) {
        return std::sqrt(
                         std::pow( (_inpt_table_prms / tabl_esmtr.estimateTblPrms(ptncl_prms)) - 1.0, 2)
                         .sum() / (_inpt_table_prms.size() - 1)
                         );
        }
        virtual std::valarray<double> optimize(TableEstimator& tabl_esmtr) const = 0;
        virtual ~AbstrOptimazer() {}
    protected:
        std::valarray<double> _inpt_table_prms;
        std::pair<std::valarray<double>, std::valarray<double>> _inpt_ptncl_bounds;
    };

    class NelderMeadOptimizer : public AbstrOptimazer {
    public:
        explicit NelderMeadOptimizer(std::valarray<double> in_tbl_prms,
                                     std::pair<std::valarray<double>, std::valarray<double>> ptncl_prms_bounds,
                                     double alpha = 1.0, double betta = 0.5, double gamma = 2, double sigma = 0.5) : AbstrOptimazer(in_tbl_prms, ptncl_prms_bounds), _alpha{alpha}, _betta{betta}, _gamma{gamma}, _sigma{sigma} {}
        std::valarray<double> optimize(TableEstimator& tabl_esmtr) const override;
    private:
        const double _alpha, _betta, _gamma, _sigma;
    };

    class Solver {
    public:
        Solver(std::unique_ptr<AbstrOptimazer> optimizer) : _optimizer{std::move(optimizer)} {}
        std::valarray<double> solve(TableEstimator& tabl_esmtr) {
            return _optimizer->optimize(tabl_esmtr);
        }
    private:
        std::unique_ptr<AbstrOptimazer> _optimizer;
    };

}

#endif
