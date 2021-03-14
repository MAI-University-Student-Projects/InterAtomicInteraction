#ifndef _PARAMETERS_H
#define _PARAMETERS_H

enum InitParamID : uint8_t { A0_ID = 0, A1_ID, KSI_ID, P_ID, Q_ID, R0_ID};
enum TableParamID : uint8_t { LATTICE_CONST_ID = 18, ECOH_ID, B_ID, C11_ID, C12_ID, C44_ID, ESOL_ID, EIN_ID, EON_ID, TOTAL_SIZE };

class Parameters {
public:
    double& operator[] (double p);
    const double& operator[] (double p) const;
    Parameters& operator+=(const Parameters& o);
    friend Parameters operator+(const Parameters& lhs, const Parameters& rhs);
private:
    double p_[TOTAL_SIZE];
};

#endif
