#include "atom.h"

double distance(const Atom& lhs, const Atom& rhs, const vector3D& period, const matrix3D& m) {
    vector3D diff = lhs._pos - rhs._pos;
    for(uint8_t i = 0; i < 3; ++i) {
        if(diff[i] > 0.5 * period[i])
            diff[i] -= period[i];
        else if (diff[i] < -0.5 * period[i])
            diff[i] += period[i];
    }
    
    return (m * diff).getLength();
}
