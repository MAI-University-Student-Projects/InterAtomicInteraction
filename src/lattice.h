#ifndef _LATTICE_H
#define _LATTICE_H

#include "matrix3D.hpp"
#include "parameters.h"
#include <vector>
#include <omp.h>

enum InteractType { AB = 0, AA, BB };

struct Atom {
    enum AtomType { A, B };
    
    constexpr explicit Atom(AtomType tp, const vector3D& vec) : _pos(vec), _type(tp) {}
    
    template<InteractType _T>
    friend vector3D minDistance(const Atom& lhs, const Atom& rhs, const vector3D& period);
    friend bool operator==(const Atom& lhs, const Atom& rhs) { return (lhs._pos == rhs._pos) && (lhs._type == rhs._type); }
    
    vector3D _pos;
    AtomType _type;
};

template<InteractType _T>
vector3D minDistance(const Atom& lhs, const Atom& rhs, const vector3D& period) {
    if constexpr (_T != AB) {
        
    }
    
    else {
        
    }
}

class Lattice {
private:
    std::vector<Atom> _atoms;
    double _a;
    
public:
    Lattice() = default;
    explicit Lattice(std::vector<Atom> atms, double lat_const) : _atoms{std::move(atms)}, _a{lat_const} { }
    Lattice(const Lattice&) = default;
    Lattice(Lattice&&) = default;
    
    template<InteractType _T>
    double fullEnergy(const Parameters& params, const matrix3D& transform) const {
        
    }
    
    ~Lattice() = default;
};

#endif
