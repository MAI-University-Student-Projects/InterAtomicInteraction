#ifndef _ATOM_HPP
#define _ATOM_HPP

#include "matrix3D.hpp"
#include "parameters.h"

enum InteractType : uint8_t { AA = 0, AB, BB };

struct Atom {
    enum AtomType : uint8_t { A, B };
    
    constexpr explicit Atom(AtomType tp, const vector3D& vec) : _pos{vec}, _type{tp} {}
    
    friend double distance(const Atom& lhs, const Atom& rhs, const vector3D& period, const matrix3D& m);
    
    friend constexpr InteractType interact_type(const Atom& lhs, const Atom& rhs) {
        return static_cast<InteractType>(lhs._type + rhs._type);
    }
    
    AtomType _type;
    vector3D _pos;
};

#endif
