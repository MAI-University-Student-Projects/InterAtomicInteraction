set(_HEADERS lattice.h parameters.h)
set(_SOURCES vector3D.hpp matrix3D.hpp atom.hpp lattice.cpp parameters.cpp)

add_library(atomic_lib ${_HEADERS} ${_SOURCES})

add_executable(inter_atomic main.cpp)

target_link_libraries(inter_atomic atomic_lib)
