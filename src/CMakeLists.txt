set(_SOURCES vector3D.hpp)

# add_library(atomic_lib ${_SOURCES})

add_executable(inter_atomic main.cpp ${_SOURCES})

# target_link_libraries(inter_atomic atomic_lib)
