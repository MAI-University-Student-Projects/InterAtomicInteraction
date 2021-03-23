find_package(OpenMP REQUIRED)

set(_HEADERS optimizer.h table_estimator.h lattice.h)
set(_SOURCES optimizer.cpp table_estimator.cpp lattice.cpp atom.hpp matrix3D.hpp vector3D.hpp)

add_library(atomic_lib ${_HEADERS} ${_SOURCES})
target_link_libraries(atomic_lib PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(atomic_lib PUBLIC ${OpenMP_INCLUDE_PATH})

add_executable(inter_atomic main.cpp)

target_link_libraries(inter_atomic atomic_lib)
