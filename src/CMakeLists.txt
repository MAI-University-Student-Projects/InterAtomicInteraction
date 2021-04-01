find_package(OpenMP REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

set(_HEADERS optimizer.h table_estimator.h lattice.h)
set(_SOURCES optimizer.cpp table_estimator.cpp lattice.cpp atom.hpp matrix3D.hpp vector3D.hpp)

add_library(atomic_lib ${_HEADERS} ${_SOURCES})
target_link_libraries(atomic_lib PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(atomic_lib PUBLIC ${OpenMP_INCLUDE_PATH})

add_executable(inter_atomic main.cpp)

target_link_libraries(inter_atomic PRIVATE nlohmann_json::nlohmann_json atomic_lib)

set_target_properties(inter_atomic atomic_lib
PROPERTIES
RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../src
ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../src
)
