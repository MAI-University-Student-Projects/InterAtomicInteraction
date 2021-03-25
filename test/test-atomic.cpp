#include <catch2/catch_all.hpp>
#include "lattice.h"

using namespace Catch::literals;
using namespace inter_atomic;

//TEST_CASE(name of test case, tags for selecting what test cases to run)
TEST_CASE("vector3D tests", "[math][vector3D]") {
    vector3D vec_main(1.0, 1.2, 1.0);
    vector3D vec_oth(0.8, 0.4, 4.43);
    constexpr vector3D summ_const = 1.8_i + 1.6_j + 5.43_k;
    constexpr vector3D diff_const = 3.43_k - 0.2_i - 0.8_j;
    SECTION("Basic operations") {
        REQUIRE(vec_main[0] == (1.0_a).margin(1e-12));
        REQUIRE(vec_main[1] == (1.2_a).margin(1e-12));
        REQUIRE(vec_main[2] == (1.0_a).margin(1e-12));
        vec_main[2] = 5;
        REQUIRE(vec_main[4] == (5.0_a).margin(1e-12));
        
        REQUIRE_FALSE((summ_const == diff_const));
    }
    //REQUIRE - stops running if failed
    //CHECK - continue running even if failed
    SECTION("Summary operations") {
        REQUIRE((vec_main + vec_oth == summ_const));
        CHECK_FALSE((vec_main == summ_const));
    }
    
    SECTION("Difference operations") {
        REQUIRE((vec_oth - vec_main == diff_const));
        CHECK_FALSE((vec_oth == diff_const));
    }
    
    SECTION("Product operations") {
        constexpr vector3D prod_res = vector3D(4.916, -3.63, -0.56);
        REQUIRE((vec_main * 2 == 2.0_i + 2.4_j + 2.0_k));
        
        REQUIRE((cross_product(vec_main, vec_oth) == 4.916_i - 3.63_j - 0.56_k));
        CHECK(vec_main * vec_oth == (5.71_a).margin(1e-12));
        CHECK(diff_const.getLength() == Catch::Approx(std::sqrt(12.4449)).margin(1e-12));
    }
}

TEST_CASE("matrix3D tests", "[math][Matrix3D]") {
    constexpr matrix3D mtx3x3{ vector3D{1.0, 2.0, 3.0}, vector3D{4.0, 5.0, 6.0}, vector3D{7.0, 8.0, 9.0} };
    constexpr vector3D vec3d = 3.0_i + 2.0_j + 4.0_k;
    constexpr vector3D res_prod = 19.0_i + 46.0_j + 73.0_k;
    CHECK((1.75_identity[2] == 1.75_k));
    constexpr matrix3D c44 = elasticity44(0.5);
    CHECK((c44[0][1] == (0.5_a).margin(1e-12) && c44[2][2] == (1.33_a).margin(1e-2)));
    CHECK((mtx3x3 * vec3d == res_prod));
}

TEST_CASE("atom tests", "[inter_atomic][Atom]") {
    constexpr Atom atom_a{Atom::AtomType::A, vector3D{0.5, 0.5, 2}};
    constexpr Atom atom_b{Atom::AtomType::B, vector3D{0, 1.0, 0.5}};
    constexpr Atom atom_period_z{Atom::AtomType::A, vector3D{0, 1.0, 3}};
    constexpr std::array period = {3, 3, 3};
    REQUIRE((interact_type(atom_period_z, atom_b) == Bond_Tp::AB && interact_type(atom_a, atom_period_z) == Bond_Tp::AA));
    double ab = distance(atom_a, atom_b, period, 1.0_identity, false);
    REQUIRE((ab == Catch::Approx(std::sqrt(2.75)).margin(1e-12) && ab == distance(atom_b, atom_a, period, 1.0_identity, false)));
    CHECK(distance(atom_period_z, atom_b, period, 1.0_identity) == Catch::Approx(std::sqrt(0.25)).margin(1e-12));
    CHECK(distance(atom_period_z, atom_a, period, elasticity44(3.0)) == (1.419_a).margin(1e-3));
}

TEST_CASE("lattice basics tests", "[inter_atomic][Lattice]") {
    constexpr std::array period = {2, 2, 2};
    double a = 2.0;
    Lattice lttc{Atom::AtomType::A, period, a};
    REQUIRE(lttc.size() == 32);
    REQUIRE(lttc.volume() == (64_a).margin(1e-12));
    REQUIRE((lttc[lttc.size() - 1]._pos == vector3D{1.5, 1.0, 1.5}));
    REQUIRE((lttc.at(1, 1, 1, 3)._pos == lttc[lttc.size() - 1]._pos));
    lttc[lttc.size() - 1]._type = Atom::AtomType::B;
    REQUIRE(interact_type(lttc[lttc.size() - 1], lttc[0]) == Bond_Tp::AB);
    lttc[lttc.size() - 1]._type = Atom::AtomType::A;
    REQUIRE(distance(lttc[lttc.size() - 1], lttc[1], period, 1.0_identity) == (1.224_a).margin(1e-3));
    lttc.set_constant(3.0);
    REQUIRE_FALSE(lttc.get_constant() == a);
}

TEST_CASE("lattice energy tests", "[inter_atomic][Lattice]") {
    Lattice lttc{Atom::AtomType::A, {1, 1, 1}, 2};
    std::valarray ptncl_prms = { 0.1, -0.1, 0.8, 8.0, 3.0, 2.0 }; // { A0_ID = 0, A1_ID, KSI_ID, P_ID, Q_ID, R0_ID, PTCL_SIZE };
    double coh_energy = lttc.cohesiveEnergy(ptncl_prms, 1.0_identity);
    REQUIRE(coh_energy == (1.618_a).margin(1e-3));
    REQUIRE(lttc.fullEnergy(ptncl_prms, 1.0_identity) == coh_energy * lttc.size());
    REQUIRE(lttc.fullEnergy(ptncl_prms, 1.0_identity, 3) == (1.7365_a).margin(1e-4));
}

