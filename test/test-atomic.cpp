#include <catch2/catch_all.hpp>
#include "matrix3D.hpp"

using namespace Catch::literals;

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
        vec_main[0] = 5;
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
        
        REQUIRE(cross_product(vec_main, vec_oth) == 4.916_i - 3.63_j - 0.56_k);
        CHECK(vec_main * vec_oth == (5.71_a).margin(1e-12));
        CHECK(diff_const.getLength() == Catch::Approx(std::sqrt(12.4449)).margin(1e-12));
    }
}

TEST_CASE("matrix3D tests", "[math][Matrix3D]") {
    constexpr matrix3D mtx3x3{ vector3D{1, 2, 3}, vector3D{4, 5, 6}, vector3D{7, 8, 9} };
    constexpr vector3D vec3d = 3.0_i + 2.0_j + 4.0_k;
    constexpr vector3D res_prod = 19.0_i + 46.0_j + 73.0_k;
    CHECK((1.75_identity[2] == 1.75_k));
    CHECK((mtx3x3 * vec3d) == res_prod);
}
