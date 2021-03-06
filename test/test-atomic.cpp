#include <catch2/catch_all.hpp>
#include "vector3D.hpp"

using namespace Catch::literals;

//TEST_CASE(name of test case, tags for selecting what test cases to run)
TEST_CASE("vector3D tests", "[struct][vector3D]") {
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
