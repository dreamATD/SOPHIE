//
// Created by 69029 on 4/1/2020.
//


#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <CommonHeader.h>
#include <Utils.h>

const valueType B = 10000000000000000;

TEST_CASE("powerTest", "[arithmetic operation]") {
    valueType x, y, z;
    getRandomInRange(x, (valueType) 1, B);
    getRandomInRange(y, (valueType) 1, B);
    getRandomInRange(z, (valueType) 1, B);
    __int128 nx = x, ny = y, nz = z, ans = 1;
    __int128 tmp = nx;
    ans = 1;
    for (__int128 i = ny; i; i >>= 1) {
        if (i & 1) ans = ans * tmp % nz;
        tmp = tmp * tmp % nz;
    }
    uint64_t ans2 = (uint64_t) ans;
    REQUIRE(power(x, y, z) == ans2);

}

TEST_CASE("addTest", "[arithmetic operation]") {
    valueType x, y, z;
    getRandomInRange(x, (valueType) 1, B);
    getRandomInRange(y, (valueType) 1, B);
    getRandomInRange(z, (valueType) 1, B);
    __int128 nx = x, ny = y, nz = z, ans;
    valueType res;
    ans = (nx + ny) % nz;
    addModAssign(res, x, y, z);
    uint64_t ans2 = (uint64_t) ans;
    REQUIRE(addMod(x, y, z) == ans2);
    REQUIRE(res == ans2);
}

TEST_CASE("subTest", "[arithmetic operation]") {
    valueType x, y, z;
    getRandomInRange(x, (valueType) 1, B);
    getRandomInRange(y, (valueType) 1, B);
    getRandomInRange(z, (valueType) 1, B);
    __int128 nx = x, ny = y, nz = z, ans = 1;
    valueType res;
    ans = (nx % nz - ny % nz + nz) % nz;
    subModAssign(res, x, y, z);
    uint64_t ans2 = (uint64_t) ans;
    REQUIRE(subMod(x, y, z) == ans2);
    REQUIRE(res == ans2);
}

TEST_CASE("mulTest", "[arithmetic operation]") {
    valueType x, y, z;
    getRandomInRange(x, (valueType) 1, B);
    getRandomInRange(y, (valueType) 1, B);
    getRandomInRange(z, (valueType) 1, B);
    __int128 nx = x, ny = y, nz = z, ans = 1;
    valueType res;
    ans = (nx * ny) % nz;
    mulModAssign(res, x, y, z);
    uint64_t ans2 = (uint64_t) ans;
    REQUIRE(mulMod(x, y, z) == ans2);
    REQUIRE(res == ans2);
}




