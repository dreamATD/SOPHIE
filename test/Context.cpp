//
// Created by 69029 on 3/31/2020.
//


#include <Context.h>
#include "catch.hpp"
#include "TestUtils.h"

#define EQ(a, b) ((a) == (b))
#define NE(a, b) ((a) != (b))

class ContextTest: public Context {
public:
    ContextTest(): Context(4, 3, 55) {}
};

TEST_CASE_METHOD(ContextTest, "ksiTest", "[context]") {
    complex<double> ksiProd(1, 0);
    for (int i = 0; i < M; ++i) {
        REQUIRE(equals(ksiPow[i] * ksiInvPow[i], 1)) ;
        if (i > 0)
            REQUIRE(equals(ksiPow[i - 1] * ksiPow[1], ksiPow[i]));
        ksiProd *= ksiPow[1];
    }

    REQUIRE(equals(ksiProd, 1));
}


TEST_CASE_METHOD(ContextTest, "revTest", "[context]") {
    for (int n = 0; n < Nh; ++n) {
        sizeType x = n, y = 0;
        for (int j = 0; j < logNh; ++j, x >>= 1) y = y << 1 ^ (x & 1);
        REQUIRE(revh[n] == y);
    }
}

TEST_CASE_METHOD(ContextTest, "primeTest", "[context]") {
    for (int i = 0; i < K; ++i) REQUIRE(mod(pVec[i], M) == 1);
    for (int j = 0; j < L; ++j) REQUIRE(mod(qVec[j], M) == 1);

    for (int i = 0; i < K; ++i)
        REQUIRE(isPrime(pVec[i]));

    for (int j = 0; j < L; ++j)
        REQUIRE(isPrime(qVec[j]));
}

TEST_CASE_METHOD(ContextTest, "pAndqTest", "[context]") {
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < L; ++j) {
            REQUIRE(mulMod(pModq[i][j], pInvModq[i][j], qVec[j]) == 1);
            REQUIRE(mulMod(qModp[j][i], qInvModp[j][i], pVec[i]) == 1);
        }

    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j)
            if (i != j)
                REQUIRE(mulMod(qVec[i], qInvModq[i][j], qVec[j]) == 1);
}

TEST_CASE_METHOD(ContextTest, "prodTest", "[context]") {
    for (int j = 0; j < L; ++j)
        REQUIRE(mulMod(PModq[j], PInvModq[j], qVec[j]) == 1);
}

TEST_CASE_METHOD(ContextTest, "hatTest", "[context]") {
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < L; ++j) {
            REQUIRE(mulMod(pHatModq[i][j], PInvModq[j], qVec[j]) == pInvModq[i][j]);
            for (int k = 0; k < L; ++k)
                if (j <= k)
                    REQUIRE(mulMod(qHatModp[k][j][i], qModp[j][i], pVec[i]) == QModp[k][i]);
        }

    for (int i = 0; i < K; ++i) {
        valueType cur = 1;
        for (int j = 0; j < K; ++j)
            if (i != j)
                mulModAssign(cur, cur, power(pVec[j], pVec[i] - 2, pVec[i]), pVec[i]);
        REQUIRE(cur == pHatInvModp[i]);
    }
}

TEST_CASE_METHOD(ContextTest, "rootAndPowTest", "[context]") {
    for (int i = 0; i < L; ++i)
        for (int j = 0; j <= i; ++j) {
            valueType cur = 1;
            for (int k = 0; k <= i; ++k)
                if (k != j)
                    mulModAssign(cur, cur, power(qVec[k], qVec[j] - 2, qVec[j]), qVec[j]);
            REQUIRE(EQ(cur, qHatInvModq[i][j]));

        }

    for (int i = 0; i < K; ++i) {
        REQUIRE(EQ(power(pRoots[i], M, pVec[i]), 1));
        for (int m = 1; m < M; ++m)
            REQUIRE(NE(power(pRoots[i], m, pVec[i]), 1));
    }

    for (int j = 0; j < L; ++j) {
        REQUIRE(EQ(power(qRoots[j], M, qVec[j]), 1));
        for (int m = 1; m < M; ++m)
            REQUIRE(NE(power(qRoots[j], m, qVec[j]), 1));
    }

    for (int i = 0; i < K; ++i) {
        REQUIRE(EQ(pPsi[i][0], 1));
        REQUIRE(EQ(pInvPsi[i][0], 1));
        for (int j = 1; j < M; ++j) {
            REQUIRE(EQ(pPsi[i][j], mulMod(pPsi[i][j - 1], pRoots[i], pVec[i])));
            REQUIRE(EQ(mulMod(pPsi[i][j], pInvPsi[i][j], pVec[i]), 1));
        }
    }

    for (int j = 0; j < L; ++j) {
        REQUIRE(EQ(qPsi[j][0], 1));
        REQUIRE(EQ(qInvPsi[j][0], 1));
        for (int i = 1; i < M; ++i) {
            REQUIRE(EQ(qPsi[j][i], mulMod(qPsi[j][i - 1], qRoots[j], qVec[j])));
            REQUIRE(EQ(mulMod(qPsi[j][i], qInvPsi[j][i], qVec[j]), 1));
        }
    }

    for (int i = 0; i < Nh; ++i) {
        for (int j = 0; j < N; ++j) {
            auto tmp1 = addMod(qPsi[0][j * rotGroup[i] % M], qInvPsi[0][j * rotGroup[i] % M], qVec[0]);
            auto tmp2 = addMod(qPsi[0][(N - j) * rotGroup[i] % M], qInvPsi[0][(N - j) * rotGroup[i] % M], qVec[0]);
            REQUIRE(addMod(tmp1, tmp2, qVec[0]) == 0);
        }
    }
}