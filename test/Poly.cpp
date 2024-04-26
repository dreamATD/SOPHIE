//
// Created by 69029 on 3/31/2020.
//

#include "TestUtils.h"
#include <Context.h>
#include "catch.hpp"
#include "Poly.h"

class PolyTest;

PolyTest evaluate(sizeType logNh, const vector<vector<valueType>> &q, sizeType L, const vector<valueType> &qVec,
                  const vector<vector<valueType>> &qPsi, const vector<sizeType> &rev, const vector<sizeType> &rotGroup,
                  const vector<vector<valueType>> &p = {}, sizeType K = 0, const vector<valueType> &pVec = {},
                  const vector<vector<valueType>> &pPsi = {});


/**
 * The polynomial class evaluates in the brute force manner.
 * */

class PolyTest: public Poly {
public:
    PolyTest(sizeType logn, vector<vector<valueType>> qdata, vector<vector<valueType>> pdata) :
            Poly(logn, move(qdata), move(pdata)) {
    }

    PolyTest(sizeType logn, sizeType l, const vector<double> &data, const vector<valueType> &qVec) :
            Poly(logn, l, data, qVec) {}

    PolyTest(sizeType logn, sizeType l, const vector<complex<double>> &data, const vector<valueType> &qVec) :
            Poly(logn, l, data, qVec) {}

    PolyTest(const PolyTest &other) = default;

    explicit PolyTest(Poly &&other) : Poly(forward<Poly>(other)) {}

    PolyTest(sizeType logn, sizeType l, sizeType k) : Poly(logn, l, k) {}

    using Poly::operator ==;

    friend ostream &operator << (ostream &out, const PolyTest &p) {
        out << '\n';
        printMat(out, "pData", p.pData);
        printMat(out, "qData", p.qData);
        return out;
    }

    friend PolyTest simpleAdd(const PolyTest &a, const PolyTest &b, const vector<valueType> &qVec,
                              const vector<valueType> &pVec) {
        PolyTest c(a);
        for (int j = 0; j < c.L; ++j)
            for (int n = 0; n < c.N; ++n)
                addModAssign(c.qData[j][n], c.qData[j][n], b.qData[j][n], qVec[j]);

        for (int i = 0; i < c.K; ++i)
            for (int n = 0; n < c.N; ++n)
                addModAssign(c.pData[i][n], c.pData[i][n], b.pData[i][n], pVec[i]);
        return c;
    }

    friend PolyTest simpleSub(const PolyTest &a, const PolyTest &b, const vector<valueType> &qVec,
                              const vector<valueType> &pVec) {
        PolyTest c(a);
        for (int j = 0; j < c.L; ++j)
            for (int n = 0; n < c.N; ++n)
                subModAssign(c.qData[j][n], c.qData[j][n], b.qData[j][n], qVec[j]);

        for (int i = 0; i < c.K; ++i)
            for (int n = 0; n < c.N; ++n)
                subModAssign(c.pData[i][n], c.pData[i][n], b.pData[i][n], pVec[i]);
        return c;
    }

    friend PolyTest simpleNeg(const PolyTest &a, const vector<valueType> &qVec, const vector<valueType> &pVec) {
        PolyTest c(a);
        for (int j = 0; j < c.L; ++j)
            for (int n = 0; n < c.N; ++n)
                negModAssign(c.qData[j][n], a.qData[j][n], qVec[j]);

        for (int i = 0; i < c.K; ++i)
            for (int n = 0; n < c.N; ++n)
                negModAssign(c.pData[i][n], a.pData[i][n], pVec[i]);

        return c;
    }

    friend PolyTest simpleMul(PolyTest a, PolyTest b,const vector<sizeType> &rotGroup,
                              const vector<sizeType> &rotGroupRev, const vector<sizeType> &rev,
                              const vector<valueType> &qVec, const vector<vector<valueType>> &qInvPsi,
                              const vector<valueType> &pVec, const vector<vector<valueType>> &pInvPsi) {
        sizeType logNh = a.logN, Nh = a.N, K = min(a.K, b.K), L = min(a.L, b.L);
        vector<vector<valueType>> pm(K, vector<valueType>(2 * Nh - 1, 0));
        vector<vector<valueType>> qm(L, vector<valueType>(2 * Nh - 1, 0));
        auto p1 = a.pData;
        auto p2 = b.pData;
        auto q1 = a.qData;
        auto q2 = b.qData;
        for (int i = 0; i < K; ++i) {
            valueType inv2 = pVec[i] + 1 >> 1;
            mulModAssign(p1[i][0], p1[i][0], inv2, pVec[i]);
            mulModAssign(p2[i][0], p2[i][0], inv2, pVec[i]);
        }

        for (int j = 0; j < L; ++j) {
            valueType inv2 = qVec[j] + 1 >> 1;
            mulModAssign(q1[j][0], q1[j][0], inv2, qVec[j]);
            mulModAssign(q2[j][0], q2[j][0], inv2, qVec[j]);
        }

        for (int i = 0; i < K; ++i) {
            valueType modular = pVec[i];
            for (int n1 = 0; n1 < Nh; ++n1)
                for (int n2 = 0; n2 < Nh; ++n2) {
                    valueType tmp = mulMod(p1[i][n1], p2[i][n2], modular);
                    int id1 = n1 + n2;
                    int id2 = n1 - n2;
                    if (id2 < 0) id2 = -id2;
                    addModAssign(pm[i][id1], pm[i][id1], tmp, modular);
                    addModAssign(pm[i][id2], pm[i][id2], tmp, modular);
                }
        }

        for (int j = 0; j < L; ++j) {
            valueType modular = qVec[j];
            for (int n1 = 0; n1 < Nh; ++n1)
                for (int n2 = 0; n2 < Nh; ++n2) {
                    valueType tmp = mulMod(q1[j][n1], q2[j][n2], modular);
                    int id1 = n1 + n2;
                    int id2 = n1 - n2;
                    if (id2 < 0) id2 = -id2;
                    addModAssign(qm[j][id1], qm[j][id1], tmp, modular);
                    addModAssign(qm[j][id2], qm[j][id2], tmp, modular);
                }
        }

        PolyTest res = evaluate(logNh, qm, L, qVec, qInvPsi, rev, rotGroup, pm, K, pVec, pInvPsi);
        res.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
        return res;
    }
};

PolyTest evaluate(sizeType logNh, const vector<vector<valueType>> &q, sizeType L, const vector<valueType> &qVec,
                  const vector<vector<valueType>> &qInvPsi, const vector<sizeType> &rev,
                  const vector<sizeType> &rotGroup, const vector<vector<valueType>> &p, sizeType K,
                  const vector<valueType> &pVec, const vector<vector<valueType>> &pInvPsi) {
    sizeType N = 1 << logNh + 1, Nh = 1 << logNh, M = N << 1;
    vector<vector<valueType>> qData(L);
    vector<vector<valueType>> pData(K);

    for (int j = 0; j < L; ++j) {
        qData[j].resize(Nh);
        for (int n = 0; n < Nh; ++n) {
            int idx = rotGroup[n];
            valueType &cur = qData[j][n], modular = qVec[j];
            cur = 0;
            for (int m = 0; m < N - 1; ++m) {
                int id = idx * m % M;
                valueType tmp = mulMod(q[j][m], addMod(qInvPsi[j][M - id], qInvPsi[j][id], modular), modular);
                addModAssign(cur, cur, tmp, modular);
            }
        }
        arrayBitReverse(qData[j], rev);
    }

    for (int i = 0; i < K; ++i) {
        pData[i].resize(Nh);
        for (int n = 0; n < Nh; ++n) {
            int idx = rotGroup[n];
            valueType &cur = pData[i][n], modular = pVec[i];
            cur = 0;
            for (int m = 0; m < N - 1; ++m) {
                int id = idx * m % M;
                valueType tmp = mulMod(p[i][m], addMod(pInvPsi[i][M - id], pInvPsi[i][id], modular), modular);
                addModAssign(cur, cur, tmp, modular);
            }
        }

        arrayBitReverse(pData[i], rev);
    }

    return PolyTest(logNh, qData, pData);
}

/**
 * The assistant class used to get the protected members and methods of Context.
 * */

class ContextForPolyTest: public Context {
public:
    ContextForPolyTest(): Context(11, 2, 55) {}
};

TEST_CASE_METHOD(ContextForPolyTest, "ntt&intt", "[poly]") {
    PolyTest p(logNh, L, K);
    p.sampleUniform(qVec, pVec, COEFF_VALUE);
    PolyTest np(p);

    np.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    np.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);

    REQUIRE(p == np);
}

TEST_CASE_METHOD(ContextForPolyTest, "non-assigned arithmetic operations", "[arithmetic operations]") {
    PolyTest a(logNh, L, K);
    PolyTest b(logNh, L, K);
    a.sampleUniform(qVec, pVec, COEFF_VALUE);
    b.sampleUniform(qVec, pVec, COEFF_VALUE);

    PolyTest ans4a = simpleAdd(a, b, qVec, pVec);
    PolyTest ans4s = simpleSub(a, b, qVec, pVec);
    PolyTest ans4n = simpleNeg(a, qVec, pVec);
    PolyTest ans4m = simpleMul(a, b, rotGroup, rotGroupRev, revh, qVec, qInvPsi, pVec, pInvPsi);

    a.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    b.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    auto pa = a.add(b, qVec, pVec);
    auto ps = a.sub(b, qVec, pVec);
    auto pn = a.neg(qVec, pVec);
    auto pm = a.mul(b, qVec, pVec);

    pa.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
    ps.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
    pn.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
    pm.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
//    ans4a.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
//    ans4s.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
//    ans4n.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
//    ans4m.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);

    REQUIRE(pa == ans4a);
    REQUIRE(ps == ans4s);
    REQUIRE(pn == ans4n);
    REQUIRE(pm == ans4m);
}

TEST_CASE_METHOD(ContextForPolyTest, "assigned arithmetic operations", "[arithmetic operations]") {
    PolyTest a(logNh, L, K);
    PolyTest b(logNh, L, K);
    a.sampleUniform(qVec, pVec, COEFF_VALUE);
    b.sampleUniform(qVec, pVec, COEFF_VALUE);

    PolyTest ans4a = simpleAdd(a, b, qVec, pVec);
    PolyTest ans4s = simpleSub(a, b, qVec, pVec);
    PolyTest ans4m = simpleMul(a, b, rotGroup, rotGroupRev, revh, qVec, qInvPsi, pVec, pInvPsi);

    const vector<sizeType> &rev = revh;
    Poly aa(a), as(a), am(a);

    b.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);

    aa.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    as.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    am.DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
    aa.addAndEqual(b, qVec, pVec);
    as.subAndEqual(b, qVec, pVec);
    am.mulAndEqual(b, qVec, pVec);
    aa.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
    as.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);
    am.iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);

    REQUIRE(aa == ans4a);
    REQUIRE(as == ans4s);
    REQUIRE(am == ans4m);
}
