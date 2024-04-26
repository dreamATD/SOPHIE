//
// Created by 69029 on 3/22/2020.
//

#include <algorithm>
#include <iostream>
#include <utility>
#include "Poly.h"
#include "Utils.h"

// TODO: improve DFT & iDFT

Poly::Poly(sizeType logn, sizeType l, sizeType k, bool sta):
        logN(logn), N(1L << logn), L(l), K(k),
        qData(L, vector<valueType>(N, 0)),
        pData(K, vector<valueType>(N, 0)),
        status(sta) {
}

Poly::Poly(sizeType logn, vector<vector<valueType>> qdata, bool sta):
        logN(logn), N(1L << logn), qData(move(qdata)), K(0), status(sta) {
    L = qData.size();
}

Poly::Poly(sizeType logn, vector<vector<valueType>> qdata, vector<vector<valueType>> pdata, bool sta):
        logN(logn), N(1L << logn), qData(move(qdata)), pData(move(pdata)), status(sta) {
    L = qData.size();
    K = pData.size();
}

// TODO: improve the mod operation.

Poly::Poly(sizeType logn, sizeType l, const vector<double> &data, const vector<valueType> &qVec, bool sta):
        logN(logn), N(1 << logn), L(l), K(0), qData(L, vector<valueType>(N)), status(sta) {

    for (int j = 0; j < L; ++j) {
        valueType q = qVec[j];
        for (int n = 0; n < N; ++n) {
            double re = data[n];
            qData[j][n] = re >= 0 ? (valueType) re % q : q - (valueType) fabs(re) % q ;
        }
    }

    status = COEFF_VALUE;
}

Poly::Poly(sizeType logn, sizeType l, const vector<complex<double>> &data, const vector<valueType> &qVec, bool sta):
        logN(logn), N(1 << logn), L(l), qData(L, vector<valueType>(N)), K(0), status(sta) {

    for (int j = 0; j < L; ++j) {
        valueType q = qVec[j];
        for (int n = 0; n < (N >> 1); ++n) {
            double re = data[n].real();
            double im = data[n].imag();
            qData[j][n] = re >= 0 ? (valueType) re : (valueType) (re + qVec[j]);
            qData[j][n + (N >> 1)] = im >= 0 ? (valueType) im : (valueType) (im + qVec[j]);
        }
    }
}

Poly &Poly::operator = (Poly &&other) noexcept {
    if (this != &other) {
        assert(N == other.N);
        L = other.L;
        K = other.K;
        status = other.status;
        qData = move(other.qData);
        pData = move(other.pData);
    }
    return *this;
}

bool Poly::operator == (const Poly &other) const {
    if (N != other.N || L != other.L || K != other.K || status != other.status) return false;
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            if (qData[j][n] != other.qData[j][n]) return false;

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            if (pData[i][n] != other.pData[i][n]) return false;
    return true;
}

Poly Poly::ZERO(sizeType logn, sizeType l, sizeType k) {
    return Poly(logn, l, k);
}

void Poly::setL(sizeType l) {
    L = l;
}

sizeType Poly::getL() const {
    return L;
}

vector<vector<valueType>> Poly::getPData() const {
    return pData;
}

void Poly::setPData(const vector<vector<valueType>> &pData) {
    Poly::pData = pData;
}

vector<vector<valueType>> Poly::getQData() const {
    return qData;
}

void Poly::setQData(const vector<vector<valueType>> &qData) {
    Poly::qData = qData;
}

const sizeType Poly::getN() const {
    return N;
}

// need to be zero poly.
void Poly::sampleHWT(sizeType h, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    while (h--) {
        sizeType idx;
        getRandomInRange(idx, (sizeType) 0, N - 1);

        if (!qData[0][idx]) {
            int hwt = randomTwo(-1, 1);
            for (int j = 0; j < L; ++j)
                qData[j][idx] = hwt > 0 ? hwt : qVec[j] + hwt;

            for (int i = 0; i < K; ++i)
                pData[i][idx] = hwt > 0 ? hwt : pVec[i] + hwt;
        }
    }
    status = COEFF_VALUE;
}

void Poly::sampleUniform(const vector<valueType> &qVec, const vector<valueType> &pVec, bool sta) {
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            getRandomInRange(qData[j][n], (valueType) 0, qVec[j] - 1);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            getRandomInRange(pData[i][n], (valueType) 0, pVec[i] - 1);
    status = sta;
}

void Poly::sampleGauss(double sigma, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    unsigned int bigNum = 0xffffffff;
    for (int n = 0; n < N; n += 2) {
        unsigned int tmp1, tmp2;
        getRandomInRange(tmp1, 1u, bigNum);
        getRandomInRange(tmp2, 1u, bigNum);
        double r1 = 1.0 * tmp1 / (1.0 + bigNum);
        double r2 = 1.0 * tmp2 / (1.0 + bigNum);
        double theta = 2 * M_PI * r1;
        double rr = sqrt(-2 * log(r2)) * sigma;

        long g1 = floor(rr * cos(theta) + 0.5);
        long g2 = floor(rr * sin(theta) + 0.5);

        for (int j = 0; j < L; ++j) {
            qData[j][n] = g1 >= 0 ? g1 : qVec[j] + g1;
            qData[j][n + 1] = g2 >= 0 ? g2 : qVec[j] + g2;
        }

        for (int i = 0; i < K; ++i) {
            pData[i][n] = g1 >= 0 ? g1 : pVec[i] + g1;
            pData[i][n + 1] = g2 >= 0 ? g2 : pVec[i] + g2;
        }
    }
    status = COEFF_VALUE;
}

void Poly::sampleZO(const vector<valueType> &qVec, const vector<valueType> &pVec) {
    for (int n = 0; n < N; ++n) {
        int zo = randomTwo(-1, 1);

        for (int j = 0; j < L; ++j)
            qData[j][n] = zo >= 0 ? zo : qVec[j] + zo;

        for (int i = 0; i < K; ++i)
            pData[i][n] = zo >= 0 ? zo : pVec[i] + zo;
    }
    status = COEFF_VALUE;
}

void Poly::DFT(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
               const vector<valueType> &pVec, const vector<vector<valueType>> &pPsi) {
    assert(status == COEFF_VALUE);
    for (int j = 0; j < L; ++j)
        ntt(qData[j], qVec[j], rotGroupRev, qPsi[j]);

    for (int i = 0; i < K; ++i)
        ntt(pData[i], pVec[i], rotGroupRev, pPsi[i]);
    status = POINT_VALUE;
}

void Poly::iDFT(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec,
                const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec,
                const vector<vector<valueType>> &pInvPsi) {
    assert(status == POINT_VALUE);
    for (int j = 0; j < L; ++j)
        intt(qData[j], qVec[j], rotGroupRev, qInvPsi[j]);

    for (int i = 0; i < K; ++i)
        intt(pData[i], pVec[i], rotGroupRev, pInvPsi[i]);
    status = COEFF_VALUE;
}

Poly Poly::add(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) const {
    assert(N == other.N && K == other.K);
    Poly res(logN, min(L, other.L), K);

    for (int j = 0; j < res.L; ++j)
        for (int n = 0; n < N; ++n)
            addModAssign(res.qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < res.K; ++i)
        for (int n = 0; n < res.N; ++n)
            addModAssign(res.pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);

    return res;
}

void Poly::addAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    assert(N == other.N && K == other.K);
    L = min(L, other.L);

    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            addModAssign(qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            addModAssign(pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);
}

Poly Poly::sub(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) const {
    assert(N == other.N && K == other.K);
    Poly res(logN, min(L, other.L), K);

    for (int j = 0; j < res.L; ++j)
        for (int n = 0; n < N; ++n)
            subModAssign(res.qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < res.K; ++i)
        for (int n = 0; n < res.N; ++n)
            subModAssign(res.pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);

    return res;
}

void Poly::subAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    assert(N == other.N && K == other.K);
    L = min(L, other.L);

    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            subModAssign(qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            subModAssign(pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);
}

void Poly::subAndEqual2(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    assert(N == other.N && K == other.K);
    L = min(L, other.L);

    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            subModAssign(qData[j][n], other.qData[j][n], qData[j][n], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            subModAssign(pData[i][n], other.pData[i][n], pData[i][n], pVec[i]);
}

Poly Poly::neg(const vector<valueType> &qVec, const vector<valueType> &pVec) const {
    Poly res(logN, L, K);
    for (int j = 0; j < L; ++j) {
        valueType modular = qVec[j];
        for (int n = 0; n < N; ++n)
            res.qData[j][n] = modular - qData[j][n];
    }

    for (int i = 0; i < K; ++i) {
        valueType modular = pVec[i];
        for (int n = 0; n < N; ++n)
            res.pData[i][n] = modular - pData[i][n];
    }

    return res;
}

void Poly::negAndEqual(const vector<valueType> &qVec, const vector<valueType> &pVec) {
    for (int j = 0; j < L; ++j) {
        valueType modular = qVec[j];
        for (int n = 0; n < N; ++n)
            qData[j][n] = modular - qData[j][n];
    }

    for (int i = 0; i < K; ++i) {
        valueType modular = pVec[i];
        for (int n = 0; n < N; ++n)
            pData[i][n] = modular - pData[i][n];
    }
}

Poly Poly::mul(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) const {
    assert(N == other.N && K == other.K && status == POINT_VALUE && other.status == POINT_VALUE);
    Poly res(logN, min(L, other.L), K);

    for (int j = 0; j < res.L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(res.qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < res.K; ++i)
        for (int n = 0; n < res.N; ++n)
            mulModAssign(res.pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);

    return res;
}

void Poly::mulAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    assert(N == other.N && K == other.K && status == POINT_VALUE && other.status == POINT_VALUE);
    L = min(L, other.L);

    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(qData[j][n], qData[j][n], other.qData[j][n], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            mulModAssign(pData[i][n], pData[i][n], other.pData[i][n], pVec[i]);
}

Poly Poly::mulAndProductToq(const Poly &other, const vector<valueType> &qVec) const {
    assert(N == other.N && status == POINT_VALUE && other.status == POINT_VALUE);
    Poly res(logN, min(L, other.L));

    for (int j = 0; j < res.L; ++j) {
        valueType modular = qVec[j];
        for (int n = 0; n < N; ++n)
            mulModAssign(res.qData[j][n], qData[j][n], other.qData[j][n], modular);
    }

    return res;
}

void Poly::mulAndEqualProductToq(const Poly &other, const vector<valueType> &qVec) {
    assert(N == other.N && status == POINT_VALUE && other.status == POINT_VALUE);

    L = min(L, other.L);
    K = 0;
    for (int j = 0; j < L; ++j) {
        valueType modular = qVec[j];
        for (int n = 0; n < N; ++n)
            mulModAssign(qData[j][n], qData[j][n], other.qData[j][n], modular);
    }
    pData.clear();


}

Poly Poly::imul(const vector<valueType> &qOther, const vector<valueType> &qVec, const vector<valueType> &pOther,
                const vector<valueType> &pVec) {
    Poly res(logN, L, K);
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(res.qData[j][n], qData[j][n], qOther[j], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            mulModAssign(res.pData[i][n], pData[i][n], pOther[i], pVec[i]);

    return res;
}

void Poly::imulAndEqual(const vector<valueType> &qOther, const vector<valueType> &qVec, const vector<valueType> &pOther,
                        const vector<valueType> &pVec) {
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(qData[j][n], qData[j][n], qOther[j], qVec[j]);

    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            mulModAssign(pData[i][n], pData[i][n], pOther[i], pVec[i]);
}

vector<double> Poly::convertToDouble(valueType scale, valueType modular) {
    vector<double> res(N);
    const vector<valueType> &mdata = qData[0];
    valueType mh = modular >> 1u;
    for (int n = 0; n < N; ++n) {
        res[n] = (mdata[n] <= mh ? (double) mdata[n] : -(double) (modular - mdata[n]));
        res[n] /= scale;
    }
    return res;
}

void Poly::zeroExpand(sizeType k) {
    K = k;
    pData = vector<vector<valueType>>(K, vector<valueType>(N, 0));
}

void Poly::modUp(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
                 const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec,
                 const vector<vector<valueType>> &pPsi, sizeType k, const vector<valueType> &qHatInvModq,
                 const vector<vector<valueType>> &qHatModp) {
    iDFT(rotGroupRev, qVec, qInvPsi);

    valueType tmp[L][N];
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(tmp[j][n], qData[j][n], qHatInvModq[j], qVec[j]);

    if (!K) pData = vector<vector<valueType>>(k, vector<valueType>(N));
    K = k;
    for (int i = 0; i < K; ++i) {
        for (int n = 0; n < N; ++n) {
            valueType cur = 0;

            for (int j = 0; j < L; ++j)
                mulAddModAssign(cur, tmp[j][n], qHatModp[j][i], cur, pVec[i]);

            pData[i][n] = cur;
        }
    }

    DFT(rotGroupRev, qVec, qPsi, pVec, pPsi);
}

void Poly::modDown(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
                   const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec,
                   const vector<vector<valueType>> &pInvPsi, const vector<valueType> &pHatInvModp,
                   const vector<vector<valueType>> &pHatModq, const vector<valueType> &PInvModq) {
    iDFT(rotGroupRev, qVec, qInvPsi, pVec, pInvPsi);

    valueType tmp[K][N];
    for (int i = 0; i < K; ++i)
        for (int n = 0; n < N; ++n)
            mulModAssign(tmp[i][n], pData[i][n], pHatInvModp[i], pVec[i]);

    for (int j = 0; j < L; ++j) {
        for (int n = 0; n < N; ++n) {
            valueType cur = 0;

            for (int i = 0; i < K; ++i)
                mulAddModAssign(cur, tmp[i][n], pHatModq[i][j], cur, qVec[j]);

            mulModAssign(qData[j][n], PInvModq[j], subMod(qData[j][n], cur, qVec[j]), qVec[j]);
        }
    }

    K = 0;
    pData.clear();

    DFT(rotGroupRev, qVec, qPsi);
}

void Poly::rescale(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
                   const vector<vector<valueType>> &qInvPsi, const vector<valueType> &qInvModq, bool doDFT) {
    iDFT(rotGroupRev, qVec, qInvPsi);

    --L;
    for (int j = 0; j < L; ++j)
        for (int n = 0; n < N; ++n)
            mulModAssign(qData[j][n], qInvModq[j], subMod(qData[j][n], qData[L][n], qVec[j]), qVec[j]);
    qData.pop_back();

    if (doDFT) DFT(rotGroupRev, qVec, qPsi);
}

Poly Poly::rotateBy(sizeType nrot, const vector<valueType> &qVec, const vector<valueType> &pVec) const {
    assert(status == COEFF_VALUE);
    if (nrot == 1) return Poly(*this);

    vector<vector<valueType>> qdata;
    sizeType dN = N << 1;
    sizeType qN = N << 2;
    for (int j = 0; j < L; ++j) {
        const auto &qdat = qData[j];
        valueType modular = qVec[j];
        vector<valueType> nqdat(N);
        for (int n = 0; n < N; ++n) {
            auto idx = n * nrot;
            auto pos = idx & dN - 1;
            bool flag = (pos >= N) ^ ((idx & qN - 1) >= dN);
            nqdat[min(pos, dN - pos)] = (qdat[n] && flag) ? modular - qdat[n] : qdat[n];
        }
        qdata.push_back(nqdat);
    }

    if (!K) return Poly(logN, qdata, COEFF_VALUE);
    vector<vector<valueType>> pdata;
    for (int i = 0; i < K; ++i) {
        const auto &pdat = pData[i];
        valueType modular = pVec[i];
        vector<valueType> npdat(N);
        for (int n = 0; n < N; ++n) {
            auto idx = n * nrot;
            auto pos = idx & dN - 1;
            bool flag = (pos >= N) ^ ((idx & qN - 1) >= dN);
            npdat[min(pos, dN - pos)] = (pdat[n] && flag) ? modular - pdat[n] : pdat[n];
        }
        pdata.push_back(npdat);
    }

    return Poly(logN, qdata, pdata, COEFF_VALUE);
}

void Poly::rotateByAndEqual(sizeType nslots, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    if (nslots == 1) return;
    *this = rotateBy(nslots, qVec, pVec);

}

void Poly::print(ostream &out, const vector<valueType> &qVec, const vector<valueType> &pVec) {
    out << '\n';
    out << "status: " << (POINT_VALUE ? "POINT_VALUE\n" : "COEFF_VALUE\n");
    out << "p:\n";
    for (int i = 0; i < K; ++i) {
        for (int n = 0; n < N; ++n)
            out << (pData[i][n] <= pVec[i] / 2 ? pData[i][n] : -(double) (pVec[i] - pData[i][n])) << ", ";
        out << '\n';
    }

    out << "q:\n";
    for (int j = 0; j < L; ++j) {
        for (int n = 0; n < N; ++n)
            out << (qData[j][n] <= qVec[j] / 2 ? qData[j][n] : -(double) (qVec[j] - qData[j][n])) << ", ";
        out << '\n';
    }
}

void Poly::setK(sizeType k) {
    K = k;
    pData.clear();
}


