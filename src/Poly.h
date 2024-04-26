//
// Created by 69029 on 3/22/2020.
//

#ifndef SRC_POLY_H
#define SRC_POLY_H

#include <cstdint>
#include <zconf.h>
#include <complex>
#include <algorithm>
#include "Utils.h"
#include "CommonHeader.h"
#define POINT_VALUE true
#define COEFF_VALUE false


using namespace std;

class Poly {

protected:

    bool status;

    const sizeType N;
    const sizeType logN;

    sizeType L;
    vector<vector<valueType>> qData;

    sizeType K;

    vector<vector<valueType>> pData;

public:
    Poly(sizeType logn, sizeType l, sizeType k = 0, bool sta = POINT_VALUE);

    Poly(sizeType logn, vector<vector<valueType>> qdata, bool sta = POINT_VALUE);

    Poly(sizeType logn, vector<vector<valueType>> qdata, vector<vector<valueType>> pdata, bool sta = POINT_VALUE);

    Poly(sizeType logn, sizeType l, const vector<double> &data, const vector<valueType> &qVec, bool sta = POINT_VALUE);

    Poly(sizeType logn, sizeType l, const vector<complex<double>> &data, const vector<valueType> &qVec,
         bool sta = POINT_VALUE);

    Poly(const Poly &other) = default;

    Poly(Poly &&other) = default;

    Poly &operator = (Poly &&other) noexcept;

    bool operator == (const Poly &other) const;

    friend ostream &operator << (ostream &out, const Poly &p) {
        out << '\n';
        printMat(out, "pData", p.pData);
        printMat(out, "qData", p.qData);
        return out;
    }

    static Poly ZERO(sizeType logn, sizeType l, sizeType k = 0);

    void setL(sizeType l);

    sizeType getL() const;

    void setK(sizeType k);

    vector<vector<valueType>> getPData() const;

    void setPData(const vector<vector<valueType>> &pData);

    vector<vector<valueType>> getQData() const;

    void setQData(const vector<vector<valueType>> &qData);

    const sizeType getN() const;

    void sampleHWT(sizeType h, const vector<valueType> &qVec, const vector<valueType> &pVec);

    void sampleUniform(const vector<valueType> &qVec, const vector<valueType> &pVec = {}, bool sta = POINT_VALUE);

    void sampleGauss(double sigma, const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    void sampleZO(const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    void DFT(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
             const vector<valueType> &pVec = {}, const vector<vector<valueType>> &pPsi = {});

    void iDFT(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec,
              const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec = {},
              const vector<vector<valueType>> &pPsiInv = {});

    Poly add(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {}) const;
    void addAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    Poly sub(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {}) const;
    void subAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {});
    void subAndEqual2(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    Poly neg(const vector<valueType> &qVec, const vector<valueType> &pVec = {}) const;
    void negAndEqual(const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    Poly mul(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {}) const;
    void mulAndEqual(const Poly &other, const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    Poly mulAndProductToq(const Poly &other, const vector<valueType> &qVec) const;

    void mulAndEqualProductToq(const Poly &other, const vector<valueType> &qVec);

    Poly imul(const vector<valueType> &qOther, const vector<valueType> &qVec, const vector<valueType> &pOther = {},
              const vector<valueType> &pVec = {});
    void imulAndEqual(const vector<valueType> &qOther, const vector<valueType> &qVec,
                      const vector<valueType> &pOther = {}, const vector<valueType> &pVec = {});

    vector<double> convertToDouble(valueType scale, valueType modular);

    void zeroExpand(sizeType k);

    void modUp(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
               const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec,
               const vector<vector<valueType>> &pPsi, sizeType k, const vector<valueType> &qHatInvModq,
               const vector<vector<valueType>> &qHatModp);

    void modDown(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
                 const vector<vector<valueType>> &qInvPsi, const vector<valueType> &pVec,
                 const vector<vector<valueType>> &pInvPsi, const vector<valueType> &pHatInvModp,
                 const vector<vector<valueType>> &pHatModq, const vector<valueType> &PInvModq);

    void rescale(const vector<sizeType> &rotGroupRev, const vector<valueType> &qVec, const vector<vector<valueType>> &qPsi,
                 const vector<vector<valueType>> &qInvPsi, const vector<valueType> &qInvModq, bool doDFT = true);

    Poly rotateBy(sizeType nrot, const vector<valueType> &qVec, const vector<valueType> &pVec = {}) const;

    void rotateByAndEqual(sizeType nslots, const vector<valueType> &qVec, const vector<valueType> &pVec = {});

    void print(ostream &out, const vector<valueType> &qVec, const vector<valueType> &pVec = {});
};


#endif //SRC_POLY_H