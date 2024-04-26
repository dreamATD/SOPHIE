//
// Created by 69029 on 3/22/2020.
//

#ifndef SRC_CONTEXT_H
#define SRC_CONTEXT_H

#include "CommonHeader.h"

class SecretKey;
class Key;
class PublicKey;
class EvalKey;
class RotKey;
class Plaintext;
class Ciphertext;
class Encoder;
class Decoder;
class Encryptor;
class Decryptor;
class Evaluator;

class Context {
    friend SecretKey;
    friend Key;
    friend PublicKey;
    friend EvalKey;
    friend RotKey;
    friend Plaintext;
    friend Ciphertext;
    friend Encoder;
    friend Decoder;
    friend Encryptor;
    friend Decryptor;
    friend Evaluator;

protected:
    const sizeType logN;
    const sizeType logNh;
    const sizeType logM;
    const sizeType N;
    const sizeType M;
    const sizeType Nh;
    const sizeType K;
    const sizeType L;

    const sizeType Q0_BIT_SIZE;
    const sizeType logp;

    const sizeType h;
    const double sigma;

    valueType scale;

    vector<complex<double>> ksiPow;
    vector<complex<double>> ksiInvPow;

    vector<sizeType> rotGroup;
    vector<sizeType> rotGroupRev;

    vector<sizeType> revh;

    vector<valueType> pVec;
    vector<valueType> qVec;

    vector<vector<valueType>> pModq;
    vector<vector<valueType>> qModp;

    vector<vector<valueType>> pInvModq;
    vector<vector<valueType>> qInvModp;

    vector<vector<valueType>> qInvModq;

    vector<valueType> PModq;
    vector<valueType> PInvModq;
    vector<vector<valueType>> QModp;

    vector<vector<valueType>> pHatModq;
    vector<vector<vector<valueType>>> qHatModp;

    vector<valueType> pHatInvModp;
    vector<vector<valueType>> qHatInvModq;

    vector<valueType> pRoots;
    vector<valueType> qRoots;

    vector<vector<valueType>> pPsi;
    vector<vector<valueType>> qPsi;
    vector<vector<valueType>> pInvPsi;
    vector<vector<valueType>> qInvPsi;

public:
    Context(sizeType logn, sizeType lev, sizeType lgp, sizeType qSize0 = 61, sizeType _h = 64, double _sigma = 3.2);

    vector<vector<valueType>> getQInvPsi() const;

    vector<valueType> getPVec() const;

    vector<valueType> getQVec() const;
};
#endif //SRC_CONTEXT_H
