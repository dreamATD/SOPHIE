//
// Created by 69029 on 4/3/2020.
//

#ifndef HOMENC_ENCODER_H
#define HOMENC_ENCODER_H


#include "Context.h"

class Encoder {
private:
    const Context *context;

    const sizeType logNh;
    const sizeType logN;
    const sizeType logM;
    const sizeType Nh;
    const sizeType N;
    const sizeType L;
    const valueType scale;

    const vector<sizeType> &revh;
    const vector<complex<double>> &ksiInvPow;
    const vector<sizeType> &rotGroup;
    const vector<sizeType> &rotGroupRev;

    const vector<valueType> &qVec;
    const vector<vector<valueType>> &qPsi;

public:
    explicit Encoder(const Context &_context);

    Plaintext encode(const vector<double> &mvec);
};


#endif //HOMENC_ENCODER_H
