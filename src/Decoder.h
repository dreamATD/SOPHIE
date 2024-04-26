//
// Created by 69029 on 4/3/2020.
//

#include "Context.h"

#ifndef HOMENC_DECODER_H
#define HOMENC_DECIDER_H


class Decoder {
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

    const vector<complex<double>> & ksiPow;

    const vector<sizeType> &rotGroup;
    const vector<sizeType> &rotGroupRev;

    const vector<valueType> &qVec;
    const vector<vector<valueType>> &qInvPsi;

public:
    Decoder(const Context &_context);

    vector<double> decode(const Plaintext &m);
};


#endif //HOMENC_DECODER_H
