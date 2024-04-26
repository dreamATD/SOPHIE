//
// Created by 69029 on 4/3/2020.
//

#ifndef HOMENC_ENCRYPTOR_H
#define HOMENC_ENCRYPTOR_H


#include "Context.h"
#include "Ciphertext.h"
#include "Plaintext.h"
#include "Key.h"

class Encryptor {
private:
    const Context *context;

    const sizeType logNh;
    const sizeType logM;

    const sizeType L;

    const double sigma;

    const PublicKey pk;

    const vector<sizeType> &rotGroupRev;
    const vector<valueType> &qVec;
    const vector<vector<valueType>> &qPsi;

public:
    Encryptor(const Context &_context, const SecretKey &sk);

    Ciphertext encrypt(const Plaintext &m);

};


#endif //HOMENC_ENCRYPTOR_H
