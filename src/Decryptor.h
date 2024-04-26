//
// Created by 69029 on 4/3/2020.
//

#ifndef HOMENC_DECRYPTOR_H
#define HOMENC_DECRYPTOR_H


#include "Plaintext.h"
#include "Ciphertext.h"
#include "SecretKey.h"

class Decryptor {

private:
    const Context *context;

    const SecretKey &sk;

    const vector<valueType> &qVec;

public:
    Decryptor(const Context &_context, const SecretKey &sk);

    Plaintext decrypt(const Ciphertext &ct);

};


#endif //HOMENC_DECRYPTOR_H
