//
// Created by 69029 on 3/23/2020.
//

#ifndef SRC_PUBLICKEY_H
#define SRC_PUBLICKEY_H


#include "Poly.h"
#include "SecretKey.h"
#include "CommonHeader.h"

class Ciphertext;
class Encryptor;
class Evaluator;

class Key {

    friend class Ciphertext;
    friend class Encryptor;
    friend class Evaluator;

protected:
    const Context *context;

    Poly bx;
    Poly ax;

    Key();

    Key(const Context &_context, const SecretKey &sk, sizeType L, sizeType K);

    const Poly &getAx() const;

    const Poly &getBx() const;
};

class PublicKey: public Key {
public:
    PublicKey(const Context &_context, const SecretKey &sk);
};

class EvalKey: public Key {
public:
    EvalKey();

    EvalKey(const Context &_context, const SecretKey &sk);
};

class RotKey: public Key {
public:
    RotKey();

    RotKey(const Context &_context, const SecretKey &sk, sizeType nslots);
};

#endif //SRC_PUBLICKEY_H
