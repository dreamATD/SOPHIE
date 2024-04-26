//
// Created by 69029 on 3/22/2020.
//

#ifndef SRC_SECRETKEY_H
#define SRC_SECRETKEY_H


#include <cstdint>
#include "Poly.h"
#include "Context.h"
#include "CommonHeader.h"

class Key;
class EvalKey;
class ConjKey;
class RotKey;
class Ciphertext;
class Encryptor;
class Decryptor;
class Evaluator;

class SecretKey {

    friend Key;
    friend EvalKey;
    friend ConjKey;
    friend RotKey;
    friend Ciphertext;
    friend Encryptor;
    friend Decryptor;
    friend Evaluator;

private:
    const Context *context;

public:
    const Poly &getSx() const;

private:
    Poly sx;

public:
    SecretKey(const Context &_context);


};


#endif //SRC_SECRETKEY_H
