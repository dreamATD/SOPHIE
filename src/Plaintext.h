//
// Created by 69029 on 3/23/2020.
//

#ifndef SRC_PLAINTEXT_H
#define SRC_PLAINTEXT_H


#include "Poly.h"
#include "Context.h"
#include "CommonHeader.h"

class Ciphertext;
class Encryptor;
class Evaluator;

class Plaintext {

    friend Ciphertext;
    friend Encryptor;
    friend Evaluator;

protected:
    const Context *context;

    Poly mx;

public:
    explicit Plaintext(const Context &_context);

    Plaintext(const Context &_context, Poly _mx);

    const Poly &getMx() const;
};


#endif //SRC_PLAINTEXT_H
