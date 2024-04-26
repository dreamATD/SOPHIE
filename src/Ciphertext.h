//
// Created by 69029 on 3/23/2020.
//

#ifndef SRC_CIPHERTEXT_H
#define SRC_CIPHERTEXT_H


#include "Poly.h"
#include "Context.h"
#include "CommonHeader.h"
#include "Key.h"


class Encryptor;
class Decryptor;
class Evaluator;

class Ciphertext {
    friend Encryptor;
    friend Decryptor;
    friend Evaluator;

protected:
    const Context *context;

    Poly bx;
    Poly ax;

public:
    Ciphertext(const Context &_context);

    Ciphertext(const Context &_context, Poly bx, Poly ax);

    Ciphertext(const Ciphertext &other) = default;

    Ciphertext(Ciphertext &&other) noexcept;

    Ciphertext &operator = (Ciphertext &&other) noexcept;
};

class MulCiphertext : public Ciphertext {
    friend Evaluator;

protected:
    Poly cx;

public:
    MulCiphertext(const Context &_context, Poly &&cx, Poly &&bx, Poly &&ax);

};


#endif //SRC_CIPHERTEXT_H
