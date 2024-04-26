//
// Created by 69029 on 4/3/2020.
//

#include "Decryptor.h"

Decryptor::Decryptor(const Context &_context, const SecretKey &sk):
        context(&_context),
        sk(sk),
        qVec(_context.qVec) {
    assert(context == sk.context);
}

Plaintext Decryptor::decrypt(const Ciphertext &ct) {
    Poly mx(ct.ax.mulAndProductToq(sk.sx, qVec));
    mx.addAndEqual(ct.bx, qVec);
    return Plaintext(*context, mx);
}
