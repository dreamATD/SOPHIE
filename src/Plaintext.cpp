//
// Created by 69029 on 3/23/2020.
//

#include "Plaintext.h"

Plaintext::Plaintext(const Context &_context):
    context(&_context), mx(_context.logNh, _context.L) {
}

Plaintext::Plaintext(const Context &_context, Poly _mx):
    context(&_context), mx(move(_mx)) {
}

const Poly &Plaintext::getMx() const {
    return mx;
}



