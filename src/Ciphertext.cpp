//
// Created by 69029 on 3/23/2020.
//

#include <assert.h>

#include <utility>
#include "Ciphertext.h"
#include "Plaintext.h"
#include "SecretKey.h"
#include "Key.h"

Ciphertext::Ciphertext(const Context &_context):
        context(&_context),
        bx(_context.logNh, _context.L),
        ax(_context.logNh, _context.L) {
}

Ciphertext::Ciphertext(const Context &_context, Poly bx, Poly ax):
        context(&_context), bx(move(bx)), ax(move(ax)) {
}

Ciphertext::Ciphertext(Ciphertext &&other) noexcept:
        context(other.context), bx(move(other.bx)), ax(move(other.ax)) {

}

Ciphertext &Ciphertext::operator = (Ciphertext &&other) noexcept {
    if (this != &other) {
        context = other.context;
        bx = move(other.bx);
        ax = move(other.ax);
    }
    return *this;
}

MulCiphertext::MulCiphertext(const Context &_context, Poly &&cx, Poly &&bx, Poly &&ax):
    Ciphertext(_context, bx, ax),
    cx(cx) {

}
