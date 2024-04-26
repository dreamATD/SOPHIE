//
// Created by 69029 on 3/22/2020.
//

#include "SecretKey.h"

SecretKey::SecretKey(const Context &_context):
    context(&_context),
    sx(_context.logNh, _context.L, _context.K) {
    sx.sampleHWT(_context.h, _context.qVec, _context.pVec);
    sx.DFT(_context.rotGroupRev, _context.qVec, _context.qPsi, _context.pVec, _context.pPsi);
}

const Poly &SecretKey::getSx() const {
    return sx;
}
