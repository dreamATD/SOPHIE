//
// Created by 69029 on 4/3/2020.
//

#include "Encryptor.h"

Encryptor::Encryptor(const Context &_context, const SecretKey &sk):
        context(&_context),
        logNh(_context.logNh),
        logM(_context.logM),
        L(_context.L),
        sigma(_context.sigma),
        pk(_context, sk),
        rotGroupRev(_context.rotGroupRev),
        qVec(_context.qVec),
        qPsi(_context.qPsi) {
    assert(context == sk.context);
}

Ciphertext Encryptor::encrypt(const Plaintext &m) {
    assert(context == m.context);

    Poly mx(m.getMx());
    Poly v(logNh, L);
    v.sampleZO(qVec);
    v.DFT(rotGroupRev, qVec, qPsi);

    Poly e0(logNh, L);
    e0.sampleGauss(sigma, qVec);
    e0.DFT(rotGroupRev, qVec, qPsi);

    Poly e1(logNh, L);
    e1.sampleGauss(sigma, qVec);
    e1.DFT(rotGroupRev, qVec, qPsi);

    Poly bx(v.mul(pk.bx, qVec)), ax(v.mul(pk.ax, qVec));

    bx.addAndEqual(mx, qVec);
    bx.addAndEqual(e0, qVec);

    ax.addAndEqual(e1, qVec);
    return Ciphertext(*context, bx, ax);
}
