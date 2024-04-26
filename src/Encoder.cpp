//
// Created by 69029 on 4/3/2020.
//

#include "Encoder.h"
#include "Poly.h"
#include "Plaintext.h"
#include "Utils.h"

Encoder::Encoder(const Context &_context):
        context(&_context),
        logNh(_context.logNh),
        logN(_context.logN),
        logM(_context.logM),
        Nh(_context.Nh),
        N(_context.N),
        scale(_context.scale),
        L(_context.L),
        revh (_context.revh),
        ksiInvPow(_context.ksiInvPow),
        rotGroup(_context.rotGroup),
        rotGroupRev(_context.rotGroupRev),
        qVec(_context.qVec),
        qPsi(_context.qPsi) {
}

Plaintext Encoder::encode(const vector<double> &mvec) {
    vector<double> mvec2 = mvec;

    ifft(mvec2, scale, revh, ksiInvPow, rotGroup);

    Poly mx(logNh, L, mvec2, qVec);
    mx.DFT(rotGroupRev, qVec, qPsi);

    return Plaintext(*context, mx);
}
