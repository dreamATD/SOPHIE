//
// Created by 69029 on 4/3/2020.
//

#include "Decoder.h"
#include "Poly.h"
#include "Plaintext.h"
#include "Utils.h"

Decoder::Decoder(const Context &_context):
    context(&_context),
    logNh(_context.logNh),
    logN(_context.logN),
    logM(_context.logM),
    Nh(_context.Nh),
    N(_context.N),
    L(_context.L),
    scale(_context.scale),
    revh(_context.revh),
    ksiPow(_context.ksiPow),
    rotGroup(_context.rotGroup),
    rotGroupRev(_context.rotGroupRev),
    qVec(_context.qVec),
    qInvPsi(_context.qInvPsi)
    {

}

vector<double> Decoder::decode(const Plaintext &m) {
    Poly mx(m.getMx());
    mx.setL(1);
    mx.iDFT(rotGroupRev, qVec, qInvPsi);
    auto res = mx.convertToDouble(scale, qVec[0]);
    fft(res, revh, ksiPow, rotGroup);
    return res;
}


