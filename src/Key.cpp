//
// Created by 69029 on 3/23/2020.
//

#include "Key.h"
#include "Poly.h"


Key::Key(): context(nullptr), bx(0, 0, 0), ax(0, 0, 0) {

}

Key::Key(const Context &_context, const SecretKey &sk, sizeType L, sizeType K = 0):
    context(&_context),
    bx(sk.sx), 
    ax(_context.logNh, L, K) {
    ax.sampleUniform(_context.qVec, _context.pVec);
    
    Poly ex(_context.logNh, L, K);
    ex.sampleGauss(_context.sigma, _context.qVec, _context.pVec);
    ex.DFT(_context.rotGroupRev, _context.qVec, _context.qPsi, _context.pVec, _context.pPsi);
    
    if (K) bx.mulAndEqual(ax, _context.qVec, _context.pVec);
    else bx.mulAndEqualProductToq(ax, _context.qVec);
    bx.subAndEqual2(ex, _context.qVec, _context.pVec);
}

const Poly &Key::getBx() const {
    return bx;
}

const Poly &Key::getAx() const {
    return ax;
}

PublicKey::PublicKey(const Context &_context, const SecretKey &sk):
    Key(_context, sk, _context.L) {
}

EvalKey::EvalKey(): Key() {

}

EvalKey::EvalKey(const Context &_context, const SecretKey &sk):
        Key(_context, sk, _context.L, _context.K) {
    Poly sxsx(sk.sx.mulAndProductToq(sk.sx, _context.qVec));

    sxsx.iDFT(_context.rotGroupRev, _context.qVec, _context.qInvPsi, _context.pVec, _context.pInvPsi);
    sxsx.imulAndEqual(_context.PModq, _context.qVec);
    sxsx.DFT(_context.rotGroupRev, _context.qVec, _context.qPsi, _context.pVec, _context.pPsi);

    sxsx.zeroExpand(_context.K);

    bx.addAndEqual(sxsx, _context.qVec, _context.pVec);
}

RotKey::RotKey(): Key() {

}

RotKey::RotKey(const Context &_context, const SecretKey &sk, sizeType nslots):
    Key(_context, sk, _context.L, _context.K) {
    Poly sxr(sk.sx);
    sxr.setK(0);
    sxr.iDFT(_context.rotGroupRev, _context.qVec, _context.qInvPsi);
    sxr.rotateByAndEqual(_context.rotGroup[nslots], _context.qVec);

    sxr.imulAndEqual(_context.PModq, _context.qVec);
    sxr.DFT(_context.rotGroupRev, _context.qVec, _context.qPsi);

    sxr.zeroExpand(_context.K);

    bx.addAndEqual(sxr, _context.qVec, _context.pVec);
}
