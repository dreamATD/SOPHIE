//
// Created by 69029 on 4/3/2020.
//

#include "Evaluator.h"

Evaluator::Evaluator(const Context &_context):
        sk(nullptr),
        context(&_context),
        logNh(_context.logNh),
        logM(_context.logM),
        Nh(_context.Nh),
        K(_context.K),
        evk(),
        rotGroup(_context.rotGroup),
        rotGroupRev(_context.rotGroupRev),
        qVec(_context.qVec),
        pVec(_context.pVec),
        qPsi(_context.qPsi),
        pPsi(_context.pPsi),
        qInvModq(_context.qInvModq),
        qInvPsi(_context.qInvPsi),
        pInvPsi(_context.pInvPsi),
        qHatInvModq(_context.qHatInvModq),
        pHatInvModp(_context.pHatInvModp),
        qHatModp(_context.qHatModp),
        pHatModq(_context.pHatModq),
        PInvModq(_context.PInvModq) {
}

Evaluator::Evaluator(const Context &_context, const SecretKey &_sk):
        sk(&_sk),
        context(&_context),
        logNh(_context.logNh),
        logM(_context.logM),
        Nh(_context.Nh),
        K(_context.K),
        evk(_context, _sk),
        rotGroup(_context.rotGroup),
        rotGroupRev(_context.rotGroupRev),
        qVec(_context.qVec),
        pVec(_context.pVec),
        qPsi(_context.qPsi),
        pPsi(_context.pPsi),
        qInvModq(_context.qInvModq),
        qInvPsi(_context.qInvPsi),
        pInvPsi(_context.pInvPsi),
        qHatInvModq(_context.qHatInvModq),
        pHatInvModp(_context.pHatInvModp),
        qHatModp(_context.qHatModp),
        pHatModq(_context.pHatModq),
        PInvModq(_context.PInvModq) {
    assert(context == _sk.context);
    rtk.resize(Nh);
}

MulCiphertext Evaluator::mulWithoutRelin(const Ciphertext &a, const Ciphertext &b) {
    assert(context == a.context && context == b.context);

    Poly d0(a.bx.mul(b.bx, qVec));

    Poly tmp1(a.ax.mul(b.bx, qVec));
    Poly tmp2(a.bx.mul(b.ax, qVec));
    Poly d1(tmp1.add(tmp2, qVec));

    Poly d2(a.ax.mul(b.ax, qVec));

    return MulCiphertext(*context, move(d0), move(d1), move(d2));
}

Ciphertext Evaluator::relinearize(const MulCiphertext &a, const Key &key) {
    Poly d0(a.cx), d1(a.bx), d2(a.cx);
    sizeType L = a.ax.getL();

    d2.modUp(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pPsi, K, qHatInvModq[L - 1],
             qHatModp[L - 1]);

    Poly bx(d2.mul(key.bx, qVec, pVec));
    Poly ax(d2.mul(key.ax, qVec, pVec));

    bx.modDown(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pInvPsi, pHatInvModp, pHatModq, PInvModq);
    ax.modDown(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pInvPsi, pHatInvModp, pHatModq, PInvModq);

    bx.addAndEqual(d0, qVec);
    ax.addAndEqual(d1, qVec);

    return Ciphertext(*context, move(bx), move(ax));
}

Ciphertext Evaluator::relinearize(MulCiphertext &&a, const Key &key) {
    Poly &d0(a.cx), &d1(a.bx), &d2(a.ax);
    sizeType L = a.ax.getL();

    d2.modUp(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pPsi, K, qHatInvModq[L - 1],
             qHatModp[L - 1]);

    Poly bx(d2.mul(key.bx, qVec, pVec));
    Poly ax(d2.mul(key.ax, qVec, pVec));

    bx.modDown(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pInvPsi, pHatInvModp, pHatModq, PInvModq);
    ax.modDown(rotGroupRev, qVec, qPsi, qInvPsi, pVec, pInvPsi, pHatInvModp, pHatModq, PInvModq);

    bx.addAndEqual(d0, qVec);
    ax.addAndEqual(d1, qVec);

    return Ciphertext(*context, move(bx), move(ax));
}

Plaintext Evaluator::add(const Plaintext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Plaintext(*context, a.mx.add(b.mx, qVec));
}

void Evaluator::addAndEqual(Plaintext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.mx.addAndEqual(b.mx, qVec);
}

Plaintext Evaluator::sub(const Plaintext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Plaintext(*context, a.mx.sub(b.mx, qVec));
}

void Evaluator::subAndEqual(Plaintext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.mx.subAndEqual(b.mx, qVec);
}

Plaintext Evaluator::neg(const Plaintext &a) {
    assert(context == a.context);
    return Plaintext(*context, a.mx.neg(qVec));
}

void Evaluator::negAndEqual(Plaintext &res) {
    assert(context == res.context);
    res.mx.negAndEqual(qVec);
}

Plaintext Evaluator::mul(const Plaintext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Plaintext(*context, a.mx.mul(b.mx, qVec));
}

void Evaluator::mulAndEqual(Plaintext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.mx.mulAndEqual(b.mx, qVec);
}


Plaintext Evaluator::leftRotate(const Plaintext &a, sizeType nslots) {
    if (!nslots) return Plaintext(a);
    nslots = nslots % Nh;
    Poly mx(a.mx);
    mx.iDFT(rotGroupRev, qVec, qInvPsi);
    mx.rotateByAndEqual(rotGroup[nslots], qVec);
    mx.DFT(rotGroupRev, qVec, qPsi);
    return Plaintext(*context, mx);
}

void Evaluator::leftRotateAndEqual(Plaintext &res, sizeType nslots) {
    if (!nslots) return;
    nslots = nslots % Nh;
    Poly &mx = res.mx;
    mx.iDFT(rotGroupRev, qVec, qInvPsi);
    mx.rotateByAndEqual(rotGroup[nslots], qVec);
    mx.DFT(rotGroupRev, qVec, qPsi);
}

Plaintext Evaluator::rightRotate(const Plaintext &a, sizeType nslots) {
    if (!nslots) return Plaintext(a);
    nslots = nslots % Nh;
    return leftRotate(a, Nh - nslots);
}

void Evaluator::rightRotateAndEqual(Plaintext &res, sizeType nslots) {
    if (!nslots) return;
    nslots = nslots % Nh;
    leftRotateAndEqual(res, Nh - nslots);
}

Ciphertext Evaluator::add(const Ciphertext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.bx.add(b.mx, qVec), a.ax);
}

Ciphertext Evaluator::add(const Plaintext &a, const Ciphertext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.mx.add(b.bx, qVec), b.ax);
}

void Evaluator::addAndEqual(Ciphertext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.bx.addAndEqual(b.mx, qVec);
}

Ciphertext Evaluator::sub(const Ciphertext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.bx.sub(b.mx, qVec), a.ax);
}

Ciphertext Evaluator::sub(const Plaintext &a, const Ciphertext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.mx.sub(b.bx, qVec), b.ax);
}

void Evaluator::subAndEqual(Ciphertext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.bx.subAndEqual(b.mx, qVec);
}

void Evaluator::subAndEqual2(const Plaintext &a, Ciphertext &res) {
    assert(context == a.context && context == res.context);
    res.bx.subAndEqual2(a.mx, qVec);
}

Ciphertext Evaluator::mul(const Ciphertext &a, const Plaintext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.bx.mul(b.mx, qVec), a.ax.mul(b.mx, qVec));
}

void Evaluator::mulAndEqual(Ciphertext &res, const Plaintext &b) {
    assert(context == res.context && context == b.context);
    res.bx.mulAndEqual(b.mx, qVec);
    res.ax.mulAndEqual(b.mx, qVec);
}

void Evaluator::rescaleAndEqual(Plaintext &a) {
    a.mx.rescale(rotGroupRev, qVec, qPsi, qInvPsi, qInvModq[a.mx.getL() - 1]);
}

Ciphertext Evaluator::add(const Ciphertext &a, const Ciphertext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.bx.add(b.bx, qVec), a.ax.add(b.ax, qVec));
}

void Evaluator::addAndEqual(Ciphertext &res, const Ciphertext &b) {
    assert(context == res.context && context == b.context);

    res.bx.addAndEqual(b.bx, qVec);
    res.ax.addAndEqual(b.ax, qVec);
}

Ciphertext Evaluator::sub(const Ciphertext &a, const Ciphertext &b) {
    assert(context == a.context && context == b.context);
    return Ciphertext(*context, a.bx.sub(b.bx, qVec), a.ax.sub(b.ax, qVec));
}

void Evaluator::subAndEqual(Ciphertext &res, const Ciphertext &b) {
    assert(context == res.context && context == b.context);

    res.bx.subAndEqual(b.bx, qVec);
    res.ax.subAndEqual(b.ax, qVec);
}

Ciphertext Evaluator::neg(const Ciphertext &a) {
    assert(context == a.context);
    return Ciphertext(*context, a.bx.neg(qVec), a.ax.neg(qVec));
}

void Evaluator::negAndEqual(Ciphertext &res) {
    assert(context == res.context);

    res.bx.negAndEqual(qVec);
    res.ax.negAndEqual(qVec);
}

Ciphertext Evaluator::mul(const Ciphertext &a, const Ciphertext &b) {
    return relinearize(mulWithoutRelin(a, b), evk);
}

void Evaluator::mulAndEqual(Ciphertext &res, const Ciphertext &b) {
    res = mul(res, b);
}

Ciphertext Evaluator::leftRotate(const Ciphertext &a, sizeType nslots) {
    if (!nslots) return a;
    nslots = nslots % Nh;
    RotKey *key = rtk[nslots];
    if (key == nullptr)
        key = rtk[nslots] = new RotKey(*context, *sk, nslots);

    Poly bx(a.bx), ax(a.ax);

    bx.iDFT(rotGroupRev, qVec, qInvPsi);
    bx.rotateByAndEqual(rotGroup[nslots], qVec);
    bx.DFT(rotGroupRev, qVec, qPsi);

    ax.iDFT(rotGroupRev, qVec, qInvPsi);
    ax.rotateByAndEqual(rotGroup[nslots], qVec);
    ax.DFT(rotGroupRev, qVec, qPsi);

    return relinearize(MulCiphertext(*context, move(bx), Poly::ZERO(logNh, bx.getL()), move(ax)),
                       *key);
}

void Evaluator::leftRotateAndEqual(Ciphertext &res, sizeType nslots) {
    if (!nslots) return;
    nslots = nslots % Nh;
    RotKey *key = rtk[nslots];
    if (key == nullptr)
        key = rtk[nslots] = new RotKey(*context, *sk, nslots);

    res.bx.iDFT(rotGroupRev, qVec, qInvPsi);
    res.bx.rotateByAndEqual(rotGroup[nslots], qVec);
    res.bx.DFT(rotGroupRev, qVec, qPsi);

    res.ax.iDFT(rotGroupRev, qVec, qInvPsi);
    res.ax.rotateByAndEqual(rotGroup[nslots], qVec);
    res.ax.DFT(rotGroupRev, qVec, qPsi);
    res = relinearize(MulCiphertext(*context, move(res.bx), Poly::ZERO(logNh, res.bx.getL()), move(res.ax)),
                      *key);
}

Ciphertext Evaluator::rightRotate(const Ciphertext &a, sizeType nslots) {
    if (!nslots) return Ciphertext(a);
    return leftRotate(a, Nh - nslots);
}

void Evaluator::rightRotateAndEqual(Ciphertext &res, sizeType nslots) {
    if (!nslots) return;
    leftRotateAndEqual(res, Nh - nslots);
}

void Evaluator::rescaleAndEqual(Ciphertext &a) {
    a.bx.rescale(rotGroupRev, qVec, qPsi, qInvPsi, qInvModq[a.bx.getL() - 1]);
    a.ax.rescale(rotGroupRev, qVec, qPsi, qInvPsi, qInvModq[a.ax.getL() - 1]);
}


