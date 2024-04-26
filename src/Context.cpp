//
// Created by 69029 on 3/22/2020.
//

#include <complex>
#include "Context.h"
#include "Utils.h"

// TODO: order primes.

Context::Context(sizeType logn, sizeType lev, sizeType lgp, sizeType qSize0, sizeType _h, double _sigma):
        logN(logn), logNh(logn - 1), logM(logn + 1),
        N(( (sizeType) 1 ) << logn), M(( (sizeType) 1 ) << logn + 1), Nh(( (sizeType) 1 ) << logn - 1),
        K(lev), L(lev),
        logp(lgp), scale(((valueType) 1) << lgp),
        Q0_BIT_SIZE(qSize0), h(_h), sigma(_sigma),
        ksiPow(M + 1), ksiInvPow(M + 1), rotGroup(Nh), rotGroupRev(Nh),
        revh(Nh),
        pVec(L), qVec(K),
        pModq(K, vector<valueType>(L)), qModp(L, vector<valueType>(K)),
        pInvModq(K, vector<valueType>(L)), qInvModp(L, vector<valueType>(K)),
        qInvModq(L, vector<valueType>(L)),
        PModq(L), PInvModq(L),QModp(L, vector<valueType>(K)),
        pHatModq(K, vector<valueType> (L)),
        qHatModp(L, vector<vector<valueType>>(L, vector<valueType>(K))),
        pHatInvModp(K), qHatInvModq(L, vector<valueType>(L)),
        pRoots(K), qRoots(L),
        pPsi(K, vector<valueType>(M + 1)),
        qPsi(L, vector<valueType>(M + 1)),
        pInvPsi(K, vector<valueType>(M + 1)),
        qInvPsi(L, vector<valueType>(M + 1)) {

    int bnd = 0;
    while (true) {
        ++bnd;

        valueType prime = (((valueType) 1) << Q0_BIT_SIZE) + (valueType) bnd * M + 1;
        if (isPrime(prime)) {
            qVec[0] = prime;
            break;
        }
    }

    bnd = 0;

    for (int j = 1; j < L; ++j) {

        while (true) {
            ++bnd;

            valueType prime1 = (((valueType) 1) << logp) + (valueType) bnd * M + 1;
            if (isPrime(prime1)) {
                qVec[j] = prime1;
                break;
            }

            valueType prime2 = (((valueType) 1) << logp) - (valueType) bnd * M + 1;
            if (isPrime(prime2)) {
                qVec[j] = prime2;
                break;
            }
        }
    }

    for (int i = 0; i < K; ++i) {

        while (true) {
            ++bnd;

            valueType prime1 = (((valueType) 1) << logp) + (valueType) bnd * M + 1;
            if (isPrime(prime1)) {
                pVec[i] = prime1;
                break;
            }

            valueType prime2 = (((valueType) 1) << logp) - (valueType) bnd * M + 1;
            if (isPrime(prime2)) {
                pVec[i] = prime2;
                break;
            }
        }
    }

    for (int i = 0; i < K; ++i)
        for (int j = 0; j < L; ++j) {
            modAssign(pModq[i][j], pVec[i], qVec[j]);
            pInvModq[i][j] = power(pVec[i], qVec[j] - 2, qVec[j]);
        }

    for (int j = 0; j < L; ++j)
        for (int i = 0; i < K; ++i) {
            modAssign(qModp[j][i], qVec[j], pVec[i]);
            qInvModp[j][i] = power(qVec[j], pVec[i] - 2, pVec[i]);
        }

    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j)
            if (i == j) qInvModq[i][j] = 0;
            else qInvModq[i][j] = power(qVec[i], qVec[j] - 2, qVec[j]);

    for (int j = 0; j < L; ++j) {
        valueType &ocur = PModq[j], &icur = PInvModq[j];
        ocur = 1;
        icur = 1;
        for (int i = 0; i < K; ++i) {
            mulModAssign(ocur, ocur, pModq[i][j], qVec[j]);
            mulModAssign(icur, icur, pInvModq[i][j], qVec[j]);
        }
    }

    for (int i = 0; i < K; ++i) {
        QModp[0][i] = qModp[0][i];
        for (int j = 1; j < L; ++j) {
            mulModAssign(QModp[j][i], QModp[j - 1][i], qModp[j][i], pVec[i]);
        }
    }

    for (int i = 0; i < K; ++i)
        for (int j = 0; j < L; ++j)
            mulModAssign(pHatModq[i][j], PModq[j], pInvModq[i][j], qVec[j]);

    for (int l = 0; l < L; ++l)
        for (int j = 0; j <= l; ++j)
            for (int i = 0; i < K; ++i)
                mulModAssign(qHatModp[l][j][i], QModp[l][i], qInvModp[j][i], pVec[i]);

    for (int i = 0; i < K; ++i) {
        valueType cur = 1;
        for (int j = 0; j < K; ++j)
            if (i != j)
                mulModAssign(cur, cur, pVec[j], pVec[i]);

        pHatInvModp[i] = power(cur, pVec[i] - 2, pVec[i]);
    }

    for (int l = 0; l < L; ++l) {
        for (int j = 0; j <= l; ++j) {
            valueType cur = 1;
            for (int i = 0; i <= l; ++i)
                if (j != i)
                    mulModAssign(cur, cur, qVec[i], qVec[j]);

            qHatInvModq[l][j] = power(cur, qVec[j] - 2, qVec[j]);
        }
    }

    for (int i = 0; i <= M; ++i) {
        double theta = 2.0 * M_PI * i / M;
        ksiPow[i].real(cos(theta));
        ksiPow[i].imag(sin(theta));
        ksiInvPow[i].real(cos(-theta));
        ksiInvPow[i].imag(sin(-theta));
    }

    revh[0] = 0;
    for (int i = 1; i < Nh; ++i)
        revh [i] = revh[i >> 1] >> 1 ^ (i & 1) << logNh - 1;

    sizeType five = 1;
    for (int i = 0; i < Nh; ++i) {
        rotGroup[i] = five;
        rotGroupRev[i] = five;
        five = five * 5 % M;
    }
    arrayBitReverse(rotGroupRev, revh);

    for (int i = 0; i < K; ++i) {
        valueType r = power(findPrimitives(pVec[i]), (pVec[i] - 1) / M, pVec[i]);
        pRoots[i] = r;
        pPsi[i][0] = pInvPsi[i][0] = 1;
        pPsi[i][M] = pInvPsi[i][M] = 1;
        for (int m = 1; m < M; ++m) {
            mulModAssign(pPsi[i][m], pPsi[i][m - 1], r, pVec[i]);
            pInvPsi[i][M - m] = pPsi[i][m];
        }
    }

    for (int j = 0; j < L; ++j) {
        valueType r = power(findPrimitives(qVec[j]), (qVec[j] - 1) / M, qVec[j]);
        qRoots[j] = r;
        qPsi[j][0] = qInvPsi[j][0] = 1;
        qPsi[j][M] = qInvPsi[j][M] = 1;
        for (int m = 1; m < M; ++m) {
            mulModAssign(qPsi[j][m], qPsi[j][m - 1], r, qVec[j]);
            qInvPsi[j][M - m] = qPsi[j][m];
        }
    }
}

vector<valueType> Context::getQVec() const {
    return qVec;
}

vector<valueType> Context::getPVec() const {
    return pVec;
}

vector<vector<valueType>> Context::getQInvPsi() const {
    return qInvPsi;
}

