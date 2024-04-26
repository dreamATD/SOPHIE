//
// Created by 69029 on 4/3/2020.
//

#ifndef HOMENC_EVALUATOR_H
#define HOMENC_EVALUATOR_H


#include <vector>
#include "Ciphertext.h"
#include "Plaintext.h"
#include "Key.h"

using namespace std;

class Evaluator {

private:
    const SecretKey *sk;

protected:
    const Context *context;

    const sizeType logNh;
    const sizeType logM;
    const sizeType Nh;

    const sizeType K;

    const EvalKey evk;
    vector<RotKey *> rtk;

    const vector<valueType> &qVec;
    const vector<valueType> &pVec;

    const vector<sizeType> &rotGroup;
    const vector<sizeType> &rotGroupRev;

    const vector<vector<valueType>> &qPsi;
    const vector<vector<valueType>> &pPsi;

    const vector<vector<valueType>> &qInvPsi;
    const vector<vector<valueType>> &pInvPsi;

    const vector<vector<valueType>> &qInvModq;

    const vector<vector<valueType>> &qHatInvModq;
    const vector<valueType> &pHatInvModp;

    const vector<vector<vector<valueType>>> &qHatModp;
    const vector<vector<valueType>> &pHatModq;

    const vector<valueType> &PInvModq;

    MulCiphertext mulWithoutRelin(const Ciphertext &a, const Ciphertext &b);

    Ciphertext relinearize(const MulCiphertext &a, const Key &key);

    Ciphertext relinearize(MulCiphertext &&a, const Key &key);

public:
    Evaluator(const Context &_context);

    Evaluator(const Context &_context, const SecretKey &_sk);

    /**
     * For Plaintext & Plaintext.
     * */

    Plaintext add(const Plaintext &a, const Plaintext &b);
    void addAndEqual(Plaintext &res, const Plaintext &b);

    Plaintext sub(const Plaintext &a, const Plaintext &b);
    void subAndEqual(Plaintext &res, const Plaintext &b);

    Plaintext neg(const Plaintext &a);
    void negAndEqual(Plaintext &res);

    Plaintext mul(const Plaintext &a, const Plaintext &b);
    void mulAndEqual(Plaintext &res, const Plaintext &b);

    Plaintext leftRotate(const Plaintext &a, sizeType nslots);
    void leftRotateAndEqual(Plaintext &res, sizeType nslots);

    Plaintext rightRotate(const Plaintext &a, sizeType nslots);
    void rightRotateAndEqual(Plaintext &res, sizeType nslots);

    void rescaleAndEqual(Plaintext &a);

    /**
     * For Plaintext & Ciphertext.
     * */

    Ciphertext add(const Ciphertext &a, const Plaintext &b);
    Ciphertext add(const Plaintext &a, const Ciphertext &b);
    void addAndEqual(Ciphertext &res, const Plaintext &b);

    Ciphertext sub(const Ciphertext &a, const Plaintext &b);
    Ciphertext sub(const Plaintext &a, const Ciphertext &b);
    void subAndEqual(Ciphertext &res, const Plaintext &b);
    void subAndEqual2(const Plaintext &a, Ciphertext &res);

    Ciphertext mul(const Ciphertext &a, const Plaintext &b);
    void mulAndEqual(Ciphertext &res, const Plaintext &b);

    /**
     * For ciphertext & ciphertext.
     * */

    Ciphertext add(const Ciphertext &a, const Ciphertext &b);
    void addAndEqual(Ciphertext &res, const Ciphertext &b);

    Ciphertext sub(const Ciphertext &a, const Ciphertext &b);
    void subAndEqual(Ciphertext &res, const Ciphertext &b);

    Ciphertext neg(const Ciphertext &a);
    void negAndEqual(Ciphertext &res);

    Ciphertext mul(const Ciphertext &a, const Ciphertext &b);
    void mulAndEqual(Ciphertext &res, const Ciphertext &b);

    Ciphertext leftRotate(const Ciphertext &a, sizeType nslots);
    void leftRotateAndEqual(Ciphertext &res, sizeType nslots);

    Ciphertext rightRotate(const Ciphertext &a, sizeType nslots);
    void rightRotateAndEqual(Ciphertext &res, sizeType nslots);

    void rescaleAndEqual(Ciphertext &a);
};


#endif //HOMENC_EVALUATOR_H
