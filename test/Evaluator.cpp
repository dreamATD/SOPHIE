//
// Created by 69029 on 4/11/2020.
//

#include <CommonHeader.h>
#include <Context.h>
#include <Utils.h>
#include <Encoder.h>
#include <Plaintext.h>
#include <Evaluator.h>
#include <Decoder.h>
#include <Decryptor.h>
#include <Encryptor.h>
#include "catch.hpp"
#include "TestUtils.h"

TEST_CASE("Plaintext OP Plaintext", "[scheme]") {
    cout << "\nTesting evaluation between plaintext and plaintext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);
    Evaluator evaluator(context);

    Plaintext ma = evaluator.add(m1, m2);
    Plaintext ms = evaluator.sub(m1, m2);
    Plaintext mm = evaluator.mul(m1, m2);

    evaluator.rescaleAndEqual(mm);

    auto mvec4a = decoder.decode(ma);
    auto mvec4s = decoder.decode(ms);
    auto mvec4m = decoder.decode(mm);

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Ciphertext OP Plaintext", "[scheme]") {
    cout << "\nTesting evaluation between ciphertext and plaintext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Ciphertext c1 = encryptor.encrypt(m1);

    Evaluator evaluator(context);

    Ciphertext ca = evaluator.add(c1, m2);
    Ciphertext cs = evaluator.sub(c1, m2);
    Ciphertext cm = evaluator.mul(c1, m2);
    evaluator.rescaleAndEqual(cm);

    auto mvec4a = decoder.decode(decryptor.decrypt(ca));
    auto mvec4s = decoder.decode(decryptor.decrypt(cs));
    auto mvec4m = decoder.decode(decryptor.decrypt(cm));

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Ciphertext OP Ciphertext", "[scheme]") {
    cout << "\nTesting evaluation between ciphertext and ciphertext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Ciphertext c1 = encryptor.encrypt(m1);
    Ciphertext c2 = encryptor.encrypt(m2);

    Evaluator evaluator(context, sk);

    Ciphertext ca = evaluator.add(c1, c2);
    Ciphertext cs = evaluator.sub(c1, c2);
    Ciphertext cm = evaluator.mul(c1, c2);
    evaluator.rescaleAndEqual(cm);

    auto mvec4a = decoder.decode(decryptor.decrypt(ca));
    auto mvec4s = decoder.decode(decryptor.decrypt(cs));
    auto mvec4m = decoder.decode(decryptor.decrypt(cm));

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Plaintext OP&EQual Plaintext", "[scheme]") {
    cout << "\nTesting plaintext OP&EQUAL plaintext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);
    Evaluator evaluator(context);

    Plaintext ma(m1);
    Plaintext ms(m1);
    Plaintext mm(m1);
    evaluator.addAndEqual(ma, m2);
    evaluator.subAndEqual(ms, m2);
    evaluator.mulAndEqual(mm, m2);

    evaluator.rescaleAndEqual(mm);

    auto mvec4a = decoder.decode(ma);
    auto mvec4s = decoder.decode(ms);
    auto mvec4m = decoder.decode(mm);

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Ciphertext OP&EQUAL Plaintext", "[scheme]") {
    cout << "\nTesting ciphertext OP&EQUAL plaintext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Ciphertext c1 = encryptor.encrypt(m1);

    Evaluator evaluator(context);

    Ciphertext ca(c1);
    Ciphertext cs(c1);
    Ciphertext cm(c1);
    evaluator.addAndEqual(ca, m2);
    evaluator.subAndEqual(cs, m2);
    evaluator.mulAndEqual(cm, m2);
    evaluator.rescaleAndEqual(cm);

    auto mvec4a = decoder.decode(decryptor.decrypt(ca));
    auto mvec4s = decoder.decode(decryptor.decrypt(cs));
    auto mvec4m = decoder.decode(decryptor.decrypt(cm));

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Ciphertext OP&EQUAL Ciphertext", "[scheme]") {
    cout << "\nTesting ciphertext OP&EQUAL ciphertext ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec1(Nh), mvec2(Nh), ans4a(Nh), ans4s(Nh), ans4m(Nh);

    for (int i = 0; i < Nh; ++i) {
        getRandomInRange(mvec1[i], (valueType) 0, (valueType) B);
        getRandomInRange(mvec2[i], (valueType) 0, (valueType) B);
        ans4a[i] = mvec1[i] + mvec2[i];
        ans4s[i] = mvec1[i] - mvec2[i];
        ans4m[i] = mvec1[i] * mvec2[i];
    }

    printVec(cout, "- message 1", mvec1);
    printVec(cout, "- message 2", mvec2);

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m1 = encoder.encode(mvec1);
    Plaintext m2 = encoder.encode(mvec2);

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Ciphertext c1 = encryptor.encrypt(m1);
    Ciphertext c2 = encryptor.encrypt(m2);

    Evaluator evaluator(context, sk);

    Ciphertext ca(c1);
    Ciphertext cs(c1);
    Ciphertext cm(c1);
    evaluator.addAndEqual(ca, c2);
    evaluator.subAndEqual(cs, c2);
    evaluator.mulAndEqual(cm, c2);
    evaluator.rescaleAndEqual(cm);

    auto mvec4a = decoder.decode(decryptor.decrypt(ca));
    auto mvec4s = decoder.decode(decryptor.decrypt(cs));
    auto mvec4m = decoder.decode(decryptor.decrypt(cm));

    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);
    printVec(cout, "- result for sub", mvec4s);
    printVec(cout, "- answer for sub", ans4s);
    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);

    REQUIRE(equals(mvec4a, ans4a));
    REQUIRE(equals(mvec4s, ans4s));
    REQUIRE(equals(mvec4m, ans4m));

    cout << '\n';
}

TEST_CASE("Plaintext rotation", "[scheme]") {
    cout << "\nTesting plaintext rotation ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55);
    vector<double> mvec(Nh), ans(Nh);

    for (int i = 0; i < Nh; ++i)
        getRandomInRange(mvec[i], (valueType) 0, (valueType) B);

    sizeType nslots = 1;
    getRandomInRange(nslots, (sizeType) 0, Nh);
    rotate_copy(mvec.begin(), mvec.begin() + nslots, mvec.end(), ans.begin());

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m = encoder.encode(mvec);

    Evaluator evaluator(context);
    Plaintext lm = evaluator.leftRotate(m, nslots);
    Plaintext rm = evaluator.rightRotate(m, Nh - nslots);

    auto lmvec = decoder.decode(lm);
    auto rmvec = decoder.decode(rm);

    cout << "number of slots: " << nslots << '\n';
    printVec(cout, "- original          message", mvec);
    printVec(cout, "- result for  left rotation", lmvec);
    printVec(cout, "- result for right rotation", rmvec);
    printVec(cout, "- answer for       rotation", ans);

    REQUIRE(equals(lmvec, ans));
    REQUIRE(equals(rmvec, ans));


    cout << '\n';
}

TEST_CASE("Ciphertext rotation", "[scheme]") {
    cout << "\nTesting Ciphertext rotation ...\n";
    const sizeType logN = 16;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 20;

    Context context(logN, 15, 55);
    vector<double> mvec(Nh), ans(Nh);

    for (int i = 0; i < Nh; ++i)
        getRandomInRange(mvec[i], (valueType) 0, (valueType) B);

    sizeType nslots = 1;
    getRandomInRange(nslots, (sizeType) 0, Nh);
    rotate_copy(mvec.begin(), mvec.begin() + nslots, mvec.end(), ans.begin());

    Encoder encoder(context);
    Decoder decoder(context);
    Plaintext m = encoder.encode(mvec);

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Ciphertext c = encryptor.encrypt(m);

    Evaluator evaluator(context, sk);

    clock_t b1 = clock();
    Ciphertext lc = evaluator.leftRotate(c, nslots);
    clock_t b2 = clock();
    Ciphertext rc = evaluator.rightRotate(c, Nh - nslots);

    auto lm = decryptor.decrypt(lc);
    auto rm = decryptor.decrypt(rc);

    auto lmvec = decoder.decode(lm);
    auto rmvec = decoder.decode(rm);

    cout << "number of slots: " << nslots << '\n';
    printVec(cout, "- original          message", mvec);
    printVec(cout, "- result for  left rotation", lmvec);
    printVec(cout, "- result for right rotation", rmvec);
    printVec(cout, "- answer for       rotation", ans);

    REQUIRE(equals(lmvec, ans));
    REQUIRE(equals(rmvec, ans));

    cout << "Rot time: " << (b2 - b1) / 1000.0 << endl;

    cout << '\n';
}

TEST_CASE("Multi-Ciphertext OP", "[scheme]") {
    cout << "\nTesting evaluation among ciphertexts ...\n";

    const sizeType logN = 16;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 100;
    const valueType L = 15;

    Context context(logN, L, 55);
    vector<double> mvec[L];
    vector<double> ans4m(Nh), ans4a(Nh);
    for (int j = 0; j < Nh; ++j) {
        ans4m[j] = 1;
        ans4a[j] = 0;
    }

    for (int j = 0; j < L; ++j) {
        mvec[j].resize(Nh);
        for (int i = 0; i < Nh; ++i) {
            getRandomInRange(mvec[j][i], (valueType) 0, (valueType) B);
            mvec[j][i] /= 60;
            ans4m[i] *= mvec[j][i];
            ans4a[i] += mvec[j][i];
        }
    }

    Encoder encoder(context);
    Decoder decoder(context);
    vector<Plaintext> m;
    for (int j = 0; j < L; ++j)
        m.emplace_back(encoder.encode(mvec[j]));

    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    vector<Ciphertext> c;

    clock_t b1 = clock();
    encryptor.encrypt(m[0]);
    clock_t b2 = clock();

    for (int j = 0; j < L; ++j)
        c.emplace_back(encryptor.encrypt(m[j]));

    Evaluator evaluator(context, sk);

    clock_t b3 = clock();
    Ciphertext ca = evaluator.add(c[0], c[1]);
    for (int j = 2; j < L; ++j)
        evaluator.addAndEqual(ca, c[j]);
    clock_t b4 = clock();

    Ciphertext cm = evaluator.mul(c[0], c[1]);
    evaluator.rescaleAndEqual(cm);

    for (int j = 2; j < L; ++j) {
        evaluator.mulAndEqual(cm, c[j]);
        evaluator.rescaleAndEqual(cm);
    }
    clock_t b5 = clock();

    Decryptor decryptor(context, sk);

    clock_t b6 = clock();
    auto ma = decryptor.decrypt(ca);
    clock_t b7 = clock();

    auto mm = decryptor.decrypt(cm);

    auto mvec4m = decoder.decode(mm);
    auto mvec4a = decoder.decode(ma);

    printVec(cout, "- result for mul", mvec4m);
    printVec(cout, "- answer for mul", ans4m);
    printVec(cout, "- result for add", mvec4a);
    printVec(cout, "- answer for add", ans4a);

    REQUIRE(equals(mvec4m, ans4m));
    REQUIRE(equals(mvec4a, ans4a));

    cout << "Enc time: " << (b2 - b1) / 1000.0 << endl;
    cout << "Dec time: " << (b7 - b6) / 1000.0 << endl;
    cout << "Add time: " << (b4 - b3) / 1000.0 << endl;
    cout << "Mul time: " << (b5 - b4) / 1000.0 << endl;
}