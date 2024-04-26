//
// Created by 69029 on 4/10/2020.
//

#include <CommonHeader.h>
#include <Context.h>
#include <Encoder.h>
#include <Decoder.h>
#include <Plaintext.h>
#include "catch.hpp"
#include "Utils.h"
#include "TestUtils.h"

class ContextForEncoderTest : public Context {
public:
    ContextForEncoderTest(): Context(15, 15, 55) {}
};

/**
 * Test whether ifft(fft(mvec)) can retrive mvec.
 * */

TEST_CASE_METHOD(ContextForEncoderTest, "fft&ifft for real", "[fft&ifft]") {
    vector<double> dat(Nh);
    for (int i = 0; i < Nh; ++i)
        getRandomInRange(dat[i], (valueType) 0, (valueType) 20);

    vector<double> ndat = dat;
    fft(ndat, revh, ksiPow, rotGroup);
    ifft(ndat, 1, revh, ksiInvPow, rotGroup);

    printVec(cout, "original array", dat);
    printVec(cout, "transformed array", ndat);

    REQUIRE(equals(dat, ndat));
}

/**
 * Test whether decode(encode(mvec)) can retrive mvec.
 * */

TEST_CASE("Encode&Decode", "[scheme]") {
    cout << "\nTesting encoding and decoding ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 30);
    vector<double> mvec(Nh);
    for (int i = 0; i < Nh; ++i)
        getRandomInRange(mvec[i], (valueType) 1, (valueType) B);

    printVec(cout, "- message before encoding", mvec);

    Encoder encoder(context);
    Plaintext m = encoder.encode(mvec);

    cout << "- plaintext after encoding" << m.getMx();

    Decoder decoder(context);
    auto dmvec = decoder.decode(m);

    printVec(cout, "- message after decoding", dmvec);

    REQUIRE(equals(mvec, dmvec));
    cout << '\n';
}