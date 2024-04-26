//
// Created by 69029 on 4/16/2020.
//

#include <CommonHeader.h>
#include <Context.h>
#include <Utils.h>
#include <Encoder.h>
#include <Decoder.h>
#include <Encryptor.h>
#include <Decryptor.h>
#include "catch.hpp"
#include "TestUtils.h"

TEST_CASE("Encryptor&Decryptor", "[scheme]") {
    cout << "\nTesting encryption and decryption ...\n";
    const sizeType logN = 15;
    const sizeType logNh = logN - 1;
    const sizeType logM = logN + 1;
    const sizeType N = 1 << logN;
    const sizeType Nh = 1 << logNh;
    const sizeType M = 1 << logM;
    const valueType B = 10;

    Context context(logN, 15, 55, 61, 64);
    vector<double> mvec(Nh);
    double err;
    for (int i = 0; i < Nh; ++i)
        getRandomInRange(mvec[i], (valueType) 0, (valueType) B);

    Encoder encoder(context);
    Decoder decoder(context);
    SecretKey sk(context);
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    Plaintext p = encoder.encode(mvec);
    Ciphertext c = encryptor.encrypt(p);
    p = decryptor.decrypt(c);
    auto dmvec = decoder.decode(p);

    err = 0;
    for (int n = 0; n < Nh; ++n)
        err = max(err, fabs(mvec[n] - dmvec[n]));

    printVec(cout, "- message before encryption", mvec);
    printVec(cout, "- message after decryption", dmvec);
    cout << "- max error: " << err << '\n';

    REQUIRE(equals(mvec, dmvec));
    cout << '\n';
}
