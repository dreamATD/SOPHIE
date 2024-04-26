//
// Created by 69029 on 3/31/2020.
//

#ifndef HOMENC_TESTUTILS_H
#define HOMENC_TESTUTILS_H

#include "../src/Utils.h"

complex<double> power(complex<double> a, sizeType b);

bool equals(const complex<double> &a, const complex<double> &b);

bool equals(const complex<double> &a, double b);

bool equals(double a, double b);

bool equals(const vector<double> &a, const vector<double> &b);


#endif //HOMENC_TESTUTILS_H
