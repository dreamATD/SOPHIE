//
// Created by 69029 on 4/9/2020.
//

#include "TestUtils.h"

complex<double> power(complex<double> a, sizeType b) {
    complex<double> res(1, 0);
    for (; b; b >>= 1, a = a * a) if (b & 1) res = res * a;
    return res;
}

// TODO: eps needs fixing.
bool equals(double a, double b) {
    static const double eps = 1e-6;
    return fabs(a - b) < eps;
}

bool equals(const complex<double> &a, const complex<double> &b) {
    static const double eps = 1e-7;
    return equals(a.real(), b.real()) && equals(a.imag(), b.imag());
}

bool equals(const complex<double> &a, double b) {
    return equals(a.real(), b) && equals(a.imag(), 0);
}

bool equals(const vector<double> &a, const vector<double> &b) {
    if (a.size() != b.size()) return false;
    auto N = a.size();
    for (int i = 0; i < N; ++i)
        if (!equals(a[i], b[i])) {
            return false;
        }
    return true;
}
