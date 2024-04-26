//
// Created by 69029 on 3/22/2020.
//

#include "Utils.h"
#include "CommonHeader.h"


bool millerTest(valueType d, valueType n) {
    valueType a;
    getRandomInRange(a, (valueType) 2, n - 2);

    valueType x = power(a, d, n);

    if (x == 1 || x == n - 1) return true;

    while (d != n - 1) {
        mulModAssign(x, x, x, n);
        d <<= 1u;

        if (x == 1) return false;
        if (x == n - 1) return true;
    }

    return false;
}

bool isPrime(valueType n) {
    int k = 200;
    if (n <= 1 || n == 4) return false;
    if (n <= 3) return true;

    valueType d = n - 1;
    while ((d & 1u) == 0) d >>= 1u;

    while (k--)
        if (!millerTest(d, n))
            return false;

    return true;
}

vector<valueType> factor(valueType x) {
    vector<valueType> res;

    if (!(x & 1u)) {
        res.push_back(2);
        while (!(x & 1u)) x >>= 1;
    }

    valueType sq = round(sqrt(x));
    for (valueType i = 3; i <= sq; ++i) {
        if (x % i == 0) res.push_back(i);
        while (x % i == 0) x /= i;
    }

    if (x > 2) res.push_back(x);
    return res;
}

valueType findPrimitives(valueType modular) {
    valueType phi = modular - 1;
    auto fac = factor(phi);

    for (valueType j = 2; j < modular; ++j) {
        bool flag = true;

        for (auto x: fac)
            if (power(j, phi / x, modular) == 1) {
                flag = false;
                break;
            }

        if (flag) return j;
    }

    exit(-1);
}

valueType roundToInt(double a, valueType modular) {
    return a >= 0 ? (valueType) a : (valueType) a - modular;
}

double convertToDouble(valueType a, valueType modular) {
    valueType mh = modular >> 1u;
    return (a <= mh) ? (double) a : (double) a - mh;
}

void mulModAssign(valueType &res, valueType a, valueType b, valueType modular)  {
    res = multiply(a, b, modular);
}

void addModAssign(valueType &res, valueType a, valueType b, valueType modular) {
    a %= modular;
    b %= modular;
    valueType da = modular - a;
    if (b < da) res = a + b;
    else res = b - da;
}

void subModAssign(valueType &res, valueType a, valueType b, valueType modular) {
    a %= modular;
    b %= modular;
    if (b <= a) res = (a - b);
    else addModAssign(res, a, modular - b, modular);
    __int128 ans = a;
    ans = (ans + modular - b) % modular;
}

void negModAssign(valueType &res, valueType a, valueType modular) {
    a %= modular;
    res = modular - a;
}

void mulAddModAssign(valueType &res, valueType a, valueType b, valueType c, valueType modular) {
    res = addMod(mulMod(a, b, modular), c, modular);
}

void modAssign(valueType &res, valueType a, valueType modular) {
    res = a % modular;
}

valueType mulMod(valueType a, valueType b, valueType modular) {
    valueType res;
    mulModAssign(res, a, b, modular);
    return res;
}

valueType addMod(valueType a, valueType b, valueType modular) {
    valueType res;
    addModAssign(res, a, b, modular);
    return res;
}

valueType subMod(valueType a, valueType b, valueType modular) {
    valueType res;
    subModAssign(res, a, b, modular);
    return res;
}

valueType mod(valueType a, valueType modular) {
    valueType res;
    modAssign(res, a, modular);
    return res;
}

valueType mulAddMod(valueType a, valueType b, valueType c, valueType modular) {
    valueType res;
    mulAddModAssign(res, a, b, c, modular);
    return res;
}

valueType negMod(valueType a, valueType modular) {
    a %= modular;
    return modular - a;
}

valueType power(valueType x, valueType y, valueType p)  {
    valueType res = 1;

    for(; y; y >>= 1u) {
        if (y & 1u) mulModAssign(res, res, x, p);
        mulModAssign(x, x, x, p);
    }

    return res;
}

valueType multiply(valueType x, valueType y, valueType p) {
    valueType res = 0;

    for (; y; y >>= 1u) {
        if (y & 1u) addModAssign(res, res, x, p);
        addModAssign(x, x, x, p);
    }

    return res;
}

void CTTransformToStd(vector<complex<double>> &a, const vector<sizeType> &rev, const vector<complex<double>> &psi,
                      const vector<sizeType> &rotGroup) {
    sizeType Nh = a.size();
    arrayBitReverse(a, rev);

    // Assume M / N = 2, M / Nh = 4;
    for (int len = 2, sc = (int)round(log2(Nh)) - 1; len <= Nh; len <<= 1, --sc) {
        int lenh = len >> 1;
        int msk = (len << 2) - 1;
        for (int i = 0; i < Nh; i += len) {
            auto j2 = i + lenh;
            for (int j = i; j < j2; ++j) {
                auto u = a[j];
                auto v = a[j + lenh] * psi[(rotGroup[j - i] & msk) << sc];
                a[j] = u + v;
                a[j + lenh] = u - v;
            }
        }
    }
}

void CTTransformToRev(vector<valueType> &a, valueType modular, const vector<valueType> &psi,
                      const vector<sizeType> &rotGroupRev) {
    sizeType Nh = a.size();
    sizeType msk = (Nh << 2) - 1; // mod M = & msk

    // Assume M / N = 2, M / Nh = 4;
    for (int len = Nh; len > 1; len >>= 1) {
        int lenh = len >> 1;
        for (int i = 0; i < Nh; i = i + len) {
            auto s = psi[rotGroupRev[i] * lenh & msk];
            int j2 = i + lenh;
            for (int j = i; j < j2; ++j) {
                auto u = a[j];
                auto v = mulMod(a[j + lenh], s, modular);
                a[j] = addMod(u, v, modular);
                a[j + lenh] = subMod(u, v, modular);
            }
        }
    }
}

void GSTransformFromStd(vector<complex<double>> &a, const vector<sizeType> &rev, const vector<complex<double>> &ipsi,
                        const vector<sizeType> &rotGroup) {
    sizeType Nh = a.size();
    for (int len = Nh, sc = 0; len; len >>= 1, ++sc) {
        int lenh = len >> 1;
        int msk = (len << 2) - 1;

        // Assume M / N = 2, M / Nh = 4;
        for (int i = 0; i < Nh; i += len) {
            int j2 = i + lenh;
            for (int j = i; j < j2; ++j) {
                auto u = a[j] + a[j + lenh];
                auto v = a[j] - a[j + lenh];
                v *= ipsi[(rotGroup[j - i] & msk) << sc];
                a[j] = u;
                a[j + lenh] = v;
            }
        }
    }

    arrayBitReverse(a, rev);
}

void GSTransformFromRev(vector<valueType> &a, valueType modular, const vector<valueType> &ipsi,
                        const vector<sizeType> &rotGroupRev) {
    int Nh = a.size();
    sizeType msk = (Nh << 2) - 1;     // mod M = & msk

    // Assume M / N = 2, M / Nh = 4;
    for (int len = 2; len <= Nh; len <<= 1) {
        int lenh = len >> 1;
        for (int i = 0; i < Nh; i = i + len) {
            auto s = ipsi[rotGroupRev[i] * lenh & msk];
            auto j2 = i + lenh;
            for (int j = i; j < j2; ++j) {
                auto u = addMod(a[j], a[j + lenh], modular);
                auto v = subMod(a[j], a[j + lenh], modular);
                mulModAssign(v, v, s, modular);
                a[j] = u;
                a[j + lenh] = v;
            }
        }
    }
}

void fft(vector<double> &a, const vector<sizeType> &rev, const vector<complex<double>> &psi,
         const vector<sizeType> &rotGroup) {
    static vector<complex<double>> na;

    sizeType Nh = a.size();
    na.resize(Nh);
    na[0] = a[0];
    auto s = psi[(Nh << 2) - Nh];
    for (int n = 1; n < Nh; ++n)
        na[n] = a[n] + a[Nh - n] * s;


    CTTransformToStd(na, rev, psi, rotGroup);

    for (int n = 0; n < Nh; ++n)
        a[n] = na[n].real();
}

void ifft(vector<double> &a, double scale, const vector<sizeType> &rev, const vector<complex<double>> &ipsi,
          const vector<sizeType> &rotGroup) {
    static vector<complex<double>> na;

    sizeType Nh = a.size();
    na.resize(Nh);
    for (int n = 0; n < Nh; ++n)
        na[n] = a[n];

    GSTransformFromStd(na, rev, ipsi, rotGroup);

    double sc = 0.5 * scale / Nh;
    auto s = ipsi[(Nh << 2) - Nh];
    a[0] = na[0].real() * scale / Nh;
    for (int n = 1; n < Nh; ++n)
        a[n] = (na[n] + na[Nh - n] * s).real() * sc;
}

void ntt(vector<valueType> &a, valueType modular, const vector<sizeType> &rotGroupRev, const vector<valueType> &psi) {
    sizeType Nh = a.size();
    auto s = psi[(Nh << 2) - Nh];
    for (int n = 1; (n << 1) <= Nh; ++n) {
        valueType tmp1 = addMod(a[n], mulMod(a[Nh - n], s, modular), modular);
        valueType tmp2 = addMod(a[Nh - n], mulMod(a[n], s, modular), modular);
        a[n] = tmp1;
        a[Nh - n] = tmp2;
    }
    CTTransformToRev(a, modular, psi, rotGroupRev);
}

void intt(vector<valueType> &a, valueType modular, const vector<sizeType> &rotGroupRev, const vector<valueType> &ipsi) {
    GSTransformFromRev(a, modular, ipsi, rotGroupRev);

    sizeType Nh = a.size();
    auto s = ipsi[(Nh << 2) - Nh];
    valueType inv2N = power(addMod(Nh, Nh, modular), modular - 2, modular);
    mulModAssign(a[0], a[0], addMod(inv2N, inv2N, modular), modular);
    for (int n = 1; (n << 1) <= Nh; ++n) {
        valueType tmp1 = addMod(a[n], mulMod(a[Nh - n], s, modular), modular);
        valueType tmp2 = addMod(a[Nh - n], mulMod(a[n], s, modular), modular);

        mulModAssign(a[n], inv2N, tmp1, modular);
        mulModAssign(a[Nh - n], inv2N, tmp2, modular);
    }
}
