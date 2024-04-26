//
// Created by 69029 on 3/22/2020.
//

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <cstdint>
#include <vector>
#include <emp-tool/emp-tool.h>
#include <random>
#include "CommonHeader.h"

using namespace std;

bool millerTest(valueType d, valueType n);

bool isPrime(valueType n);

vector<valueType> factor(valueType x);

valueType findPrimitives(valueType modular);

valueType roundToInt(double a, valueType modular);

double convertToDouble(valueType a, valueType modular);
template<class T>
void getRandomInRange(T &res, T left, T right) {
    static default_random_engine generator((unsigned) time (NULL));
    uniform_int_distribution<T> distribution(left, right);
    res = distribution(generator);
}

template<class T>
T randomTwo(T a, T b) {
    short x;
    getRandomInRange(x, (short) 0, (short) 1);
    return x ? a : b;
}

template <class S, class T>
void getRandomInRange(S &res, T left, T right) {
    T tmp;
    getRandomInRange(tmp, left, right);
    res = tmp;
}

/**
 * FFT and related functions for encoder
 */

template <class T>
void arrayBitReverse(vector<T> &data, const vector<sizeType> &rev) {
    assert(data.size() == rev.size());
    for (int n = 0; n < data.size(); ++n)
        if (n < rev[n]) swap(data[n], data[rev[n]]);
}

void CTTransformToStd(vector<complex<double>> &a, const vector<sizeType> &rev,
                      const vector<complex<double>> &psi, const vector<sizeType> &rotGroup);

void CTTransformToRev(vector<valueType> &a, valueType modular, const vector<valueType> &psi,
                      const vector<sizeType> &rotGroupRev);

void GSTransformFromStd(vector<complex<double>> &a, const vector<sizeType> &rev,
                        const vector<complex<double>> &ipsi, const vector<sizeType> &rotGroup);

void GSTransformFromRev(vector<valueType> &a, valueType modular, const vector<valueType> &ipsi,
                        const vector<sizeType> &rotGroupRev);

void fft(vector<double> &a, const vector<sizeType> &rev, const vector<complex<double>> &ipsi,
         const vector<sizeType> &rotGroup);

void ifft(vector<double> &a, double scale, const vector<sizeType> &rev, const vector<complex<double>> &ipsi,
          const vector<sizeType> &rotGroup);

void ntt(vector<valueType> &a, valueType modular, const vector<sizeType> &rotGroupRev, const vector<valueType> &psi);

void intt(vector<valueType> &a, valueType modular, const vector<sizeType> &rotGroupRev, const vector<valueType> &ipsi);

void mulModAssign(valueType &res, valueType a, valueType b, valueType modular);

void addModAssign(valueType &res, valueType a, valueType b, valueType modular);

void subModAssign(valueType &res, valueType a, valueType b, valueType modular);

void negModAssign(valueType &res, valueType a, valueType modular);

void mulAddModAssign(valueType &res, valueType a, valueType b, valueType c, valueType modular);

void modAssign(valueType &res, valueType a, valueType modular);

valueType mulMod(valueType a, valueType b, valueType modular);

valueType addMod(valueType a, valueType b, valueType modular);

valueType subMod(valueType a, valueType b, valueType modular);

valueType negMod(valueType a, valueType modular);

valueType mulAddMod(valueType a, valueType b, valueType c, valueType modular);

valueType mod(valueType a, valueType modular);

valueType multiply(valueType x, valueType y, valueType p);

valueType power(valueType x, valueType y, valueType p);

// TODO: to be removed.

template <class T>
void printVec(ostream &out, const string &str, const vector<T> &vec) {
    cout << str << ':';

    if (vec.empty()) return;
    const sizeType N = min((size_t) 16, vec.size());
    for_each(vec.begin(), vec.begin() + N, [](T a) {cout << ' ' << a;});
    cout << (16 < vec.size() ? " ...\n" : "\n");
}

template <class T>
void printMat(ostream &out, const string &str, const vector<vector<T>> &mat) {
    out << "  " << str << ": ";

    if (!mat.empty() && !mat[0].empty()) {
        const sizeType N = min((size_t) 16, mat.size());
        const sizeType M = min((size_t) 16, mat[0].size());
        for (auto v = mat.begin(); v != mat.begin() + N; ++v) {
            for (auto x = v->begin(); x != v->begin() + M; ++x) {
                out.width(19);
                out << (*x) << ", ";
            }
            if (v->size() > 8) out << " ...";
            if (v != mat.begin() + N - 1) {
                out << '\n';
                out.width(4 + str.length());
                out << ' ';
            }
        }
        if (mat.size() > 8) out << "\n....";
    }
    out << '\n';
}

#endif //SRC_UTILS_H
