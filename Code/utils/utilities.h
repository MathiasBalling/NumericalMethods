/*
 * utilities.h
 *
 *  Created on: Feb 2, 2015
 *      Author: jpb
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "nr3.h"
#include <iostream>
#include <string>

namespace util {
void print(MatDoub mat, string symbol = "") {
  if (symbol.compare(""))
    std::cout << symbol << "	Matrix " << mat.nrows() << "x" << mat.ncols()
              << ":" << endl;

  for (int m = 0; m < mat.nrows(); m++) {
    for (int n = 0; n < mat.ncols(); n++) {
      std::cout << setw(15) << mat[m][n] << "\t";
    }
    std::cout << endl;
  }
  std::cout << endl;
}

void printDiag(MatDoub mat, string symbol = "") {
  if (symbol.compare(""))
    std::cout << symbol << "	MatrixDiag " << mat.nrows() << "x"
              << mat.ncols() << ":" << endl;
  double nmax = mat.nrows() < mat.nrows() ? mat.nrows() : mat.nrows();
  for (int n = 0; n < nmax; n++) {
    std::cout << setw(15) << mat[n][n] << "\t";
  }
  std::cout << endl;
}

MatDoub diag(VecDoub &V) {
  double m = V.size();
  MatDoub M(m, m);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      M[i][j] = 0;

  for (int i = 0; i < m; i++)
    M[i][i] = V[i];

  return M;
}

void print(VecDoub vec, string symbol = "") {
  if (symbol.compare(""))
    std::cout << symbol << "	Vector " << vec.size() << "D:" << endl;

  for (int m = 0; m < vec.size(); m++) {
    std::cout << setw(15) << vec[m];
  }
  std::cout << endl;
}

MatDoub Transpose(const MatDoub &Mat) {
  MatDoub res(Mat.ncols(), Mat.nrows());
  for (int n = 0; n < res.nrows(); n++) {
    for (int m = 0; m < res.ncols(); m++) {
      res[n][m] = Mat[m][n];
    }
  }
  return res;
}

MatDoub T(const MatDoub &Mat) { return Transpose(Mat); }

double norm(VecDoub &v) {
  double res = 0;
  for (int i = 0; i < v.size(); i++) {
    res += v[i] * v[i];
  }
  return sqrt(res);
}
} // namespace util

MatDoub operator*(const MatDoub &A1, const MatDoub &A2) {
  if (A1.ncols() != A2.nrows()) {
    cerr << "in prod: the number of rows in A1 is not equal to the number of "
            "cols in A2"
         << endl;
  }

  MatDoub res(A1.nrows(), A2.ncols());
  for (int n = 0; n < A1.nrows(); n++) {
    for (int m = 0; m < A2.ncols(); m++) {
      double temp = 0;
      for (int i = 0; i < A1.ncols(); i++) {
        temp += A1[n][i] * A2[i][m];
      }
      res[n][m] = temp;
    }
  }
  return res;
}

VecDoub operator*(const MatDoub &A, const VecDoub &b) {
  if (A.ncols() != b.size()) {
    cerr << "in prod: the number of rows in A is not equal to the size of "
            "vector b"
         << endl;
  }

  VecDoub res(A.nrows());
  for (int n = 0; n < A.nrows(); n++) {
    double temp = 0;
    for (int m = 0; m < A.ncols(); m++) {
      temp += A[n][m] * b[m];
    }
    res[n] = temp;
  }
  return res;
}

VecDoub operator-(const VecDoub &v1, const VecDoub &v2) {
  if (v1.size() != v2.size()) {
    std::cout << "in minus: the size of  cclv1 is not equal to the size of "
                 "vector b"
              << std::endl;
  }
  VecDoub res(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    res[i] = v1[i] - v2[i];
  }
  return res;
}

VecDoub operator+(const VecDoub &a, const VecDoub &b) {
  if (a.size() != b.size()) {
    std::cout << "in prod: the number of rows in A is not equal to the size of "
                 "vector b"
              << std::endl;
  }
  VecDoub res(a.size());
  for (int i = 0; i < a.size(); i++) {
    res[i] = a[i] + b[i];
  }
  return res;
}

VecDoub operator/(const VecDoub &v, double s) {
  VecDoub res(v.size());
  for (int i = 0; i < v.size(); i++) {
    res[i] = v[i] / s;
  }
  return res;
}

#endif /* UTILITIES_H_ */
