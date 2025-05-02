#include "cholesky.h"
#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"
#include <print>

void solveFilip(VecDoub x, VecDoub y, int parameters) {
  // FilipData
  // util::print(x, "xFilip");
  // util::print(y, "yFilip");

  int size = x.size();
  MatDoub A(size, parameters);
  VecDoub b(size);
  // NO sigma! (sigma=1)
  double sigma = 1.;
  for (int i = 0; i < size; i++) {
    // a1 + a2 x + a3 x^2
    for (int j = 0; j < parameters; j++) {
      A[i][j] = pow(x[i], j) / sigma;
    }
    b[i] = y[i] / sigma;
  }
  MatDoub A_T = util::Transpose(A);
  MatDoub C = A_T * A;
  VecDoub c = A_T * b;
  util::print(b);

  // Using LU decomposition
  try {
    auto ludcmp_solver = LUdcmp(C);
    VecDoub a_ludcmp(parameters);
    ludcmp_solver.solve(c, a_ludcmp);
    util::print(a_ludcmp, "LU decomposition");
  } catch (...) {
    std::println("Error in LU decomposition");
  }

  // Using Cholesky
  try {
    VecDoub a_cholesky(parameters);
    auto cholesky_solver = Cholesky(C);
    cholesky_solver.solve(c, a_cholesky);
    util::print(a_cholesky, "Cholesky");

  } catch (...) {
    std::println("Error in Cholesky");
  }
}

void solvePont(VecDoub x, VecDoub y, int parameters) {
  // PontiusData
  // util::print(x, "xPont");
  // util::print(y, "yPont");

  int size = x.size();
  MatDoub A(size, parameters);
  VecDoub b(size);
  // NO sigma! (sigma=1)
  double sigma = 1.;
  for (int i = 0; i < size; i++) {
    // a1 + a2 x + a3 x^2
    for (int j = 0; j < parameters; j++) {
      A[i][j] = pow(x[i], j) / sigma;
    }
    b[i] = y[i] / sigma;
  }
  MatDoub A_T = util::Transpose(A);
  MatDoub C = A_T * A;
  VecDoub c = A_T * b;
  util::print(b);

  // Using LU decomposition
  try {
    auto ludcmp_solver = LUdcmp(C);
    VecDoub a_ludcmp(parameters);
    ludcmp_solver.solve(c, a_ludcmp);
    util::print(a_ludcmp, "LU decomposition");
  } catch (...) {
    std::println("Error in LU decomposition");
  }

  // Using Cholesky
  try {
    VecDoub a_cholesky(parameters);
    auto cholesky_solver = Cholesky(C);
    cholesky_solver.solve(c, a_cholesky);
    util::print(a_cholesky, "Cholesky");

  } catch (...) {
    std::println("Error in Cholesky");
  }
}

int main() {
  // Load data from FilipData.dat
  VecDoub xFilip(82);
  VecDoub yFilip(82);
  std::ifstream Filip("../../../2/FilipData.dat");
  if (!Filip) {
    std::println(stderr, "Error opening FilipData.dat");
    return 1;
  }
  for (int i = 0; i < 82; i++) {
    Filip >> yFilip[i];
    Filip >> xFilip[i];
  }

  // Load data from PontiusData.dat
  VecDoub xPont(40);
  VecDoub yPont(40);
  std::ifstream Pont("../../../2/PontiusData.dat");
  if (!Pont) {
    std::println(stderr, "Error opening PontiusData.dat");
    return 1;
  }
  for (int i = 0; i < 40; i++) {
    Pont >> yPont[i];
    Pont >> xPont[i];
  }

  solvePont(xPont, yPont, 3);
  solveFilip(xFilip, yFilip, 11);

  return 0;
}
