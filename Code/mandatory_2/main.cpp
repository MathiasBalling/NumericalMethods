#include "nr3.h"
#include "roots_multidim.h"
#include "utilities.h"
#include <cassert>
#include <cmath>
#include <print>

// Material constants:
const double v = 120.0;                 // kg
const double k = 2.5;                   // m
const double w = 4.0;                   // kg/m
const double alpha = 2.0 * pow(10, -7); // kg^-1

// Other constants:
const double d = 30.0;
double n = 5.0;

VecDoub vecfunc(VecDoub_I q) {
  // q=[L0, L, p, x, theta, phi/varphi, a, H]
  assert(q.size() == 8);

  // Zero vectors
  VecDoub f(8);
  f[0] = (q[6] * (cosh(q[3] / q[6]) - 1.0)) - q[2];
  f[1] = (2.0 * q[6] * sinh(q[3] / q[6])) - q[1];
  f[2] = (2.0 * q[3] + 2.0 * k * cos(q[4])) - d;
  f[3] = (q[2] + k * sin(q[4])) - n;
  f[4] = (sinh(q[3] / q[6])) - tan(q[5]);
  f[5] = ((1.0 + (v / (w * q[0]))) * tan(q[5])) - tan(q[4]);
  f[6] = (q[0] * (1.0 + alpha * q[7])) - q[1];
  f[7] = ((w * q[0]) / (2.0 * sin(q[5]))) - q[7];

  return f;
}

void print_roots_table(std::vector<VecDoub> x_k) {
  std::vector<VecDoub> d_k;
  d_k.push_back(VecDoub(2));
  for (size_t i = 1; i < x_k.size(); i++) {
    d_k.push_back(x_k[i] - x_k[i - 1]);
  }

  // Calculate the convergence constant C
  std::vector<double> C;
  C.push_back(NAN);
  C.push_back(NAN);
  for (size_t i = 2; i < d_k.size(); i++) {
    // Take norm of d_k
    C.push_back(util::norm(d_k[i]) / (pow(util::norm(d_k[i - 1]), 2)));
  }

  // Calculate the error
  std::vector<double> e_k;
  for (size_t i = 0; i < d_k.size(); i++) {
    e_k.push_back(C.at(i) * pow(util::norm(d_k[i]), 2));
  }

  // Print the table
  // Table header
  std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|", "k", "x_k", "d_k", "C",
               "e");
  // Table header separator
  std::println("|{:-^10}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|", "", "", "", "", "");

  // Table body
  for (size_t i = 0; i < x_k.size(); i++) {
    std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|", i + 1,
                 x_k.at(i)[0], d_k.at(i)[0], C.at(i), e_k.at(i));
    for (int j = 1; j < x_k.at(i).size(); j++) {
      std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21}|{:^21}|", i + 1,
                   x_k.at(i)[j], d_k.at(i)[j], "", "");
    }
  }
}

template <class T>
void newton_multi(VecDoub_IO &x, Bool &check, T &vecfunc, int max_its = 200,
                  bool print_table = false) {
  std::vector<VecDoub> x_k; // For table
  const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
  const Doub TOLX = numeric_limits<Doub>::epsilon();
  Int i, j, its, n = x.size();
  Doub den, f, fold, stpmax, sum, temp, test;
  VecDoub g(n), p(n), xold(n);
  MatDoub fjac(n, n);
  NRfmin<T> fmin(vecfunc);
  NRfdjac<T> fdjac(vecfunc);
  VecDoub &fvec = fmin.fvec;
  f = fmin(x);
  test = 0.0;
  for (i = 0; i < n; i++)
    if (abs(fvec[i]) > test)
      test = abs(fvec[i]);
  if (test < 0.01 * TOLF) {
    check = false;
    return;
  }
  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += SQR(x[i]);
  stpmax = STPMX * MAX(sqrt(sum), Doub(n));
  for (its = 0; its < max_its; its++) {
    fjac = fdjac(x, fvec);
    for (i = 0; i < n; i++) {
      sum = 0.0;
      for (j = 0; j < n; j++)
        sum += fjac[j][i] * fvec[j];
      g[i] = sum;
    }
    for (i = 0; i < n; i++)
      xold[i] = x[i];
    fold = f;
    for (i = 0; i < n; i++)
      p[i] = -fvec[i];
    LUdcmp alu(fjac);
    alu.solve(p, p);
    lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (abs(fvec[i]) > test)
        test = abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      x_k.push_back(x);
      if (print_table) {
        print_roots_table(x_k);
      } else {
        println("Done in {} iterations.", its);
      }
      return;
    }
    if (check) {
      test = 0.0;
      den = MAX(f, 0.5 * n);
      for (i = 0; i < n; i++) {
        temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
        if (temp > test)
          test = temp;
      }
      check = (test < TOLMIN);
      x_k.push_back(x);
      if (print_table) {
        print_roots_table(x_k);
      } else {
        println("Done in {} iterations.", its);
      }
      return;
    }
    test = 0.0;
    for (i = 0; i < n; i++) {
      temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX) {
      x_k.push_back(x);
      if (print_table) {
        print_roots_table(x_k);
      } else {
        println("Done in {} iterations.", its);
      }
      return;
    }
    x_k.push_back(x);
  }
  throw("MAXITS exceeded in newt");
}

int main() {

  {
    println("\n{:/^120}", " n=5.0 ");
    n = 5.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {29.0, 29.1, 5.0, 15.0, 1.0, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }

  {
    println("\n{:/^120}", " n=2.0 ");
    n = 2.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {27.0, 27.1, 2.0, 14.0, 1.0, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }

  {
    println("\n{:/^120}", " n=1.0 ");
    n = 1.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {26.0, 26.1, 1.0, 13.0, 0.5, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }

  {
    println("\n{:/^120}", " n=0.5 ");
    n = 0.5;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {26.0, 26.1, 0.5, 13.0, 0.2, 0.1, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }

  {
    println("\n{:/^120}", " n=0.2 ");
    n = 0.2;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {25.5, 25.6, 0.2, 12.0, 0.1, 0.05, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }

  {
    println("\n{:/^120}", " n=0.1 ");
    n = 0.1;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {25.5, 25.6, 0.1, 12.0, 0.1, 0.05, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    println("Local minimum: {}", check);
  }
}
