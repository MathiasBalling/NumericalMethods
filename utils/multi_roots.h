#ifndef MULTI_ROOTS_H
#define MULTI_ROOTS_H
#include "nr3.h"
#include "roots_multidim.h"
#include "utilities.h"

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
                  Doub tolerance = 1.0e-8, bool print_table = false) {
  std::vector<VecDoub> x_k; // For table
  const Doub TOLF = tolerance, TOLMIN = 1.0e-12, STPMX = 100.0;
  const Doub TOLX = std::numeric_limits<Doub>::epsilon();
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
      }
      // else {
      //   println("Done in {} iterations.", its);
      // }
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
      }
      // else {
      //   println("Done in {} iterations.", its);
      // }
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
      }
      // else {
      //   println("Done in {} iterations.", its);
      // }
      return;
    }
    x_k.push_back(x);
  }
  throw("MAXITS exceeded in newt");
}
#endif // !MULTI_ROOTS_H
