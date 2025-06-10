#ifndef MULTI_ROOTS_H
#define MULTI_ROOTS_H
#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"

static bool using_backtracking = false;
static double lambda_backtracking = 0.0;

void print_roots_table(std::vector<std::pair<VecDoub, double>> x_k_and_lambda,
                       int max_its = 200, double accuracy = 1.0e-8) {
  std::vector<VecDoub> x_k;
  std::vector<double> lambda_k;
  for (const auto &pair : x_k_and_lambda) {
    x_k.push_back(pair.first);
    lambda_k.push_back(pair.second);
  }

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
  std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", "k", "x_k", "d_k",
               "C", "e", "lambda");
  // Table header separator
  std::println("|{:-^10}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|", "", "", "",
               "", "", "");

  // Table body
  for (size_t i = 0; i < x_k.size(); i++) {
    if (i == 0) {
      if (lambda_k[i] != 0.0) {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", i + 1, "",
                     "", "", "", lambda_k[i]);
      } else {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", i + 1, "",
                     "", "", "", "No bt");
      }
      for (int j = 0; j < x_k.at(i).size(); j++) {
        std::println("|{:^10}|{:^21.12}|{:^21}|{:^21}|{:^21}|{:^21.12}|",
                     std::format("{}({})", (i + 1), j), x_k.at(i)[j], "", "",
                     "", "");
      }
    } else if (i == 1) {
      if (lambda_k[i] != 0.0) {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", i + 1, "",
                     "", "", "", lambda_k[i]);
      } else {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", i + 1, "",
                     "", "", "", "No bt");
      }
      for (int j = 0; j < x_k.at(i).size(); j++) {
        std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21}|{:^21}|{:^21}|",
                     std::format("{}({})", (i + 1), j), x_k.at(i)[j],
                     d_k.at(i)[j], "", "", "");
      }
    } else {
      if (lambda_k[i] != 0.0) {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21.12}|{:^21.12}|{:^21.12}|",
                     i + 1, "", "", C.at(i), e_k.at(i), lambda_k[i]);
      } else {
        std::println("|{:^10}|{:^21}|{:^21}|{:^21.12}|{:^21.12}|{:^21}|", i + 1,
                     "", "", C.at(i), e_k.at(i), "No bt");
      }
      for (int j = 0; j < x_k.at(i).size(); j++) {
        std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21}|{:^21}|{:^21}|",
                     std::format("{}({})", (i + 1), j), x_k.at(i)[j],
                     d_k.at(i)[j], "", "", "");
      }
      if (i + 2 > max_its || e_k.at(i) < accuracy) {
        return;
      }
    }
  }
}

template <class T>
void lnsrch_custom(VecDoub_I &xold, const Doub fold, VecDoub_I &g,
                   VecDoub_IO &p, VecDoub_O &x, Doub &f, const Doub stpmax,
                   Bool &check, T &func) {
  using_backtracking = false; // Reset using_backtracking

  const Doub ALF = 1.0e-4, TOLX = std::numeric_limits<Doub>::epsilon();
  Doub a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
  Doub rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
  Int i, n = xold.size();
  check = false;
  for (i = 0; i < n; i++)
    sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 0; i < n; i++)
      p[i] *= stpmax / sum;
  for (i = 0; i < n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0)
    throw("Roundoff problem in lnsrch.");
  test = 0.0;
  for (i = 0; i < n; i++) {
    temp = abs(p[i]) / MAX(abs(xold[i]), 1.0);
    if (temp > test)
      test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;) {
    for (i = 0; i < n; i++)
      x[i] = xold[i] + alam * p[i];
    f = func(x);
    if (alam < alamin) {
      for (i = 0; i < n; i++)
        x[i] = xold[i];
      check = true;
      return;
    } else if (f <= fold + ALF * alam * slope) {
      if (using_backtracking) {
        lambda_backtracking = alam;
      }
      return;
    } else {
      // Using backtracking
      using_backtracking = true;
      if (alam == 1.0)
        tmplam = -slope / (2.0 * (f - fold - slope));
      else {
        rhs1 = f - fold - alam * slope;
        rhs2 = f2 - fold - alam2 * slope;
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) /
            (alam - alam2);
        if (a == 0.0)
          tmplam = -slope / (2.0 * b);
        else {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0.0)
            tmplam = 0.5 * alam;
          else if (b <= 0.0)
            tmplam = (-b + sqrt(disc)) / (3.0 * a);
          else
            tmplam = -slope / (b + sqrt(disc));
        }
        if (tmplam > 0.5 * alam)
          tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = f;
    alam = MAX(tmplam, 0.1 * alam);
  }
}

template <class T> struct NRfdjac {
  const Doub EPS;
  T &func;
  NRfdjac(T &funcc) : EPS(1.0e-8), func(funcc) {}
  MatDoub operator()(VecDoub_I &x, VecDoub_I &fvec) {
    Int n = x.size();
    MatDoub df(n, n);
    VecDoub xh = x;
    for (Int j = 0; j < n; j++) {
      Doub temp = xh[j];
      Doub h = EPS * abs(temp);
      if (h == 0.0)
        h = EPS;
      xh[j] = temp + h;
      h = xh[j] - temp;
      VecDoub f = func(xh);
      xh[j] = temp;
      for (Int i = 0; i < n; i++)
        df[i][j] = (f[i] - fvec[i]) / h;
    }
    return df;
  }
};

template <class T> struct NRfmin {
  VecDoub fvec;
  T &func;
  Int n;
  NRfmin(T &funcc) : func(funcc) {}
  Doub operator()(VecDoub_I &x) {
    n = x.size();
    Doub sum = 0;
    fvec = func(x);
    for (Int i = 0; i < n; i++)
      sum += SQR(fvec[i]);
    return 0.5 * sum;
  }
};

// `newt` function with data table
template <class T>
void newton_multi_table(VecDoub_IO &x, Bool &check, T &vecfunc,
                        int max_its = 200, Doub accuracy = 1.0e-8) {
  std::vector<std::pair<VecDoub, double>>
      x_k; // For table x and backtracking lambda
  const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
  const Doub TOLX = std::numeric_limits<Doub>::epsilon();
  const Doub MAXITS = 200;
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
  for (its = 0; its < MAXITS; its++) {
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
    lnsrch_custom(xold, fold, g, p, x, f, stpmax, check, fmin);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (abs(fvec[i]) > test)
        test = abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      if (using_backtracking) {
        x_k.push_back({x, lambda_backtracking}); // For table
      } else {
        x_k.push_back({x, 0.0}); // For table
      }
      print_roots_table(x_k, max_its, accuracy); // For table
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
      if (using_backtracking) {
        x_k.push_back({x, lambda_backtracking}); // For table
      } else {
        x_k.push_back({x, 0.0}); // For table
      }
      print_roots_table(x_k, max_its, accuracy); // For table
      return;
    }
    test = 0.0;
    for (i = 0; i < n; i++) {
      temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX) {
      if (using_backtracking) {
        x_k.push_back({x, lambda_backtracking}); // For table
      } else {
        x_k.push_back({x, 0.0}); // For table
      }
      print_roots_table(x_k, max_its, accuracy); // For table
      return;
    }
    if (using_backtracking) {
      x_k.push_back({x, lambda_backtracking}); // For table
    } else {
      x_k.push_back({x, 0.0}); // For table
    }
  }
  throw("MAXITS exceeded in newt");
}
#endif // !MULTI_ROOTS_H
