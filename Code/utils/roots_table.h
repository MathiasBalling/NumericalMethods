#ifndef ROOTS_TABLE_H
#define ROOTS_TABLE_H
#include "nr3.h"
#include <cmath>
#include <print>

namespace roots_table {
enum class RootType { Bisection, Ridder, Secant, FalsePositive, Newton };
void print_table(const std::vector<double> &x_k, RootType method) {
  std::vector<double> d_k;
  d_k.push_back(NAN);
  for (size_t i = 1; i < x_k.size(); i++) {
    d_k.push_back(x_k[i] - x_k[i - 1]);
  }

  // Calculate the convergence constant C
  std::vector<double> C;
  for (size_t i = 0; i < d_k.size(); i++) {
    switch (method) {
    case RootType::Bisection:
      C.push_back(0.5);
      break;
    case RootType::FalsePositive:
      C.push_back(abs(d_k[i]) / abs(d_k[i - 1]));
      break;
    case RootType::Secant:
      C.push_back(abs(d_k[i]) / (pow(abs(d_k[i - 1]), 1.62)));
      break;
    case RootType::Ridder:
      C.push_back(abs(d_k[i]) / (pow(abs(d_k[i - 1]), 3)));
      break;
    case RootType::Newton:
      C.push_back(abs(d_k[i]) / (pow(abs(d_k[i - 1]), 2)));
      break;
    default:
      throw("Unknown method");
      break;
    }
  }

  // Calculate the error
  std::vector<double> e_k;
  for (size_t i = 0; i < d_k.size(); i++) {
    switch (method) {
    case RootType::Bisection:
      e_k.push_back(abs(d_k[i]));
      break;
    case RootType::FalsePositive:
      e_k.push_back((-C.at(i) / (1 - C.at(i))) * d_k[i]);
      break;
    case RootType::Secant:
      e_k.push_back(C.at(i) * pow(abs(d_k[i]), 1.62));
      break;
    case RootType::Ridder:
      e_k.push_back(C.at(i) * pow(abs(d_k[i]), 3));
      break;
    case RootType::Newton:
      e_k.push_back(C.at(i) * pow(abs(d_k[i]), 2));
      break;
    default:
      throw("Unknown method");
      break;
    }
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
                 x_k.at(i), d_k.at(i), C.at(i), e_k.at(i));
  }
}

template <class T>
Doub root_bisection(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  std::vector<double> x_k; // For table
  const Int MAX_ITER = 50;
  Doub dx, xmid, rtb;
  Doub f = func(x1);
  Doub fmid = func(x2);
  if (f * fmid >= 0.0)
    throw("Root must be bracketed for bisection in rtbis");

  rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (Int j = 0; j < MAX_ITER; j++) {
    fmid = func(xmid = rtb + (dx *= 0.5));
    x_k.push_back(xmid); // For table
    if (fmid <= 0.0)
      rtb = xmid;
    if (abs(dx) < xacc || fmid == 0.0) {
      print_table(x_k, RootType::Bisection); // Print table
      return rtb;
    }
  }
  throw("Too many bisections in rtbis");
}

template <class T>
Doub root_ridder(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  std::vector<double> x_k; // For table
  const Int MAXIT = 60;
  Doub fl = func(x1);
  Doub fh = func(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    Doub xl = x1;
    Doub xh = x2;
    Doub ans = -9.99e99;
    for (Int j = 0; j < MAXIT; j++) {
      Doub xm = 0.5 * (xl + xh);
      Doub fm = func(xm);
      Doub s = sqrt(fm * fm - fl * fh);
      if (s == 0.0) {
        print_table(x_k, RootType::Ridder); // Print table
        return ans;
      }
      Doub xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (abs(xnew - ans) <= xacc) {
        print_table(x_k, RootType::Ridder); // Print table
        return ans;
      }
      ans = xnew;
      x_k.push_back(ans); // For table
      Doub fnew = func(ans);
      if (fnew == 0.0) {
        print_table(x_k, RootType::Ridder); // Print table
        return ans;
      }
      if (SIGN(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (SIGN(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (SIGN(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else
        throw("never get here.");
      if (abs(xh - xl) <= xacc) {
        print_table(x_k, RootType::Ridder); // Print table
        return ans;
      }
    }
    throw("zriddr exceed maximum iterations");
  } else {
    if (fl == 0.0)
      return x1;
    if (fh == 0.0)
      return x2;
    throw("root must be bracketed in zriddr.");
  }
}

template <class T>
Doub root_secant(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  std::vector<double> x_k;
  const int MAX_ITERTER = 30;
  Doub xl, rts;
  Doub fl = func(x1);
  Doub f = func(x2);
  if (abs(fl) < abs(f)) {
    rts = x1;
    xl = x2;
    SWAP(fl, f);
  } else {
    xl = x1;
    rts = x2;
  }
  x_k.push_back(rts); // For table
  for (int j = 0; j < MAX_ITERTER; j++) {
    Doub dx = (xl - rts) * f / (f - fl);
    xl = rts;
    fl = f;
    rts += dx;
    f = func(rts);
    x_k.push_back(rts); // For table
    if (abs(dx) < xacc || f == 0.0) {
      print_table(x_k, RootType::Secant); // Print table
      return rts;
    }
  }
  throw("Maximum number of iterations exceeded in rtsec");
}

template <class T>
Doub root_false_position(T &func, const Doub x1, const Doub x2,
                         const Doub xacc) {
  std::vector<double> x_k; // For table
  const Int MAXIT = 30;
  Doub xl, xh, del;
  Doub fl = func(x1);
  Doub fh = func(x2);
  if (fl * fh > 0.0)
    throw("Root must be bracketed in rtflsp");
  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    xl = x2;
    xh = x1;
    SWAP(fl, fh);
  }
  Doub dx = xh - xl;
  for (Int j = 0; j < MAXIT; j++) {
    Doub rtf = xl + dx * fl / (fl - fh);
    Doub f = func(rtf);
    x_k.push_back(rtf); // For table
    if (f < 0.0) {
      del = xl - rtf;
      xl = rtf;
      fl = f;
    } else {
      del = xh - rtf;
      xh = rtf;
      fh = f;
    }
    dx = xh - xl;
    if (abs(del) < xacc || f == 0.0) {
      print_table(x_k, RootType::FalsePositive); // Print table
      return rtf;
    }
  }
  throw("Maximum number of iterations exceeded in rtflsp");
}

template <class T>
Doub root_newton(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
  std::vector<double> x_k; // For table
  const Int MAX_ITER = 20;
  Doub rtn = 0.5 * (x1 + x2);
  x_k.push_back(rtn); // For table
  for (Int j = 0; j < MAX_ITER; j++) {
    Doub f = funcd(rtn);
    Doub df = funcd.df(rtn);
    Doub dx = f / df;
    rtn -= dx;
    x_k.push_back(rtn); // For table
    if ((x1 - rtn) * (rtn - x2) < 0.0)
      throw("Jumped out of brackets in rtnewt");
    if (abs(dx) < xacc) {
      print_table(x_k, RootType::Newton); // Print table
      return rtn;
    }
  }
  throw("Maximum number of iterations exceeded in rtnewt");
}
} // namespace roots_table
#endif // !ROOTS_TABLE_H
