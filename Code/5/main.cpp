#include "roots.h"
#include <numbers>
#include <print>

double bisection(double (*func)(double), double xl, double xh,
                 double accuracy) {
  const size_t max_iter = 50;
  Doub dx, x_mid, rtb;
  Doub f_xl = func(xl);
  Doub f_xh = func(xh);
  if (f_xl * f_xh >= 0.0)
    throw("Root must be bracketed for bisection in rtbis");

  std::println("|{:^10}|{:^20}|{:^20}|{:^20}|", "k", "x_k", "d_k", "rtb");

  rtb = f_xl < 0.0 ? (dx = xh - xl, xl) : (dx = xl - xh, xh);
  for (int i = 0; i < max_iter; i++) {
    auto xmid_old = x_mid;
    f_xh = func(x_mid = rtb + (dx *= 0.5));
    std::println("|{:^10}|{:^20.10}|{:^20.10}|{:^20.10}|", i + 1, x_mid,
                 x_mid - xmid_old, rtb);
    if (f_xh <= 0.0)
      rtb = x_mid;
    if (abs(dx) < accuracy || f_xh == 0.0)
      return rtb;
  }
  throw("Too many bisections in rtbis");
}

template <class T>
Doub ridder(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  const int MAXIT = 60;
  Doub fl = func(x1);
  Doub fh = func(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    Doub xl = x1;
    Doub xh = x2;
    Doub ans = -9.99e99;

    // For printing table
    std::println("|{:^10}|{:^20}|{:^20}|", "k", "x_k", "d_k");

    Doub xm = 0.5 * (xl + xh);
    Doub xm_old = xm;
    for (int j = 0; j < MAXIT; j++) {
      xm_old = xm;
      xm = 0.5 * (xl + xh);
      Doub f_mid = func(xm);
      Doub s = sqrt(f_mid * f_mid - fl * fh);
      if (j > 0) {
        std::println("|{:^10}|{:^20.10}|{:^20.10}|", j, xm, xm - xm_old);
      }
      if (s == 0.0)
        return ans;
      Doub xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * f_mid / s);
      if (abs(xnew - ans) <= xacc)
        return ans;
      ans = xnew;
      Doub fnew = func(ans);
      if (fnew == 0.0)
        return ans;
      if (SIGN(f_mid, fnew) != f_mid) {
        xl = xm;
        fl = f_mid;
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
      if (abs(xh - xl) <= xacc)
        return ans;
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
Doub secant(T &func, const Doub x1, const Doub x2, const Doub xacc) {
  const int MAXIT = 30;
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
  for (int j = 0; j < MAXIT; j++) {
    Doub dx = (xl - rts) * f / (f - fl);
    xl = rts;
    fl = f;
    rts += dx;
    f = func(rts);
    if (abs(dx) < xacc || f == 0.0)
      return rts;
  }
  throw("Maximum number of iterations exceeded in rtsec");
}

double f(double x) { return x - cos(x); }

int main(int argc, char *argv[]) {
  // Bisection
  println("\n{:/^100}", " Bisection ");
  auto res_bisection = bisection(f, 0, std::numbers::pi / 2, 0.00000001);
  // From roots.h:
  std::println("Mine={} roots.h={}", res_bisection,
               rtbis(f, 0, std::numbers::pi / 2, 0.00000001));

  // Ridder
  println("\n{:/^100}", " Ridder ");
  auto res_ridder = ridder(f, 0, std::numbers::pi / 2, 0.000000000000001);
  // From roots.h:
  std::println("Mine={} roots.h={}", res_ridder,
               zriddr(f, 0, std::numbers::pi / 2, 0.000000000000001));

  // Secant
  println("\n{:/^100}", " Secant ");
  auto res_secant = secant(f, 0, std::numbers::pi / 2, 0.00000001);
  // From roots.h:
  std::println("Mine={} roots.h:{}", res_secant,
               rtsec(f, 0, std::numbers::pi / 2, 0.00000001));

  return 0;
}
