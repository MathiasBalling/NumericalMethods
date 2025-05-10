#include "nr3.h"
#include "utils/multi_roots.h" // Modified newton method
#include "utils/utilities.h"
#include <cassert>
#include <cmath>
#include <print>

// Material constants:
const double v = 120.0;                 // kg
const double k = 2.5;                   // m
const double w = 4.0;                   // kg/m
const double alpha = 2.0 * pow(10, -7); // kg^-1

// Other constants:
const double d = 30.0; // m
double n = 5.0;        // m (Updated in main)

VecDoub vecfunc(VecDoub_I q) {
  assert(q.size() == 8);

  // q=[L0, L, p, x, theta, phi/varphi, a, H]
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

int main() {

  {
    std::println("\n{:/^120}", " n=5.0 ");
    n = 5.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {29.0, 29.1, 5.0, 15.0, 1.0, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }

  {
    std::println("\n{:/^120}", " n=2.0 ");
    n = 2.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {27.0, 27.1, 2.0, 14.0, 1.0, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }

  {
    std::println("\n{:/^120}", " n=1.0 ");
    n = 1.0;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {26.0, 26.1, 1.0, 13.0, 0.5, 0.5, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }

  {
    std::println("\n{:/^120}", " n=0.5 ");
    n = 0.5;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {26.0, 26.1, 0.5, 13.0, 0.2, 0.1, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }

  {
    std::println("\n{:/^120}", " n=0.2 ");
    n = 0.2;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {25.5, 25.6, 0.2, 12.0, 0.1, 0.05, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }

  {
    std::println("\n{:/^120}", " n=0.1 ");
    n = 0.1;
    // q=[L0, L, p, x, theta, phi/varphi, a, H]
    Doub qvals[8] = {25.5, 25.6, 0.1, 12.0, 0.1, 0.05, 40.0, 100.0};
    VecDoub_IO q(8, qvals);
    util::print(q, "q initial guess");
    bool check;
    auto func = vecfunc;
    newton_multi(q, check, func);
    util::print(q, "q roots with newton");
    std::println("Local minimum: {}", check);
  }
}
