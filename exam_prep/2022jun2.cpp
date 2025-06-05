#include "nr3.h"
#include "roots_multidim.h"
#include "utils/multi_roots.h"
#include "utils/utilities.h"
#include <cassert>
#include <print>

// Constants:
const double a1 = 1.1, a2 = 2.1, a3 = 0.8;
const double b1 = 0.4, b2 = 1.3, b3 = 0.5;
const double epsilon = 1e-8;

VecDoub r_A(const double u_A) {
  VecDoub r(3);
  r[0] = a1 * pow(cos(1 + u_A), 3);
  r[1] = a2 * pow(u_A, 2);
  r[2] = a3 * u_A * sin(u_A);
  return r;
}

VecDoub r_A_prime(const double u_A) {
  return (r_A(u_A + epsilon) - r_A(u_A)) / epsilon;
}

VecDoub r_B(const double u_B) {
  VecDoub r(3);
  r[0] = b1 * (u_B + exp(-pow(u_B, 2)));
  r[1] = b2 * pow(u_B, 3);
  r[2] = b3 * cos(u_B);
  return r;
}

VecDoub r_B_prime(const double u_B) {
  return (r_B(u_B + epsilon) - r_B(u_B)) / epsilon;
}

VecDoub vecfunc(VecDoub_I u) {
  assert(u.size() == 2);
  double u_A = u[0];
  double u_B = u[1];

  // Zero vectors
  VecDoub f(2);
  f[0] = r_A_prime(u_A) * (r_A(u_A) - r_B(u_B));
  f[1] = -r_B_prime(u_B) * (r_A(u_A) - r_B(u_B));
  return f;
}

int main() {
  VecDoub u(2);
  u[0] = 1.0;
  u[1] = -1.0;

  {
    // 1)
    VecDoub f = vecfunc(u);
    util::print(f, "f");
  }
  {
    // 2)
    bool check;
    auto res = u;
    newton_multi_table(res, check, vecfunc, 10, 0);
    std::println("check: {}", check);
    //|    it    |        value        |      diff
    //|    10    |   0.061582338781    | -9.09223552562e-05
    //|    10    |   -0.457511911535   |  0.000133120369046
  }
  {
    // 3)
    bool check;
    auto res = u;
    newton_multi_table(res, check, vecfunc, 200, 1.0e-4);
    std::println("check: {}", check);
    // Look at table. After 8 its the accuracy is better than 1e-4
  }
  return 0;
}
