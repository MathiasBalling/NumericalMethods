#include "nr3.h"
#include "utils/multi_roots.h"
#include "utils/utilities.h"
#include <cassert>
#include <print>

VecDoub vecfunc(VecDoub_I x) {
  assert(x.size() == 4);

  VecDoub f(4);

  f[0] = 3 * x[0] + x[1] * sin(x[2]) - cos(x[0]) + cos(pow(x[1], 2)) + 4.2;
  f[1] = 3 * x[1] + x[0] * x[2] * x[3] + sin(x[1]) - 5.1;
  f[2] = -pow(x[1], 2) + x[2] * pow(x[3], 2) + 3 * x[2] + 5.2;
  f[3] = x[0] + 3 * x[3] + sin(pow(x[2], 2) * pow(x[3], 2)) + cos(x[1]) - 2.3;
  return f;
}
int main() {

  {
    // 1)
    std::println("\n{:/^100}", " i ");
    double vals[4] = {-0.7, 1.2, 2.3, -4.1};
    VecDoub x(4, vals);
    VecDoub res = vecfunc(x);
    util::print(res, "res");
  }

  {
    // 2 & 3)
    std::println("\n{:/^100}", " ii & iii ");
    double vals[4] = {0.0, 0.0, 0.0, 0.0};
    VecDoub x(4, vals);
    bool check;
    // TODO: Check for backtracking
    newton_multi_table(x, check, vecfunc, 7, 0);
    std::println("Local minimum = {}", check);
  }

  return 0;
}
