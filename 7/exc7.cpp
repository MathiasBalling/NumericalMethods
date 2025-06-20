#include "nr3.h"
// #include "roots_multidim.h"
#include "utils/multi_roots.h"
#include "utils/utilities.h"
#include <cassert>
#include <print>

VecDoub vecfunc(VecDoub_I x) {
  assert(x.size() == 2);

  VecDoub f(2);
  f[0] = x[0] + 2 * sin(x[1] - x[0]) - exp(-sin(x[1] + x[0]));
  f[1] = x[0] * cos(x[1]) + sin(x[0]) - 1;
  return f;
}

int main() {
  ////////////////////// i //////////////////////
  std::println("\n{:/^100}", " i ");
  {
    Doub x[] = {1.0, 1.0};
    auto f = vecfunc(VecDoub(2, x));
    util::print(f);
  }

  ////////////////////// ii //////////////////////
  std::println("\n{:/^100}", " ii ");
  {
    std::println("I would use newton for this problem because...");
  }

  ////////////////////// iii //////////////////////
  std::println("\n{:/^100}", " iii ");
  {
    Doub xvals[] = {1.0, 2.0};
    VecDoub_IO x(2, xvals);
    bool check;
    auto func = vecfunc;
    newton_multi_table(x, check, func, 200, 1.0e-8);
    util::print(x, "x roots with newton");
    std::println("Local minimum = {}", check);
  }

  ////////////////////// iv //////////////////////
  std::println("\n{:/^100}", " iv ");
  {
  }

  ////////////////////// v //////////////////////
  std::println("\n{:/^100}", " v ");
  {
  }

  ////////////////////// vi //////////////////////
  std::println("\n{:/^100}", " vi ");
  {
    std::println("Globally convergent Newton method page 477-482");
  }

  ////////////////////// vii //////////////////////
  std::println("\n{:/^100}", " vii ");
  {
    Doub xvals[] = {1.0, 1.0};
    VecDoub x(2, xvals);
    bool check;
    auto func = vecfunc;
    // newt(x, check, func);
    util::print(x, "x roots with newton");
    std::println("Local minimum = {}", check);
  }

  return 0;
}
