#include "nr3.h"
#include "utils/parabolic_pde.h"
#include <numbers>
#include <print>

double u(double x, double t) {
  if (t == 0.0) {
    return pow(x, 2);
  }
  if (x == 0.0) {
    return 0.0;
  }
  if (x == 1.0) {
    return 1 + sin(t);
  }

  assert(false && "Not implemented");
  return 0;
}

double f(double x, double t) { return sin(std::numbers::pi * x) * exp(-t); }
int main() {
  // 5b
  {
    // 1)
    std::println("\n{:/^100}", " i ");
  }

  {
    // 2)
    std::println("\n{:/^100}", " ii ");
    double alpha = 4.0;
    double x_target = 0.5;
    double t_target = 10.0;
    size_t N = (1.0 - 0.0) / 0.25;
    std::println("N: {}", N);

    parabolic_pde_table(N, alpha, x_target, t_target, f, u, 1e-7, 100000);
  }
}
