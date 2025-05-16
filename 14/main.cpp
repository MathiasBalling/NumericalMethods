#include "utils/parabolic_pde.h"
#include <cassert>
#include <print>

double u(double x, double t) {
  if (t == 0.0) {
    return pow(x, 4);
  }
  if (x == 0.0) {
    return 0.0;
  }
  if (x == 1.0) {
    return 1.0;
  }

  assert(false && "Not implemented");
  return 0;
}
double f(double x, double t) { return x * (1 - x) * cos(t) * exp(-t / 10.0); }

int main() {
  double alpha = 1.0;
  size_t N = 2;
  double x_target = 0.5;
  double t_target = 20.0;
  parabolic_pde_table(N, alpha, x_target, t_target, f, u, 1e-4, 1000);
}
