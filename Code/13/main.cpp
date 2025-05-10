#include "pde.h"

// Function to be integrated
double f(double x, double y) { return 1 + x + y; }

// Boundary conditions
double u_boundary(double x, double y) {
  (void)x;
  (void)y;
  return 0;
}

int main() {
  const size_t N = 4;
  const double lambda = 0.0;
  const auto res = pde(N, lambda, 0.0, 1.0, f, u_boundary);
  std::println("res: {}", res);
}
