#include "utils/pde.h"

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
  pde_table(N, lambda, f, u_boundary, 1e-5, 100);
}
