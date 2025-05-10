#ifndef PDE_H
#define PDE_H
#include "banded.h"
#include "nr3.h"
#include <cassert>
#include <cmath>
#include <cstddef>

// Index functions
// j is column, k is row
size_t idx(size_t j, size_t k, size_t N) {
  assert(k < N);
  assert(j < N);
  return N * k + j;
}

// Returns j and k
// j is column, k is row
std::pair<size_t, size_t> inv_idx(size_t i, size_t N) {
  size_t k = i / N;
  size_t j = i % N;
  return {j, k};
}
// Dirichlet condition: Fixed value at boundary.
// Neumann condition: Fixed derivative at boundary.
// This function uses the Dirichlet condition at the boundary.
double pde(const size_t N, const double lambda, const double low,
           const double high, double f(double x, double y),
           double u_boundary(double x, double y)) {

  // Boundary values
  auto a0 = [&](double x) { return u_boundary(x, 0.0); };
  auto a1 = [&](double x) { return u_boundary(x, 1.0); };
  auto b0 = [&](double y) { return u_boundary(0.0, y); };
  auto b1 = [&](double y) { return u_boundary(1.0, y); };

  const int n = N - 1; // We won't use the boundary points a x or y = high
  // N=4 gives (N-1)^2 = 9 points inside the domain
  const int total_points = n * n;

  MatDoub A(total_points, total_points, 0.0);
  VecDoub phi(total_points, 0.0);
  VecDoub omega(total_points, 0.0);
  const double h = (high - low) / (double)N;
  const double h2 = h * h;
  // row is k
  for (size_t k = 0; k < n; k++) {
    // col is j
    for (size_t j = 0; j < n; j++) {
      const size_t index = idx(j, k, n);
      const double x = low + (double)(j + 1) * h;
      const double y = low + (double)(k + 1) * h;

      A[index][index] = 4 + pow(h, 2) * lambda;
      phi[index] = h2 * f(x, y); // Boundary condition values added later

      // Left neighbor
      if (j > 0) {
        A[index][idx(j - 1, k, n)] = -1;
      } else {
        phi[index] += b0(y); // left boundary
      }

      // Right neighbor
      if (j < n - 1) {
        A[index][idx(j + 1, k, n)] = -1;
      } else {
        phi[index] += b1(y); // right boundary
      }

      // Bottom neighbor
      if (k > 0) {
        A[index][idx(j, k - 1, n)] = -1;
      } else {
        phi[index] += a0(x); // bottom boundary
      }

      // Top neighbor
      if (k < n - 1) {
        A[index][idx(j, k + 1, n)] = -1;
      } else {
        phi[index] += a1(x); // top boundary
      }
    }
  }

  // Convert to band-diagonal system
  MatDoub A_bandec(total_points, 2 * n + 1);
  for (int i = 0; i < total_points; i++) {
    A_bandec[i][n] = A[i][i];
    for (int j = -n; j <= n; j++) {
      const int A_index = i + j;
      const int A_bandec_index = j + n;
      if (A_index >= 0 && A_index <= total_points) {
        A_bandec[i][A_bandec_index] = A[i][A_index];
      }
    }
  }

  // Setup Bandec
  Bandec band(A_bandec, n, n);
  // Vector for storing solution
  VecDoub sol(total_points);
  // Solve
  band.solve(phi, sol);

  // find u(0.5, 0.5)
  size_t center_idx = idx((n - 1) / 2, (n - 1) / 2, n);
  return sol[center_idx];
}

#endif // PDE_H
