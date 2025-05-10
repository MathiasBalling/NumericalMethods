#ifndef PDE_H
#define PDE_H
#include "banded.h"
#include "nr3.h"
#include "utils/errors.h"
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

void pde_table(const int starting_steps, const double lambda,
               double f(double x, double y),
               double u_boundary(double x, double y), double accuracy = 1e-5,
               int max_N = 200) {
  const double low = 0.0;
  const double high = 1.0;
  std::vector<double> A_k;
  const double expected_order = 2.0;
  bool should_stop = false;
  size_t N = starting_steps;
  while (!should_stop) {
    const double res = pde(N, lambda, low, high, f, u_boundary);
    A_k.push_back(res);

    if (A_k.size() >= 2) {
      const auto error =
          richardson_extrapolation_error_current(A_k, pow(2, expected_order));
      if (error < accuracy || N >= max_N) {
        should_stop = true;
      }
    }
    N *= 2;
  }

  // Calculate the differences
  std::vector<double> A_diff_k = {NAN};
  for (size_t i = 1; i < A_k.size(); i++) {
    double diff = A_k[i - 1] - A_k[i];
    A_diff_k.push_back(diff);
  }

  // Calculate the orders
  const auto alpha_k_computed = alpha_k_order_computed(A_k);

  // Calculate the richardson extrapolation
  const auto rich_error = richardson_extrapolation_error(
      A_k, pow(2, expected_order)); // pow(2, expected_order) as we use N*=2

  const auto order_estimate = compute_order_estimate(alpha_k_computed);

  // Print the table
  // Table header
  std::println("|{:^6}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|", "N", "A(N)",
               "A(N/2)-A(N)", "alpha^k", "Rich error", "Order est.");
  // Table header separator
  std::println("|{:-^6}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|", "", "", "",
               "", "", "");

  // Table body
  auto N_A = starting_steps;
  for (size_t i = 0; i < A_k.size(); i++) {
    if (i == 0) {
      std::println("|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|",
                   N_A, A_k.at(i), "", "", "", "");
    } else if (i == 1) {
      std::println("|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|",
                   N_A, A_k.at(i), A_diff_k.at(i), "", rich_error.at(i), "");
    } else {
      std::println("|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|",
                   N_A, A_k.at(i), A_diff_k.at(i), alpha_k_computed.at(i),
                   rich_error.at(i), order_estimate.at(i));
    }
    N_A *= 2;
  }
}

#endif // PDE_H
