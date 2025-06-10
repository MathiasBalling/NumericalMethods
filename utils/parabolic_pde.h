#ifndef PARABOLIC_PDE_H
#define PARABOLIC_PDE_H
#include "nr3.h"
#include "tridag.h"
#include "utils/errors.h"
#include <cmath>

static size_t f_comps_current = 0;

double target_value(VecDoub u_estimates, const double x_target) {
  const auto N = u_estimates.size();
  const double dx = (1.0 - 0.0) / (double)N;
  const int index = (int)((x_target - 0.0) / dx);
  return u_estimates[index];
}

VecDoub parabolic_pde(const size_t N, const double alpha, const double t_target,
                      double f(double x, double t),
                      double u(double x, double t)) {
  double h = (1.0 - 0.0) / (double)N;
  // Both for Δx and Δt
  double dt = h, dx = h;
  double r = alpha * (dt / pow(dx, 2));

  // Boundary values
  auto g = [&](double x) { return u(x, 0.0); };
  auto a = [&](double t) { return u(0.0, t); };
  auto b = [&](double t) { return u(1.0, t); };

  VecDoub J_a(N - 1), J_b(N - 1), J_c(N - 1), J_u(N - 1), J_r(N - 1);
  VecDoub u_cur(N + 1);
  for (int j = 0; j < N + 1; j++) {
    const double x = j * dx;
    u_cur[j] = g(x);
  }

  size_t t_steps = t_target / dt;
  for (size_t n = 0; n <= t_steps; n++) {
    double t_n = n * dt;
    for (int j = 1; j < N; j++) {
      const double x = j * dx;
      J_a[j - 1] = -0.5 * r;
      J_b[j - 1] = 1 + r;
      J_c[j - 1] = -0.5 * r;
      J_r[j - 1] = 0.5 * r * u_cur[j - 1] + (1 - r) * u_cur[j] +
                   0.5 * r * u_cur[j + 1] +
                   (dt / 2.0) * (f(x, t_n + dt) + f(x, t_n));
      f_comps_current += 2;
    }

    tridag(J_a, J_b, J_c, J_r, J_u);

    // Update u_cur for next step
    for (int j = 0; j < N - 1; j++) {
      u_cur[j + 1] = J_u[j];
    }
    u_cur[0] = a(t_n + dt);
    u_cur[N] = b(t_n + dt);
  }

  return u_cur;
}

void parabolic_pde_table(const size_t N_start, const double alpha,
                         const double x_target, const double t_target,
                         double f(double x, double t),
                         double u(double x, double t),
                         const double accuracy = 1e-4,
                         const size_t max_N = 200) {
  const double expected_order = 2;
  std::vector<double> A_k;
  std::vector<size_t> f_comps;
  bool should_stop = false;
  size_t N_cur = N_start;
  while (!should_stop) {
    const auto u_estimates = parabolic_pde(N_cur, alpha, t_target, f, u);

    // Handle f computations
    f_comps.push_back(f_comps_current);
    f_comps_current = 0;

    const auto estimate = target_value(u_estimates, x_target);
    A_k.push_back(estimate);
    N_cur *= 2;
    if (A_k.size() > 2) {
      const auto error =
          richardson_extrapolation_error_current(A_k, pow(2, expected_order));
      std::println("error: {}", error);
      should_stop = abs(error) < accuracy || N_cur >= max_N;
    }
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
  std::println("|{:^6}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|{:^10}|", "N", "A(N)",
               "A(N/2)-A(N)", "alpha^k", "Rich error", "Order est.", "f comps");
  // Table header separator
  std::println("|{:-^6}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^10}|", "",
               "", "", "", "", "", "");

  // Table body
  auto N = N_start;
  for (size_t i = 0; i < A_k.size(); i++) {
    if (i == 0) {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|", N,
          A_k.at(i), "", "", "", "", f_comps.at(i));
    } else if (i == 1) {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|", N,
          A_k.at(i), A_diff_k.at(i), "", rich_error.at(i), "", f_comps.at(i));
    } else {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|", N,
          A_k.at(i), A_diff_k.at(i), alpha_k_computed.at(i), rich_error.at(i),
          order_estimate.at(i), f_comps.at(i));
    }
    N *= 2;
  }
}

#endif // PARABOLIC_PDE_H
