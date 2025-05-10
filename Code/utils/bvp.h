#ifndef BVP_H
#define BVP_H
#include "nr3.h"
#include "tridag.h"
#include "utils/errors.h"

std::vector<double> target_values(std::vector<VecDoub> y_estimates,
                                  const double target_x, const double x_low,
                                  const double x_high) {
  const auto size = y_estimates.size();
  std::vector<double> estimates;
  for (size_t i = 0; i < size; i++) {
    const size_t N = y_estimates[i].size() - 1;
    const double h = (x_high - x_low) / (double)N;
    const size_t index = (target_x - x_low) / h;
    estimates.push_back(y_estimates[i][index]);
  }
  return estimates;
}

// y(target_x) = ?
std::vector<VecDoub>
finite_difference_method(const int starting_N, const double target_x,
                         const double a, const double b, const double alpha,
                         const double beta,
                         double F(double y_prime, double y, double x),
                         double F_y(double y_prime, double y, double x),
                         double F_y_prime(double y_prime, double y, double x),
                         double accuracy, size_t max_N) {
  std::vector<VecDoub> y_estimates;
  int N = starting_N;

  // Initial guess y_i using linear interpolation
  // N=2
  // y(alpha,x,beta)
  VecDoub y(N + 1);
  for (int i = 0; i < N + 1; i++) {
    y[i] = alpha + (beta - alpha) / (double)N * (double)i;
  }

  bool should_stop = false;
  while (!should_stop) {
    const double h = (b - a) / (double)N;
    // Define J (using tridag)
    VecDoub J_a(N - 1), J_b(N - 1), J_c(N - 1), J_u_solver(N - 1),
        J_r_solver(N - 1);
    {
      // First element (using alpha)
      double x_1 = a + h;
      J_a[0] = 0.0; // Not used
      J_b[0] = 2.0 + pow(h, 2) * F_y((y[2] - alpha) / (2.0 * h), y[1], x_1);
      J_c[0] =
          -1.0 + (h / 2.0) * F_y_prime((y[2] - alpha) / (2.0 * h), y[1], x_1);

      // Last element (using beta)
      double x_N_1 = b - h;
      J_a[N - 2] = -1. - (h / 2.0) * F_y_prime((beta - y[N - 1]) / (2.0 * h),
                                               y[N - 1], x_N_1);
      J_b[N - 2] =
          2.0 + pow(h, 2) * F_y((beta - y[N - 1]) / (2.0 * h), y[N - 1], x_N_1);
      J_c[N - 2] = 0.0; // Not used

      // Middle elements
      for (int i = 1; i < N - 2; i++) {
        double x_i = a + (i + 1) * h;
        J_a[i] =
            -1. -
            h / 2.0 * F_y_prime((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
        J_b[i] =
            2.0 + pow(h, 2) * F_y((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
        J_c[i] =
            -1.0 +
            h / 2.0 * F_y_prime((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
      }
    }

    // Calculate phi(y)=0
    {
      double x_1 = a + h;
      J_r_solver[0] = -alpha + 2.0 * y[1] - y[2] +
                      pow(h, 2) * F((y[2] - alpha) / (2.0 * h), y[1], x_1);
      double x_N_1 = b - h;
      J_r_solver[N - 2] =
          -y[N - 2] + 2.0 * y[N - 1] - beta +
          pow(h, 2) * F((beta - y[N - 2]) / (2.0 * h), y[N - 1], x_N_1);
      for (int i = 2; i < N - 1; i++) {
        double x_i = a + i * h;
        J_r_solver[i - 1] =
            -y[i - 1] + 2.0 * y[i] - y[i + 1] +
            pow(h, 2) * F((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
      }

      // J_r = -phi(y)
      for (int i = 0; i < N - 1; i++) {
        J_r_solver[i] = -J_r_solver[i];
      }
    }

    // Use tridag to calculate delta y
    tridag(J_a, J_b, J_c, J_r_solver, J_u_solver);

    // util::print(y, "y");
    for (int i = 1; i < N; i++) {
      y[i] += J_u_solver[i - 1];
    }

    // Save the current y
    y_estimates.push_back(y);

    // Check if we should stop
    if (y_estimates.size() > 2) {
      const auto target_estimates = target_values(y_estimates, target_x, a, b);

      // Expected order is 2, then since N*=2, alpha^k is then pow(2,2)
      const auto error =
          richardson_extrapolation_error_current(target_estimates, pow(2, 2));
      // Use richardson error to stop or use max_N
      if (abs(error) < accuracy || N >= max_N) {
        should_stop = true;
      }
    }

    // Half the step size and fill in missing y values with linear interpolation
    // for the next round
    N *= 2;
    VecDoub y_new(N + 1);
    y_new[0] = alpha;
    y_new[N] = beta;
    for (int i = 1; i < N; i++) {
      if (i % 2 == 0) {
        y_new[i] = y[i / 2];
      } else {
        y_new[i] = 0.5 * (y[(i - 1) / 2] + y[(i + 1) / 2]);
      }
    }

    y = y_new;
  }

  return y_estimates;
}

void finite_difference_method_table(
    const int starting_N, const double target_x, const double a, const double b,
    const double alpha, const double beta,
    double F(double y_prime, double y, double x),
    double F_y(double y_prime, double y, double x),
    double F_y_prime(double y_prime, double y, double x),
    double accuracy = 1e-4, size_t max_N = 1000) {
  const double expected_order = 2.0;
  const auto y_estimates =
      finite_difference_method(starting_N, target_x, a, b, alpha, beta, F, F_y,
                               F_y_prime, accuracy, max_N);
  const auto A_k = target_values(y_estimates, target_x, a, b);
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

  const std::vector<int> f_comps(A_k.size());

  // Print the table
  // Table header
  std::println("|{:^6}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|{:^10}|", "i", "A(i)",
               "A(i-1)-A(i)", "alpha^k", "Rich error", "Order est.", "f comps");
  // Table header separator
  std::println("|{:-^6}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^10}|", "",
               "", "", "", "", "", "");

  // Table body
  auto N = starting_N;
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

#endif // !BVP_H
