#ifndef ODE_H
#define ODE_H

#include "multi_roots_table.h"
#include "nr3.h"
#include "utilities.h"

enum class ODEMethod { Euler, Midpoint, FORTH_ORDER, TRAPEZOIDAL, LEAP_FROG };

// Too keep track of f-comps
static size_t f_comps_current = 0;
// Euler method
VecDoub
first_order_runge_kuttea_method(double low, double high, int steps, VecDoub_I y,
                                VecDoub derivs(const Doub x, VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;

  VecDoub y_n = y;
  for (double x = low; x < high; x += h) {
    f_comps_current++;
    auto y_n_next = y_n + h * derivs(x, y_n);
    y_n = y_n_next;
  }
  return y_n;
}

// Midpoint method
VecDoub second_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &y,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    f_comps_current += 2;
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto y_n_next = y_n + k2;
    y_n = y_n_next;
  }
  return y_n;
}

VecDoub fourth_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &y,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    f_comps_current += 4;
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto k3 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k2);
    auto k4 = h * derivs(x_n + h, y_n + k2);
    auto y_n_next = y_n + 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
    y_n = y_n_next;
  }
  return y_n;
}

VecDoub trapezoidal(double x_low, double x_high, int steps, VecDoub_I &y,
                    VecDoub derivs(const Doub x, VecDoub_I &y)) {
  // Trapezoidal method
  const double h = (x_high - x_low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = x_low; x_n < x_high; x_n += h) {
    // Euler step
    f_comps_current += 1;
    auto y_n_next = y_n + h * derivs(x_n, y_n);

    auto func = [&](VecDoub y_guess) {
      f_comps_current += 2;
      return y_guess - y_n -
             h / 2 * (derivs(x_n, y_n) + derivs(x_n + h, y_guess));
    };

    bool check;
    newton_multi(y_n_next, check, func, 200, 1.0e-7);
    y_n = y_n_next;
  }

  return y_n;
}

VecDoub leap_frog_method(double low, double high, int steps, VecDoub_I &y,
                         VecDoub derivs(const Doub x, VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  f_comps_current += 2;
  VecDoub y_n_last = y + h * derivs(low, y);
  VecDoub y_n = y_n_last + h * derivs(low + h, y_n_last);
  for (double x_n = low + h + h; x_n < high; x_n += h) {
    f_comps_current++;
    VecDoub y_n_next = y_n_last + 2 * h * derivs(x_n, y_n);
    y_n_last = y_n;
    y_n = y_n_next;
  }
  return y_n;
}

std::vector<VecDoub> alpha_k_order_computed(const std::vector<VecDoub> &A_k) {
  const auto vec_size = A_k[0].size();
  std::vector<VecDoub> alpha_k = {
      VecDoub(vec_size, NAN),
      VecDoub(vec_size, NAN)}; // No error on the first two

  for (size_t i = 2; i < A_k.size(); i++) {
    VecDoub alpha_i(vec_size);
    for (size_t j = 0; j < vec_size; j++) {
      const double A_1 = A_k[i - 2][j];
      const double A_2 = A_k[i - 1][j];
      const double A_3 = A_k[i][j];
      alpha_i[j] = (A_1 - A_2) / (A_2 - A_3);
    }
    alpha_k.push_back(alpha_i);
  }
  return alpha_k;
}

std::vector<VecDoub>
richardson_extrapolation_error(const std::vector<VecDoub> &A_k,
                               const double alpha_k_order_expected) {
  const auto vec_size = A_k[0].size();
  std::vector<VecDoub> A_R = {VecDoub(vec_size, NAN)}; // No error on the first
  for (size_t i = 1; i < A_k.size(); i++) {
    VecDoub error(vec_size);
    for (size_t j = 0; j < vec_size; j++) {
      const double A_1 = A_k[i - 1][j];
      const double A_2 = A_k[i][j];
      error[j] = (A_2 - A_1) / (alpha_k_order_expected - 1);
    }
    A_R.push_back(error);
  }
  return A_R;
}

VecDoub
richardson_extrapolation_error_current(const std::vector<VecDoub> &A_k,
                                       const double alpha_k_order_expected) {
  VecDoub errors = VecDoub(A_k[0].size(), std::numeric_limits<double>::max());

  for (size_t i = 0; i < A_k[0].size(); i++) {
    const double A_1 = A_k[A_k.size() - 2][i];
    const double A_2 = A_k[A_k.size() - 1][i];
    // WARN:  Only use expected if compute_order_estimate is good.
    errors[i] = (A_2 - A_1) / (alpha_k_order_expected - 1);
  }
  return errors;
}

std::vector<VecDoub>
compute_order_estimate(const std::vector<VecDoub> &alpha_k_computed) {
  const auto vec_size = alpha_k_computed[0].size();
  std::vector<VecDoub> order_estimate = {
      VecDoub(vec_size, NAN),
      VecDoub(vec_size, NAN)}; // No order for first two entries
  for (size_t i = 2; i < alpha_k_computed.size(); i++) {
    VecDoub p(vec_size);
    for (size_t j = 0; j < vec_size; j++) {
      p[j] = log2(abs(alpha_k_computed[i][j]));
    }
    order_estimate.push_back(p);
  }
  return order_estimate;
}

void print_ode_table(const double low, const double high,
                     const int starting_steps, VecDoub_I &y,
                     VecDoub derivs(const Doub x, VecDoub_I &y),
                     const ODEMethod type, const double accuracy = 1e-5,
                     const size_t max_steps = 20000) {
  std::vector<VecDoub> A_k;
  std::vector<int> f_comps;
  double expected_order = 0.0;
  int its = starting_steps;
  bool should_stop = false;
  while (!should_stop) {
    f_comps_current = 0;
    auto y_n = y;
    switch (type) {
    case ODEMethod::Euler:
      A_k.push_back(
          first_order_runge_kuttea_method(low, high, its, y_n, derivs));
      f_comps.push_back(f_comps_current);
      expected_order = 1.0;
      break;
    case ODEMethod::Midpoint:
      A_k.push_back(
          second_order_runge_kuttea_method(low, high, its, y_n, derivs));
      f_comps.push_back(f_comps_current);
      expected_order = 2.0;
      break;
    case ODEMethod::FORTH_ORDER:
      A_k.push_back(
          fourth_order_runge_kuttea_method(low, high, its, y_n, derivs));
      f_comps.push_back(f_comps_current);
      expected_order = 4.0;
      break;
    case ODEMethod::TRAPEZOIDAL:
      A_k.push_back(trapezoidal(low, high, its, y_n, derivs));
      f_comps.push_back(f_comps_current);
      expected_order = 2.0;
      break;
    case ODEMethod::LEAP_FROG:
      A_k.push_back(leap_frog_method(low, high, its, y_n, derivs));
      f_comps.push_back(f_comps_current);
      expected_order = 2.0;
      break;
    default:
      throw("Unknown integration type");
    }

    if (A_k.size() >= 2) {
      const auto errors =
          richardson_extrapolation_error_current(A_k, pow(2, expected_order));
      double biggest_error = 0.0;
      for (size_t i = 0; i < errors.size(); i++) {
        if (abs(errors[i]) > biggest_error) {
          biggest_error = abs(errors[i]);
        }
      }
      if (biggest_error < accuracy || its >= max_steps) {
        should_stop = true;
      }
    }
    its *= 2.0;
  }

  const auto vec_size = y.size();

  // Calculate the differences
  std::vector<VecDoub> A_diff_k = {VecDoub(vec_size, NAN)};
  for (size_t i = 1; i < A_k.size(); i++) {
    VecDoub diff(vec_size);
    for (size_t j = 0; j < vec_size; j++) {
      diff[j] = A_k[i - 1][j] - A_k[i][j];
    }
    A_diff_k.push_back(diff);
  }

  // Calculate the orders
  auto alpha_k_computed = alpha_k_order_computed(A_k);

  // Calculate the richardson extrapolation
  auto rich_error = richardson_extrapolation_error(
      A_k, pow(2, expected_order)); // pow(2, expected_order) as we use
                                    // N-1=1,2,4,8

  auto order_estimate = compute_order_estimate(alpha_k_computed);

  // Print the table
  for (size_t i = 0; i < vec_size; i++) {
    // Table header
    std::println("For y({})", i);
    std::println("|{:^7}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|{:^10}|", "i",
                 "A(i)", "A(i-1)-A(i)", "alpha^k", "Rich error", "Order est.",
                 "f comps");
    // Table header separator
    std::println("|{:-^7}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^10}|", "",
                 "", "", "", "", "", "");

    // Table body
    int its = starting_steps;
    for (size_t j = 0; j < A_k.size(); j++) {
      if (j == 0) {
        std::println(
            "|{:^7}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
            its, A_k[j][i], "", "", "", "", f_comps.at(j));
      } else if (j == 1) {
        std::println(
            "|{:^7}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
            its, A_k[j][i], A_diff_k[j][i], "", rich_error[j][i], "",
            f_comps.at(j));
      } else {
        std::println(
            "|{:^7}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
            its, A_k[j][i], A_diff_k[j][i], alpha_k_computed[j][i],
            rich_error[j][i], order_estimate[j][i], f_comps[j]);
      }
      its *= 2;
    }
  }
}

#endif // !ODE_H
