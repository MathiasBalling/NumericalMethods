#ifndef INTEGRALS_H
#define INTEGRALS_H

#include "utils/errors.h"
#include <cmath>
#include <cstddef>
#include <print>
#include <vector>
enum class IntegrationType { Midpoint, Trapezoidal, Simpson };

double midpoint(double (*func)(double), double limit_low, double limit_high,
                int its = 100) {
  double step = (limit_high - limit_low) / its;
  double result = 0.0;
  for (int i = 0; i < its; ++i) {
    double x1 = limit_low + i * step;
    double x2 = x1 + step;
    result += func((x1 + x2) / 2) * step;
  }
  return result;
}

double trapezoidal(double (*func)(double), double limit_low, double limit_high,
                   int its = 100) {
  double step = (limit_high - limit_low) / its;
  double result = 0.5 * func(limit_low);
  for (int i = 1; i < its; ++i) {
    double x = limit_low + i * step;
    result += func(x);
  }
  result += 0.5 * func(limit_high);

  result *= step;
  return result;
}

double simpson(double (*func)(double), double limit_low, double limit_high,
               int its = 100) {
  if (its % 2 != 0) {
    ++its; // Ensure even number of intervals
  }
  double step = (limit_high - limit_low) / its;
  double result = func(limit_low) + func(limit_high);
  for (int i = 1; i < its; ++i) {
    double x = limit_low + i * step;
    result += func(x) * (i % 2 == 0 ? 2 : 4);
  }
  result *= step / 3.0;
  return result;
}

void print_quadrature_table(double (*func)(double), double limit_low,
                            double limit_high, IntegrationType type,
                            double accuracy = 1e-10, int max_its = 10000) {
  std::vector<double> A_k;
  std::vector<int> f_comps;
  double expected_order = 0.0;
  double error = std::numeric_limits<double>::max();
  int its = 2;
  while (abs(error) > accuracy && its < max_its) {
    switch (type) {
    case IntegrationType::Midpoint:
      A_k.push_back(midpoint(func, limit_low, limit_high, its));
      f_comps.push_back(its);
      expected_order = 2.0;
      break;
    case IntegrationType::Trapezoidal:
      A_k.push_back(trapezoidal(func, limit_low, limit_high, its));
      f_comps.push_back(its + 1);
      expected_order = 2.0;
      break;
    case IntegrationType::Simpson:
      A_k.push_back(simpson(func, limit_low, limit_high, its));
      f_comps.push_back((its % 2 == 0 ? its : its + 1) + 1);
      expected_order = 4.0;
      break;
    default:
      throw("Unknown integration type");
    }
    error = richardson_extrapolation_error_current(A_k, pow(2, expected_order));
    its *= 2;
  }

  // Calculate the differences
  std::vector<double> A_diff_k = {NAN};
  for (size_t i = 1; i < A_k.size(); i++) {
    double diff = A_k[i - 1] - A_k[i];
    A_diff_k.push_back(diff);
  }

  // Calculate the orders
  auto alpha_k_computed = alpha_k_order_computed(A_k);

  // Calculate the richardson extrapolation
  auto rich_error = richardson_extrapolation_error(
      A_k,
      pow(2, expected_order)); // pow(2, expected_order) as we use N-1=1,2,4,8

  auto order_estimate = compute_order_estimate(alpha_k_computed);

  // Print the table
  // Table header
  std::println("|{:^6}|{:^21}|{:^21}|{:^21}|{:^21}|{:^21}|{:^10}|", "its",
               "A(i)", "A(i-1)-A(i)", "alpha^k", "Rich error", "Order est.",
               "f comps");
  // Table header separator
  std::println("|{:-^6}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|{:-^10}|", "",
               "", "", "", "", "", "");

  // Table body
  for (size_t i = 0; i < A_k.size(); i++) {
    int its = pow(2, i + 1);
    if (i == 0) {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
          its, A_k.at(i), "", "", "", "", f_comps.at(i));
    } else if (i == 1) {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
          its, A_k.at(i), A_diff_k.at(i), "", rich_error.at(i), "",
          f_comps.at(i));
    } else {
      std::println(
          "|{:^6}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21.12}|{:^10}|",
          its, A_k.at(i), A_diff_k.at(i), alpha_k_computed.at(i),
          rich_error.at(i), order_estimate.at(i), f_comps.at(i));
    }
  }
}
#endif // !INTEGRALS_H
