#ifndef ERRORS_H
#define ERRORS_H
#include "vector"
#include <cmath>

std::vector<double> alpha_k_order_computed(const std::vector<double> &A_k) {
  std::vector<double> alpha_k = {NAN, NAN}; // No error on the first two
  for (size_t i = 2; i < A_k.size(); i++) {
    const double A_1 = A_k[i - 2];
    const double A_2 = A_k[i - 1];
    const double A_3 = A_k[i];
    const double alpha_i = (A_1 - A_2) / (A_2 - A_3);
    alpha_k.push_back(alpha_i);
  }
  return alpha_k;
}

std::vector<double>
richardson_extrapolation_error(const std::vector<double> &A_k,
                               const double alpha_k_order_expected) {
  std::vector<double> A_R = {NAN}; // No error on the first
  for (size_t i = 1; i < A_k.size(); i++) {
    const double A_1 = A_k[i - 1];
    const double A_2 = A_k[i];
    double error = (A_2 - A_1) / (alpha_k_order_expected - 1);
    A_R.push_back(error);
  }
  return A_R;
}

double
richardson_extrapolation_error_current(const std::vector<double> &A_k,
                                       const double alpha_k_order_expected) {
  double error = std::numeric_limits<double>::max();
  if (A_k.size() < 2) {
    return error;
  } else {
    const double A_1 = A_k[A_k.size() - 2];
    const double A_2 = A_k[A_k.size() - 1];
    // INFO: Only use expected order if computed_order_estimate converges.
    error = (A_2 - A_1) / (alpha_k_order_expected - 1);
    return error;
  }
}

std::vector<double>
compute_order_estimate(const std::vector<double> &alpha_k_computed) {
  std::vector<double> order_estimate = {NAN,
                                        NAN}; // No order for first two entries
  for (size_t i = 2; i < alpha_k_computed.size(); i++) {
    double p = log2(abs(alpha_k_computed[i]));
    order_estimate.push_back(p);
  }
  return order_estimate;
}
#endif // !ERRORS_H
