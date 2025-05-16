#ifndef ERRORS_MULTI_H
#define ERRORS_MULTI_H
#include "nr3.h"

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

#endif // ERRORS_MULTI_H
