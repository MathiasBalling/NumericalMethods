#ifndef LINREG_ERROR_H
#define LINREG_ERROR_H
#include "math.h"
#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <print>
namespace linreg_error {
double calculate_random_fitting(double m, double n) {
  // row=m
  // col=n
  return sqrt((m - n) / m);
}

double calculate_residual_error(MatDoub A, VecDoub b, VecDoub x) {
  VecDoub top = A * x - b;
  return norm(top) / norm(b);
}

VecDoub calculate_standard_deviation_svd(SVD &svd_solver, double threshold) {
  auto w = svd_solver.w;
  auto V = svd_solver.v;
  auto sigma = VecDoub(w.size());
  double sum;
  for (int j = 0; j < w.size(); j++) {
    sum = 0;
    for (int i = 0; i < w.size(); i++) {
      if (w[i] > threshold) {
        sum += pow(V[j][i] / w[i], 2);
      }
    }
    sigma[j] = sqrt(sum);
  }
  return sigma;
}

double calculate_r2(const VecDoub &true_val, const VecDoub &pred_val) {
  double mean_true_val = 0.0;
  for (int i = 0; i < true_val.size(); i++) {
    mean_true_val += true_val[i];
  }
  mean_true_val /= true_val.size();

  double ss_tot = 0.0;
  double ss_res = 0.0;
  for (int i = 0; i < true_val.size(); i++) {
    ss_tot += pow((true_val[i] - mean_true_val), 2.0);
    ss_res += pow((true_val[i] - pred_val[i]), 2.0);
  }

  return 1 - (ss_res / ss_tot);
}

double calculate_rmse(const VecDoub &true_val, const VecDoub &pred_val) {
  double sum = 0.0;
  int n = true_val.size();
  for (int i = 0; i < n; i++) {
    double diff = true_val[i] - pred_val[i];
    sum += pow(diff, 2.0);
  }
  return sqrt(sum / n);
}
} // namespace linreg_error
#endif // LINREG_ERROR_H
