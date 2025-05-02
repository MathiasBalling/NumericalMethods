#include "linreg_error.h"
#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <print>

int main() {
  // Load data from FilipData.dat
  std::ifstream A_file("../../../mandatory_1/Ex1A.dat");
  std::ifstream b_file("../../../mandatory_1/Ex1b.dat");
  if (!A_file) {
    std::println(stderr, "Error opening Ex1A.dat");
    return 1;
  }
  if (!b_file) {
    std::println(stderr, "Error opening Ex1b.dat");
    return 1;
  }
  int A_rows, A_cols, b_rows, b_cols;
  A_file >> A_rows >> A_cols;
  b_file >> b_rows >> b_cols;
  assert(A_rows == b_rows);
  MatDoub A(A_rows, A_cols);
  VecDoub b(b_rows);
  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      A_file >> A[i][j];
    }
    b_file >> b[i];
  }
  // util::print(A, "A");
  // util::print(b, "b");

  ////////////////////// i //////////////////////
  std::println("\n{:/^100}", " i ");
  auto svd_solver = SVD(A);
  auto w = svd_solver.w;
  util::print(w, "w");

  ///////////////////// ii //////////////////////
  std::println("\n{:/^100}", " ii ");
  auto x = VecDoub(A_cols);
  double threshold = 0.0000000001;
  svd_solver.solve(b, x, threshold);
  util::print(x, "x");

  ///////////////////// iii //////////////////////
  std::println("\n{:/^100}", " iii ");
  // Calculate fitted values
  VecDoub b_fit = A * x;
  std::println("Using threshold: {}", threshold);
  double residual_error = linreg_error::calculate_residual_error(A, b, x);
  double r_fitting =
      linreg_error::calculate_random_fitting(A.nrows(), A.ncols());
  auto sigma_parameters =
      linreg_error::calculate_standard_deviation_svd(svd_solver, threshold);
  std::println("Residual error: {}", residual_error);
  std::println("Random fitting: {}", r_fitting);
  util::print(sigma_parameters, "Std. dev.");

  ///////////////////// iv //////////////////////
  std::println("\n{:/^100}", " iv ");
  auto r = A * x - b;
  util::print(r, "r");

  ////////////////////// v //////////////////////
  std::println("\n{:/^100}", " v ");
  double sigma;
  MatDoub A_sigma = A;
  VecDoub b_sigma = b;
  for (int i = 0; i < A.nrows(); i++) {
    sigma = std::max(1.0, abs(r[i]));
    for (int j = 0; j < A.ncols(); j++) {
      A_sigma[i][j] = A_sigma[i][j] / sigma;
    }
    b_sigma[i] = b_sigma[i] / sigma;
  }

  std::println("A[0][0]={}", A_sigma[0][0]);
  std::println("b[6]={}", b_sigma[6]);

  ///////////////////// vi //////////////////////
  std::println("\n{:/^100}", " vi ");
  auto svd_solver_sigma = SVD(A_sigma);
  VecDoub x_sigma = VecDoub(A_sigma.ncols());
  svd_solver_sigma.solve(b_sigma, x_sigma, threshold);
  util::print(x_sigma, "x_sigma");

  // Error estimation for new (Not needed)
  // VecDoub b_fit_sigma = A_sigma * x_sigma;
  // std::println("Using threshold: {}", threshold);
  // double residual_error_sigma =
  //     linreg_error::calculate_residual_error(A_sigma, b_sigma, x_sigma);
  // double r_fitting_sigma =
  //     linreg_error::calculate_random_fitting(A.nrows(), A.ncols());
  // auto sigma_parameters_sigma =
  // linreg_error::calculate_standard_deviation_svd(
  //     svd_solver_sigma, threshold);
  // std::println("Residual error: {}", residual_error_sigma);
  // std::println("Random fitting: {}", r_fitting_sigma);
  // util::print(sigma_parameters_sigma, "Std. dev.");

  return 0;
}
