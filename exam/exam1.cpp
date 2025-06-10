#include "nr3.h"
#include "svd.h"
#include "utils/linreg_error.h"
#include "utils/utilities.h"
#include <print>
#include <utility>

std::pair<VecDoub, VecDoub> solve_SVD(VecDoub x, VecDoub y, int parameters,
                                      double threshold) {
  int size = x.size();
  MatDoub A(size, parameters);
  VecDoub b(size);

  // NO sigma! (sigma=1)
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < parameters; j++) {
      A[i][j] = pow(x[i], j);
    }
    b[i] = y[i];
  }

  try {
    auto svd_solver = SVD(A);
    VecDoub x_sol(parameters);
    svd_solver.solve(b, x_sol, threshold);

    // Calculate error
    std::println("Using threshold: {}", threshold);
    auto sol_sigma =
        linreg_error::calculate_standard_deviation_svd(svd_solver, threshold);
    double residual_error = linreg_error::calculate_residual_error(A, b, x_sol);
    double r_fitting =
        linreg_error::calculate_random_fitting(A.nrows(), A.ncols());

    util::print(x_sol, "Parameters");
    std::println("Residual error: {}", residual_error);
    std::println("Random fitting: {}", r_fitting);
    util::print(sol_sigma, "Std. dev.");

    ////////////////// With updated sigma ///////////////////////////
    std::println("\n{:-^100}", " with updated sigma ");
    auto r = (A * x_sol) - b;

    double sigma;
    MatDoub A_sigma = A;
    VecDoub b_sigma = b;
    for (int i = 0; i < A.nrows(); i++) {
      // sigma = std::max(
      //     1.0,
      //     abs(r[i]));
      // FIX: Is it correct to use max(..,..) or residual_error?
      sigma = residual_error;
      for (int j = 0; j < A.ncols(); j++) {
        A_sigma[i][j] = A_sigma[i][j] / sigma;
      }
      b_sigma[i] = b_sigma[i] / sigma;
    }

    auto svd_solver_sigma = SVD(A_sigma);
    VecDoub x_sol_sigma(parameters);
    svd_solver_sigma.solve(b_sigma, x_sol_sigma, threshold);

    // Calculate error
    std::println("Using threshold: {}", threshold);

    auto sol_sigma_sigma = linreg_error::calculate_standard_deviation_svd(
        svd_solver_sigma, threshold);
    double residual_error_sigma =
        linreg_error::calculate_residual_error(A_sigma, b_sigma, x_sol_sigma);
    double r_fitting_sigma =
        linreg_error::calculate_random_fitting(A.nrows(), A.ncols());

    util::print(x_sol_sigma, "Parameters with sigma");
    std::println("Residual error sigma: {}", residual_error_sigma);
    std::println("Random fitting sigma: {}", r_fitting_sigma);
    util::print(sol_sigma_sigma, "Std. dev. with sigma");

    return std::make_pair(x_sol, x_sol_sigma);
  } catch (...) {
    std::println("Error in SVD");
    return std::make_pair(VecDoub(0), VecDoub(0));
  }
}

int main() {
  // Load data from FilipData.dat
  std::ifstream A_data("../../../exam/NUM_S25_Ex1A.dat");
  if (!A_data) {
    std::println(stderr, "Error opening FilipData.dat");
    return 1;
  }
  int A_width, A_height;
  A_data >> A_height;
  A_data >> A_width;

  MatDoub A(A_height, A_width);
  for (int i = 0; i < A_height; i++) {
    for (int j = 0; j < A_width; j++) {
      A_data >> A[i][j];
    }
  }

  std::ifstream b_data("../../../exam/NUM_S25_Ex1b.dat");
  if (!b_data) {
    std::println(stderr, "Error opening FilipData.dat");
    return 1;
  }
  int b_width, b_height;
  b_data >> b_height;
  b_data >> b_width;

  VecDoub b(b_height);
  for (int i = 0; i < b_height; i++) {
    b_data >> b[i];
  }

  // util::print(A, "A");
  // util::print(b, "b");

  {
    // 1)
    std::println("\n{:/^100}", " i ");
    SVD svd_solver(A);
    util::print(svd_solver.w, "w");
    // w , vec_size 8
    // 11.951932409770585      5.397266515466196       4.719596639028403       3.953906595603368
    // 3.6230909930002446       3.368300112475192       3.0916664682435697      4.811863543913732e-15
  }

  {
    // 2)
    std::println("\n{:/^100}", " ii ");
    // TODO: Null space
  }

  {
    // 3)
    std::println("\n{:/^100}", " iii ");

    double threshold = 0.0000000001;
    SVD svd_solver(A);
    VecDoub x_sol(A_width);
    svd_solver.solve(b, x_sol, threshold);
    util::print(x_sol, "x_sol");
  }

  {
    // 4)
    std::println("\n{:/^100}", " iv ");

    double threshold = 0.0000000001;
    SVD svd_solver(A);
    VecDoub x_sol(A_width);
    svd_solver.solve(b, x_sol, -1);

    auto sol_sigma =
        linreg_error::calculate_standard_deviation_svd(svd_solver, threshold);
    double residual_error = linreg_error::calculate_residual_error(A, b, x_sol);
    double r_fitting =
        linreg_error::calculate_random_fitting(A.nrows(), A.ncols());

    util::print(x_sol, "Parameters");
    std::println("Residual error: {}", residual_error);
    std::println("Random fitting: {}", r_fitting);
    util::print(sol_sigma, "Std. dev.");
  }
  return 0;
}
