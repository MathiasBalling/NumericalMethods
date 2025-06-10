#include "nr3.h"
#include "svd.h"
#include "utils/linreg_error.h"
#include "utils/utilities.h"
#include <print>

int main() {
  // Load data from FilipData.dat
  std::ifstream A_data("../../../exam/NUM_S25_Ex1A.dat");
  if (!A_data) {
    std::println(stderr, "Error opening NUM_S25_Ex1A.dat");
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
  }

  {
    // 2)
    std::println("\n{:/^100}", " ii ");
    double threshold = 1e-14;
    SVD svd_solver(A);
    auto nullspace = svd_solver.nullspace(threshold);
    util::print(nullspace, "nullspace");
  }

  {
    // 3)
    std::println("\n{:/^100}", " iii ");

    double threshold = 1e-14;
    SVD svd_solver(A);
    VecDoub x_sol(A_width);
    svd_solver.solve(b, x_sol, threshold);
    util::print(x_sol, "x_sol");
  }

  {
    // 4)
    std::println("\n{:/^100}", " iv ");

    double threshold = 1e-14;
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
