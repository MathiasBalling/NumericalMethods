#include "linreg_error.h"
#include "nr3.h"
#include "svd.h"
#include "utilities.h"
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
    svd_solver.solve(b, x_sol, svd_solver.eps);
    // Calculate fitted values
    // VecDoub y_fit(size);
    // for (int i = 0; i < size; i++) {
    //   y_fit[i] = 0.0;
    //   for (int j = 0; j < parameters; j++) {
    //     y_fit[i] += x_sol[j] * pow(x[i], j);
    //   }
    // }

    // Calculate error
    std::println("Usinge threshold: {}", threshold);
    auto sol_sigma =
        linreg_error::calculate_standard_deviation_svd(svd_solver, threshold);
    double residual_error = linreg_error::calculate_residual_error(A, b, x_sol);
    double r_fitting =
        linreg_error::calculate_random_fitting(A.nrows(), A.ncols());
    std::println("Residual error: {}", residual_error);
    std::println("Random fitting: {}", r_fitting);
    util::print(sol_sigma, "Std. dev.");

    ////////////////////////////////////////////////////////
    // With updated sigma
    auto r = (A * x_sol) - b;
    util::print(r, "r");
    std::println("\n{:/^100}", " with updated sigma ");
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
    auto svd_solver_sigma = SVD(A_sigma);
    VecDoub x_sol_sigma(parameters);
    svd_solver_sigma.solve(b_sigma, x_sol_sigma, svd_solver_sigma.eps);
    // Calculate fitted values
    // VecDoub y_fit_sigma(size);
    // for (int i = 0; i < size; i++) {
    //   y_fit_sigma[i] = 0.0;
    //   for (int j = 0; j < parameters; j++) {
    //     y_fit_sigma[i] += x_sol_sigma[j] * pow(x[i], j);
    //   }
    // }

    // Calculate error
    std::println("Usinge threshold: {}", threshold);
    auto sol_sigma_sigma = linreg_error::calculate_standard_deviation_svd(
        svd_solver_sigma, threshold);
    double residual_error_sigma =
        linreg_error::calculate_residual_error(A_sigma, b_sigma, x_sol_sigma);
    double r_fitting_sigma =
        linreg_error::calculate_random_fitting(A.nrows(), A.ncols());
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
  VecDoub xFilip(82);
  VecDoub yFilip(82);
  std::ifstream Filip("../../../3/FilipData.dat");
  if (!Filip) {
    std::println(stderr, "Error opening FilipData.dat");
    return 1;
  }
  for (int i = 0; i < 82; i++) {
    Filip >> yFilip[i];
    Filip >> xFilip[i];
  }
  // util::print(xFilip);
  // util::print(yFilip);

  // Load data from PontiusData.dat
  VecDoub xPont(40);
  VecDoub yPont(40);
  VecDoub sigmaPont(40);
  std::ifstream Pont("../../../3/PontiusData.dat");
  if (!Pont) {
    std::println(stderr, "Error opening PontiusData.dat");
    return 1;
  }
  for (int i = 0; i < 40; i++) {
    Pont >> yPont[i];
    Pont >> xPont[i];
  }
  // util::print(xPont);
  // util::print(yPont);

  // With sigma = 1
  std::println("\n{:/^100}", " Pont ");
  auto sol_pont = solve_SVD(xPont, yPont, 3, 0.1);
  util::print(sol_pont.first, "Pont");
  util::print(sol_pont.second, "Pont sigma");
  std::println("\n{:/^100}", " Filip ");
  auto sol_filip = solve_SVD(xFilip, yFilip, 11, 0.001);
  util::print(sol_filip.first, "Filip");
  util::print(sol_filip.second, "Filip sigma");

  return 0;
}
