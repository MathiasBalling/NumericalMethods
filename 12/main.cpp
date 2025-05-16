#include "utils/bvp.h"

// Problem:
// y_mm(x)  = 2*x+sin(y'(x))-cos(y(x)) for 0<x<2
// y(0)=0
// y(2)=1

// y''=F
double F(double y_prime, double y, double x) {
  return 2 * x + sin(y_prime) - cos(y);
};

// d/dy F = F_y
double F_y(double y_prime, double y, double x) {
  (void)y_prime; // Unused
  (void)x;       // Unused
  return sin(y);
};

// d/d(y') F = F_y_prime
double F_y_prime(double y_prime, double y, double x) {
  (void)y; // Unused
  (void)x; // Unused
  return cos(y_prime);
};

int main() {
  {
    double a = 0;
    double b = 2;
    double alpha = 0;
    double beta = 1;
    int N = 4;
    // y(1) = ?
    double target_x = 1;
    finite_difference_method_table(N, target_x, a, b, alpha, beta, F, F_y,
                                   F_y_prime, 0, 1000);
  }
}
