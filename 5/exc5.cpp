#include "roots.h"
#include "utils/roots.h"
#include <cstdlib>
#include <numbers>
#include <print>

// Only function
double f(double x) { return x - cos(x); }

// Function with derivative
struct fd {
  // Function
  Doub operator()(const Doub x) { return x - cos(x); }
  // Derivative of function
  Doub df(const Doub x) { return 1 + sin(x); }
};

int main() {
  // Bisection
  std::println("\n{:/^100}", " Bisection ");
  auto res_bisection =
      roots_table::root_bisection(f, 0, std::numbers::pi / 2, pow(10, -8));
  // From roots.h:
  std::println("Results:\nMine={}\nroots.h={}", res_bisection,
               rtbis(f, 0, std::numbers::pi / 2, pow(10, -8)));

  // False Position
  std::println("\n{:/^100}", " False Positive ");
  auto res_false_position = roots_table::root_false_position(
      f, 0, std::numbers::pi / 2, pow(10, -16));
  // From roots.h:
  std::println("Results:\nMine={}\nroots.h={}", res_false_position,
               rtflsp(f, 0, std::numbers::pi / 2, pow(10, -16)));

  // Secant
  std::println("\n{:/^100}", " Secant ");
  auto res_secant =
      roots_table::root_secant(f, 0, std::numbers::pi / 2, pow(10, -16));
  // From roots.h:
  std::println("Results:\nMine={}\nroots.h={}", res_secant,
               rtsec(f, 0, std::numbers::pi / 2, pow(10, -16)));

  // Newton
  std::println("\n{:/^100}", " Newton ");
  auto funcd = fd(); // Struct with function and derivative
  auto res_newton =
      roots_table::root_newton(funcd, 0, std::numbers::pi / 2, pow(10, -16));
  // From roots.h:
  std::println("Results:\nMine={}\nroots.h={}", res_newton,
               rtnewt(funcd, 0, std::numbers::pi / 2, pow(10, -16)));

  // Ridder
  std::println("\n{:/^100}", " Ridder ");
  auto res_ridder =
      roots_table::root_ridder(f, 0, std::numbers::pi / 2, pow(10, -16));
  // From roots.h:
  std::println("Results:\nMine={}\nroots.h={}", res_ridder,
               zriddr(f, 0, std::numbers::pi / 2, pow(10, -16)));
  return EXIT_SUCCESS;
}
