#pragma once

#include <vector>

#include "Poly.h"

// Division polynomials of the curve y^2 = x^3 + ax + b.
//
// The n-th division polynomial psi_n vanishes exactly at the points P with
// n * P = infinity. For odd n it is a polynomial in x alone; for even n it is
// y times a polynomial in x. This class stores the x-part f_n:
//   psi_n = f_n      for odd n
//   psi_n = y * f_n  for even n
class DivisionPolynomials {
public:
    // Computes f_0 .. f_max_index (max_index is raised to at least 4)
    DivisionPolynomials(const Fq& curve_a, const Fq& curve_b, int max_index);

    // Largest index that is stored
    int max_index() const;

    // f_n for 0 <= n <= max_index()
    const Poly& stored(int n) const;

    // f_n for odd n. Works for stored n and, without storing the result, for
    // every odd n up to 2 * max_index() - 3.
    Poly odd(int n) const;

private:
    // f_(2m+1) from f_(m-1) .. f_(m+2)
    Poly odd_step(int m) const;
    // f_(2m) from f_(m-2) .. f_(m+2)
    Poly even_step(int m) const;

    Poly curve_squared_; // (x^3 + ax + b)^2
    std::vector<Poly> table_;
};
