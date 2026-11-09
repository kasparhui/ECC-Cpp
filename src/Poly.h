#pragma once

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Fq.h"

// A polynomial with coefficients in Fq.
// coeffs[i] is the coefficient of x^i. The highest stored coefficient is
// never zero, so the zero polynomial has no coefficients at all.
struct Poly {
    std::vector<Fq> coeffs;

    Poly() = default;

    // Takes any coefficient list, lowest degree first; trailing zeros are removed
    explicit Poly(std::vector<Fq> c);

    static Poly constant(const Fq& value);

    static Poly constant(int value);

    // The polynomial x
    static Poly x();

    // Degree of the polynomial; -1 for the zero polynomial
    int degree() const;

    bool is_zero() const;

    // Coefficient of the highest power; the polynomial must not be zero
    const Fq& leading() const;

    // Removes zero coefficients at the top
    void trim();

    Poly operator+(const Poly& other) const;

    Poly operator-(const Poly& other) const;

    Poly operator-() const;

    Poly operator*(const Poly& other) const;

    Poly operator*(const Fq& factor) const;

    Poly operator*(int factor) const;

    bool operator==(const Poly& other) const;

    // The same polynomial divided by its leading coefficient
    Poly monic() const;
};

struct DivisionResult {
    Poly quotient;
    Poly remainder;
};

// Polynomial long division: numerator = quotient * denominator + remainder
// with degree(remainder) < degree(denominator). The denominator must not be zero.
DivisionResult divide(const Poly& numerator, const Poly& denominator);

// Greatest common divisor, with leading coefficient 1 (zero if both are zero)
Poly gcd(const Poly& first, const Poly& second);

// Product of two coefficient lists. out must have room for nx + ny - 1
// entries and must not overlap the inputs. Uses Karatsuba multiplication.
void multiply(const Fq* x, std::size_t nx, const Fq* y, std::size_t ny, Fq* out);

std::ostream& operator<<(std::ostream& os, const Poly& poly);
