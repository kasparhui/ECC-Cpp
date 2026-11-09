#pragma once

#include "PolyMod.h"

// A point of the curve whose coordinates are polynomials in x, taken modulo
// a division polynomial (or a factor of one). It stands for all torsion
// points at once: plugging in the x-coordinate of a torsion point gives the
// coordinates of the result for that point.
//
// The y-coordinate always has the form y * (polynomial in x); only that
// polynomial is stored. So the point is (x, y * this->y).
struct TorsionPoint {
    Poly x;
    Poly y;
    bool is_infinity = true;
};

// Thrown when a division is not possible because the denominator shares a
// factor with the modulus. That factor is itself a valid, smaller modulus,
// so the computation can be restarted with it.
struct FactorFound {
    Poly factor;
};

// The group law on torsion points of y^2 = x^3 + ax + b modulo a polynomial
class TorsionGroup {
public:
    // The modulus must outlive the group
    TorsionGroup(const Fq& curve_a, const Fq& curve_b, const Modulus& modulus);

    // May throw FactorFound
    TorsionPoint add(const TorsionPoint& first, const TorsionPoint& second) const;

    TorsionPoint negate(const TorsionPoint& point) const;

    // factor * point for factor >= 0. May throw FactorFound.
    TorsionPoint multiply(const TorsionPoint& point, int factor) const;

    bool equal(const TorsionPoint& first, const TorsionPoint& second) const;

private:
    TorsionPoint double_point(const TorsionPoint& point) const;
    // numerator / denominator mod modulus; throws FactorFound if impossible
    Poly divide(const Poly& numerator, const Poly& denominator) const;
    // Shared end of addition and doubling: the result for a given slope / y
    TorsionPoint from_slope(const Poly& slope, const TorsionPoint& first, const Poly& second_x) const;

    const Modulus& modulus_;
    Poly curve_a_; // the constant a as a polynomial
    Poly curve_;   // x^3 + ax + b mod modulus
};
