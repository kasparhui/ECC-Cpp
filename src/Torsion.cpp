#include "Torsion.h"

#include <stdexcept>

TorsionGroup::TorsionGroup(const Fq& curve_a, const Fq& curve_b, const Modulus& modulus)
    : modulus_(modulus),
      curve_a_(Poly::constant(curve_a)),
      curve_(modulus.reduce(Poly({curve_b, curve_a, Fq(), Fq::one()}))) {}

Poly TorsionGroup::divide(const Poly& numerator, const Poly& denominator) const {
    InverseResult inverse = modulus_.inverse(denominator);
    if (!inverse.invertible) {
        if (inverse.value.degree() >= modulus_.degree()) {
            // The denominator is zero modulo the modulus; callers rule this out
            throw std::logic_error("TorsionGroup: division by zero");
        }
        throw FactorFound{inverse.value};
    }
    return modulus_.mul(numerator, inverse.value);
}

TorsionPoint TorsionGroup::from_slope(const Poly& slope, const TorsionPoint& first, const Poly& second_x) const {
    // The real slope is y * slope, and y^2 = x^3 + ax + b, so
    //   x3 = (y * slope)^2 - x1 - x2 = curve * slope^2 - x1 - x2
    //   y3 = y * slope * (x1 - x3) - y * y1 = y * (slope * (x1 - x3) - y1)
    TorsionPoint result;
    result.is_infinity = false;
    result.x = modulus_.mul(curve_, modulus_.mul(slope, slope)) - first.x - second_x;
    result.y = modulus_.mul(slope, first.x - result.x) - first.y;
    return result;
}

TorsionPoint TorsionGroup::double_point(const TorsionPoint& point) const {
    if (point.is_infinity || point.y.is_zero()) {
        return TorsionPoint();
    }
    // Real slope: (3x^2 + a) / (2 * y * y1). Divided by y this is
    //   (3x^2 + a) / (2 * y1 * curve)
    const Poly numerator = modulus_.mul(point.x, point.x) * 3 + curve_a_;
    const Poly denominator = modulus_.mul(point.y, curve_) * 2;
    return from_slope(divide(numerator, denominator), point, point.x);
}

TorsionPoint TorsionGroup::add(const TorsionPoint& first, const TorsionPoint& second) const {
    if (first.is_infinity) return second;
    if (second.is_infinity) return first;

    if (first.x == second.x) {
        if (first.y == second.y) {
            return double_point(first);
        }
        if (first.y == -second.y) {
            return TorsionPoint(); // P + (-P)
        }
        // Same x but the y-parts are neither equal nor opposite: this can
        // only happen if the modulus splits, and y1 - y2 exposes a factor.
        const InverseResult inverse = modulus_.inverse(first.y - second.y);
        if (inverse.invertible) {
            throw std::logic_error("TorsionGroup: inconsistent points");
        }
        throw FactorFound{inverse.value};
    }

    // Real slope: y * (y2 - y1) / (x2 - x1)
    const Poly slope = divide(second.y - first.y, second.x - first.x);
    return from_slope(slope, first, second.x);
}

TorsionPoint TorsionGroup::negate(const TorsionPoint& point) const {
    TorsionPoint result = point;
    result.y = -result.y;
    return result;
}

TorsionPoint TorsionGroup::multiply(const TorsionPoint& point, int factor) const {
    if (factor < 0) {
        throw std::invalid_argument("TorsionGroup::multiply: factor must not be negative");
    }
    // double-and-add
    TorsionPoint result;
    TorsionPoint addend = point;
    while (true) {
        if (factor & 1) {
            result = add(result, addend);
        }
        factor >>= 1;
        if (factor == 0) {
            break;
        }
        addend = add(addend, addend);
    }
    return result;
}

bool TorsionGroup::equal(const TorsionPoint& first, const TorsionPoint& second) const {
    if (first.is_infinity || second.is_infinity) {
        return first.is_infinity && second.is_infinity;
    }
    return first.x == second.x && first.y == second.y;
}
