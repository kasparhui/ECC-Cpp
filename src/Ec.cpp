#include "Ec.h"

#include <ostream>

const Fq a = Fq({
    Fp(45), Fp(38), Fp(38), Fp(37), Fp(38), Fp(12), Fp(25), Fp(27),
    Fp(45), Fp(34), Fp(29), Fp(4), Fp(27), Fp(7), Fp(36), Fp(18),
    Fp(0), Fp(19), Fp(23), Fp(46), Fp(10), Fp(16), Fp(7), Fp(39),
    Fp(6), Fp(22), Fp(7), Fp(20), Fp(2), Fp(26), Fp(12), Fp(37),
    Fp(38), Fp(29), Fp(36), Fp(37), Fp(10), Fp(15), Fp(15), Fp(10),
    Fp(7), Fp(41), Fp(25), Fp(22), Fp(17), Fp(3)
});
const Fq b = Fq({
    Fp(45), Fp(17), Fp(27), Fp(27), Fp(38), Fp(26), Fp(46), Fp(39),
    Fp(10), Fp(31), Fp(24), Fp(17), Fp(6), Fp(40), Fp(7),
    Fp(24), Fp(34), Fp(41), Fp(4), Fp(20), Fp(10), Fp(38), Fp(33),
    Fp(36), Fp(38), Fp(20), Fp(13), Fp(27), Fp(0), Fp(29), Fp(28),
    Fp(29), Fp(41), Fp(32), Fp(46), Fp(26), Fp(41), Fp(40), Fp(38),
    Fp(26), Fp(46), Fp(22), Fp(5), Fp(40), Fp(29), Fp(10)
});

Point::Point(const Fq& x_coord, const Fq& y_coord) : x(x_coord), y(y_coord), is_infinity(false) {}

Point Point::operator+(const Point& other) const {
    if (is_infinity) return other;
    if (other.is_infinity) return *this;

    Fq slope;
    if (x == other.x) {
        if (y == -other.y) {
            return Point(); // P + (-P) is the point at infinity
        }
        // Same point: doubling, slope of the tangent
        slope = (x * x * 3 + a) / (y * 2);
    } else {
        // Different points: slope of the line through both
        slope = (other.y - y) / (other.x - x);
    }

    Fq x3 = slope * slope - x - other.x;
    Fq y3 = slope * (x - x3) - y;
    return Point(x3, y3);
}

Point Point::operator-() const {
    if (is_infinity) return *this;

    // The negative of (x, y) is (x, -y)
    return Point(x, -y);
}

Point Point::operator-(const Point& other) const {
    return *this + (-other);
}

Point Point::operator*(int scalar) const {
    if (is_infinity) return *this;

    // double-and-add algorithm for scalar multiplication
    Point result; // Point at infinity
    Point addend = *this;

    while (scalar > 0) {
        if (scalar & 1) {
            result = result + addend;
        }
        addend = addend + addend;
        scalar >>= 1;
    }

    return result;
}

Point Point::operator*(const std::array<unsigned int, 8>& scalar) const {
    if (is_infinity) return *this;

    // Windowed NAF (width 4) scalar multiplication.
    // The scalar is rewritten with digits 0, +-1, +-3, +-5, +-7 such that any
    // non-zero digit is followed by at least three zeros. Only about every
    // 5th digit is non-zero, so far fewer additions are needed than with
    // plain double-and-add. Negative digits are subtractions, which are as
    // cheap as additions because negating a point is free.

    // Copy the scalar, least significant word first (scalar[7] is the lowest
    // word). One extra word catches a possible carry out of the top.
    std::array<unsigned int, 9> k{};
    for (int i = 0; i < 8; ++i) {
        k[i] = scalar[7 - i];
    }

    std::array<int8_t, 258> digits{};
    int digit_count = 0;
    auto is_zero = [&k]() {
        for (unsigned int word : k) {
            if (word != 0) return false;
        }
        return true;
    };
    while (!is_zero()) {
        int digit = 0;
        if (k[0] & 1) {
            digit = k[0] & 15;              // k mod 16
            if (digit >= 8) digit -= 16;    // choose the digit in -7..7
            // k -= digit, which makes k divisible by 16
            if (digit > 0) {
                k[0] -= digit;              // low bits are >= digit: no borrow
            } else {
                unsigned long long carry = static_cast<unsigned long long>(-digit);
                for (int i = 0; i < 9 && carry != 0; ++i) {
                    carry += k[i];
                    k[i] = static_cast<unsigned int>(carry);
                    carry >>= 32;
                }
            }
        }
        digits[digit_count++] = static_cast<int8_t>(digit);
        // k >>= 1
        for (int i = 0; i < 8; ++i) {
            k[i] = (k[i] >> 1) | (k[i + 1] << 31);
        }
        k[8] >>= 1;
    }

    // Odd multiples: multiples[n] = (2n + 1) * P, i.e. P, 3P, 5P, 7P
    Point twice = *this + *this;
    std::array<Point, 4> multiples;
    multiples[0] = *this;
    for (int n = 1; n < 4; ++n) {
        multiples[n] = multiples[n - 1] + twice;
    }

    // Go from the most significant digit down: double, then add or subtract
    Point result; // Point at infinity
    for (int i = digit_count - 1; i >= 0; --i) {
        result = result + result;
        int digit = digits[i];
        if (digit > 0) {
            result = result + multiples[digit / 2];
        } else if (digit < 0) {
            result = result - multiples[(-digit) / 2];
        }
    }

    return result;
}

bool Point::operator==(const Point& other) const {
    if (is_infinity && other.is_infinity) return true;
    if (is_infinity || other.is_infinity) return false;
    return (x == other.x) && (y == other.y);
}

bool Point::is_on_curve() const {
    if (is_infinity) return true;

    // Check if the point satisfies the elliptic curve equation: y^2 = x^3 + ax + b
    Fq left = y * y;
    Fq right = (x * x * x) + (a * x) + b;
    return left == right;
}

std::ostream& operator<<(std::ostream& os, const Point& point) {
    if (point.is_infinity) {
        os << "Point at Infinity";
    } else {
        os << "Point(" << point.x << ", " << point.y << ")";
    }
    return os;
}