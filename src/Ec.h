#pragma once

#include "Fq.h"

extern const Fq a; // Coefficient 'a' of the elliptic curve equation
extern const Fq b; // Coefficient 'b' of the elliptic curve equation

// A point on the elliptic curve y^2 = x^3 + ax + b over Fq
struct Point {
    Fq x;
    Fq y;
    bool is_infinity = true;

    // The point at infinity
    Point() = default;

    Point(const Fq& x_coord, const Fq& y_coord);

    Point operator+(const Point& other) const;

    Point operator-() const;

    Point operator-(const Point& other) const;

    // scalar must not be negative
    Point operator*(int scalar) const;

    // scalar is a 256-bit number, most significant 32-bit word first
    Point operator*(const std::array<unsigned int, 8>& scalar) const;

    bool operator==(const Point& other) const;

    bool is_on_curve() const;
};

std::ostream& operator<<(std::ostream& os, const Point& point);
