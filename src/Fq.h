#pragma once

#include "Fp.h"

constexpr int d = 46; // Degree of the field extension

// An element of the finite field Fq with q = 47^46, represented as a
// polynomial of degree < 46 over Fp modulo the irreducible polynomial x^46 + 2
struct Fq {
    std::array<Fp, d> coeffs{}; // coeffs[i] is the coefficient of x^i

    Fq() = default;

    Fq(const std::array<Fp, d>& c);

    // The element 1
    static Fq one();

    bool is_zero() const;

    Fq operator+(const Fq& other) const;

    Fq operator-(const Fq& other) const;

    Fq operator-() const;

    Fq operator*(const Fq& other) const;

    Fq operator*(const int other) const;

    // exp is a 256-bit number, most significant 32-bit word first
    Fq pow(const std::array<unsigned int, 8>& exp) const;

    Fq pow(const unsigned int exp) const;

    // k-fold Frobenius map: computes (*this)^(47^k), for k >= 0
    Fq frobenius(int k) const;

    // Inverting zero is not defined; it returns zero
    Fq invert() const;

    Fq operator/(const Fq& other) const;

    bool operator==(const Fq& other) const;
};

std::ostream& operator<<(std::ostream& os, const Fq& fq);
