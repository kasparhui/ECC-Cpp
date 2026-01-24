#include "Fp.h"
Fp::Fp(const int v) : value(v % 47) {
    if (value < 0) {
        value += 47;
    }
}

Fp::Fp(const Fp& other) : value(other.value) {}

Fp Fp::operator+(const Fp& other) const {
    return Fp(value + other.value);
}

Fp Fp::operator+(const int other) const {
    return *this + Fp(other);
}

Fp Fp::operator-(const Fp& other) const {
    return Fp(value - other.value);
}

Fp Fp::operator-(const int other) const {
    return *this - Fp(other);
}

Fp Fp::operator*(const Fp& other) const {
    return Fp(value * other.value);
}

Fp Fp::operator*(const int other) const {
    return *this * Fp(other);
}

Fp Fp::operator/(const Fp& other) const {
    // Using Fermat's Little Theorem for modular inverse
    int inv = 1;
    int base = other.value;
    int exp = 45; // 47 - 2
    while (exp) {
        if (exp % 2 == 1) {
            inv = (inv * base) % 47;
        }
        base = (base * base) % 47;
        exp /= 2;
    }
    return Fp(value * inv);
}

Fp Fp::operator/(const int other) const {
    return *this / Fp(other);
}