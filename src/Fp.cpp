#include "Fp.h"
// Precomputed modular inverses for field Fp where p = 47
const std::vector<int> inv_list = {
    0, 1, 24, 16, 12, 19, 8, 27, 6, 21, 
    33, 30, 4, 29, 37, 22, 3, 36, 34, 5, 
    40, 9, 15, 45, 2, 32, 38, 7, 42, 13, 
    11, 44, 25, 10, 18, 43, 17, 14, 26, 
    41, 20, 39, 28, 35, 31, 23, 46
};

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

Fp Fp::operator-() const {
    return Fp(-value);
}

Fp Fp::operator*(const Fp& other) const {
    return Fp(value * other.value);
}

Fp Fp::operator*(const int other) const {
    return *this * Fp(other);
}

Fp Fp::operator/(const Fp& other) const {
    int inv = inv_list[other.value];
    if (inv == 0) {
        throw std::runtime_error("Division by zero in Fp");
    }
    return Fp(value * inv);
}

Fp Fp::operator/(const int other) const {
    return *this / Fp(other);
}

std::ostream& operator<<(std::ostream& os, const Fp& fp) {
    os << fp.value;
    return os;
}