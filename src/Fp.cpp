#include "Fp.h"

Fp::Fp() : value(0) {}

Fp::Fp(const int v) : value(v) {
}

Fp::Fp(const Fp& other) : value(other.value) {}

void Fp::mod() {
    value = value % p;
    if (value < 0) {
        value += p;
    }
}

Fp Fp::operator+(const Fp& other) const {
    return (value+other.value)>=p ? value + other.value - p : value + other.value;
}

Fp Fp::operator+(const int other) const {
    return (value+other)>=p ? value + other - p : value + other;
}

Fp Fp::operator-(const Fp& other) const {
    return (value-other.value)<0 ? value - other.value + p : value - other.value;
}

Fp Fp::operator-(const int other) const {
    return (value-other)<0 ? value - other + p : value - other;
}

Fp Fp::operator-() const {
    return (value != 0) ? p - value : 0;
}

Fp Fp::operator*(const Fp& other) const {
    return mul_table[value * p + other.value];
}

Fp Fp::operator*(const int other) const {
    return mul_table[value * p + other];
}

Fp Fp::operator/(const Fp& other) const {
    return div_table[value * p + other.value];
}

Fp Fp::operator/(const int other) const {
    return div_table[value * p + other];
}

bool Fp::operator==(const Fp& other) const {
    return value == other.value;
}

std::ostream& operator<<(std::ostream& os, const Fp& fp) {
    os << fp.value;
    return os;
}