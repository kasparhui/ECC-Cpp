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
    return (value>=other.value) ? value - other.value : p + value - other.value;
}

Fp Fp::operator-(const int other) const {
    return (value>=other) ? value - other : p + value - other;
}

Fp Fp::operator-() const {
    return (value != 0) ? p - value : 0;
}

Fp Fp::operator*(const Fp& other) const {
    return mul_table[value][other.value];
}

Fp Fp::operator*(const int other) const {
    return mul_table[value][other];
}

Fp Fp::operator/(const Fp& other) const {
    return div_table[value][other.value];
}

Fp Fp::operator/(const int other) const {
    return div_table[value][other];
}

bool Fp::operator==(const Fp& other) const {
    return value == other.value;
}

std::ostream& operator<<(std::ostream& os, const Fp& fp) {
    os << fp.value;
    return os;
}