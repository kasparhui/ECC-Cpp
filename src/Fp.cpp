#include "Fp.h"

Fp::Fp() : value(0) {}

Fp::Fp(const int v) : value(v) {
    mod();
}

Fp::Fp(const Fp& other) : value(other.value) {}

void Fp::mod() {
    value = value % p;
    if (value < 0) {
        value += p;
    }
}

Fp Fp::operator+(const Fp& other) const {
    return value + other.value;
}

Fp Fp::operator+(const int other) const {
    return this->value + other;
}

Fp Fp::operator-(const Fp& other) const {
    return value - other.value;
}

Fp Fp::operator-(const int other) const {
    return this->value - other;
}

Fp Fp::operator-() const {
    return -this->value;
}

Fp Fp::operator*(const Fp& other) const {
    return value * other.value;
}

Fp Fp::operator*(const int other) const {
    return this->value * other;
}

Fp Fp::operator/(const Fp& other) const {
    return value * inv_list[other.value];
}

Fp Fp::operator/(const int other) const {
    return this->value * inv_list[other];
}

bool Fp::operator==(const Fp& other) const {
    return value == other.value;
}

std::ostream& operator<<(std::ostream& os, const Fp& fp) {
    os << fp.value;
    return os;
}