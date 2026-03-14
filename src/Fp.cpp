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
    return add_table[value][other.value];
}

Fp Fp::operator+(const int other) const {
    return add_table[value][other];
}

Fp Fp::operator-(const Fp& other) const {
    return sub_table[value][other.value];
}

Fp Fp::operator-(const int other) const {
    return sub_table[value][other];
}

Fp Fp::operator-() const {
    return sub_table[0][value];
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
    os << int(fp.value);
    return os;
}