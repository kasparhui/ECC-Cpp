#include "Fq.h"

const int d = 46; // Degree of the field extension
const std::array<unsigned int, 8> q_sub_2 = {0xB66A744B, 0xFCDBBC9B, 0x96BBF58B, 0x7AE788BB,
                                            0xFF979CA4, 0x4C7FA532, 0xAEEBC56D, 0xE1751A5F};

Fq::Fq(): coeffs{}{}

Fq::Fq(const std::array<Fp, 46>& coeffs) : coeffs(coeffs) {
}

Fq::Fq(const Fq &other) : coeffs(other.coeffs) {}

void Fq::mod() {
    for (int i = 0; i < 46; ++i) {
        coeffs[i].mod();
    }
}

Fq Fq::operator+(const Fq &other) const {
    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = coeffs[i] + other.coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator-(const Fq& other) const {
    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = coeffs[i] - other.coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator-() const {
    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = -coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator*(const Fq& other) const {

    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        for (int j = 0; j < 46; ++j) {
            if (i + j >= 46) {
                result_coeffs[(i + j) - 46] = result_coeffs[(i + j) - 46] - (coeffs[i] * other.coeffs[j]) * 2;
            } else {
                result_coeffs[i + j] = result_coeffs[i + j] + (coeffs[i] * other.coeffs[j]);
            }
        }
    }
    return Fq(result_coeffs);
}

Fq Fq::operator*(const int other) const {
    std::array<Fp, 46> result_coeffs;
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = coeffs[i] * other;
    }
    return Fq(result_coeffs);
}

Fq Fq::pow(const std::array<unsigned int, 8>& exp) const {
    assert(!exp.empty());
    Fq result;
    result.coeffs[0] = Fp(1); // Initialize result to 1
    Fq power = *this;
    // liste von hinten durchgehen
    for (int idx = 8 - 1; idx >= 0; --idx) {
        unsigned int e = exp[idx];
        for (int bit = 0; bit < 32; ++bit) {
            if (e & 1) {
                result = result * power;
            }
            power = power * power;
            e >>= 1;
        }
    }
    return result;
}

Fq Fq::pow(const unsigned int exp) const {
    Fq result;
    result.coeffs[0] = Fp(1); // Initialize result to 1
    Fq power = *this;
    unsigned int e = exp;
    while (e > 0) {
        if (e & 1) {
            result = result * power;
        }
        power = power * power;
        e >>= 1;
    }
    return result;
}

Fq Fq::invert() const {
    // Using Fermat's little theorem: a^(q-2) mod q
    return this->pow(q_sub_2);
}

// using Fermat's little theorem for inversion
Fq Fq::operator/(const Fq& other) const {
    return (*this) * other.invert();
}

bool Fq::operator==(const Fq& other) const {
    for (int i = 0; i < d; ++i) {
        if (!(coeffs[i] == other.coeffs[i])) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Fq& fq) {
    os << '{';
    for (int i = 0; i < d; ++i) {
        os << fq.coeffs[i];
        if (i < d - 1) {
            os << ", ";
        }
    }
    os << '}';
    return os;
}
