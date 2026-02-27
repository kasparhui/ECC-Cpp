#include "Fq.h"

Fq::Fq(): coeffs{}{}

Fq::Fq(const std::array<Fp, 46>& coeffs) : coeffs(coeffs) {
}

Fq::Fq(const Fq &other) : coeffs(other.coeffs) {}

Fq::Fq(Fq&& other) noexcept : coeffs(std::move(other.coeffs)) {}

Fq& Fq::operator=(Fq&& other) noexcept {
    if (this != &other) {
        coeffs = std::move(other.coeffs);
    }
    return *this;
}

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
    return result_coeffs;
}

Fq Fq::operator-(const Fq& other) const {
    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = coeffs[i] - other.coeffs[i];
    }
    return result_coeffs;
}

Fq Fq::operator-() const {
    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = -coeffs[i];
    }
    return result_coeffs;
}

Fq Fq::operator*(const Fq& other) const {

    std::array<Fp, 46> result_coeffs{};
    for (int i = 0; i < 46; ++i) {
        for (int j = 0; j < 46; ++j) {
            int idx = i + j;
            if (idx >= 46) {
                result_coeffs[(idx) - 46] = result_coeffs[(idx) - 46] - (coeffs[i] * other.coeffs[j]) * 2;
            } else {
                result_coeffs[idx] = result_coeffs[idx] + (coeffs[i] * other.coeffs[j]);
            }
        }
    }
    return result_coeffs;
}

Fq Fq::operator*(const int other) const {
    std::array<Fp, 46> result_coeffs;
    for (int i = 0; i < 46; ++i) {
        result_coeffs[i] = coeffs[i] * other;
    }
    return result_coeffs;
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
                result = std::move(result) * power;
            }
            power = std::move(power) * power;
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
    for (int i = 0; i < 46; ++i) {
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
