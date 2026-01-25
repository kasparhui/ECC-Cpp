#include "Fq.h"

const int d = 46; // Degree of the field extension
const std::vector<unsigned int> q_sub_2 = {0xB66A744B, 0xFCDBBC9B, 0x96BBF58B, 0x7AE788BB,
                                            0xFF979CA4, 0x4C7FA532, 0xAEEBC56D, 0xE1751A5F};

Fq::Fq(const std::vector<Fp>& coeffs) : coeffs(coeffs) {
    assert(coeffs.size() == d);
}

Fq::Fq(const Fq &other) : coeffs(other.coeffs) {}

Fq Fq::operator+(const Fq &other) const {
    std::vector<Fp> result_coeffs(d);
    for (int i = 0; i < d; ++i) {
        result_coeffs[i] = coeffs[i] + other.coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator-(const Fq& other) const {
    std::vector<Fp> result_coeffs(d);
    for (int i = 0; i < d; ++i) {
        result_coeffs[i] = coeffs[i] - other.coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator-() const {
    std::vector<Fp> result_coeffs(d);
    for (int i = 0; i < d; ++i) {
        result_coeffs[i] = -coeffs[i];
    }
    return Fq(result_coeffs);
}

Fq Fq::operator*(const Fq& other) const {
    std::vector<Fp> result_coeffs(2*d-1, Fp(0)); // Degree can go up to 2(d-1)
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            result_coeffs[i + j] = result_coeffs[i + j] + (coeffs[i] * other.coeffs[j]);
        }
    }
    // Reduce modulo x^d + 2
    for (int k = 0; k < d; ++k) {
        result_coeffs[k] = result_coeffs[k] - result_coeffs[k + d] * 2;
    }
    result_coeffs.resize(d);
    return Fq(result_coeffs);
}

Fq Fq::operator*(const int other) const {
    std::vector<Fp> result_coeffs(d);
    for (int i = 0; i < d; ++i) {
        result_coeffs[i] = coeffs[i] * other;
    }
    return Fq(result_coeffs);
}

Fq Fq::pow(const std::vector<unsigned int>& exp) const {
    Fq result(std::vector<Fp>(d, Fp(0)));
    result.coeffs[0] = Fp(1); // Initialize result to 1
    Fq power = *this;
    // liste von hinten durchgehen
    for (int idx = exp.size() - 1; idx >= 0; --idx) {
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

Fq Fq::invert() const {
    // Using Fermat's little theorem: a^(q-2) mod q
    return this->pow(q_sub_2);
}

// using Fermat's little theorem for inversion
Fq Fq::operator/(const Fq& other) const {
    return (*this) * other.invert();
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
