#include "Fq.h"

const int d = 46; // Degree of the field extension

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
// unspecified division operator
Fq Fq::operator/(const Fq& other) const {
    return *this;
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
