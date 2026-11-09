#include "Fq.h"

#include <ostream>

Fq::Fq(const std::array<Fp, d>& c) : coeffs(c) {}

Fq Fq::one() {
    Fq result;
    result.coeffs[0] = Fp(1);
    return result;
}

bool Fq::is_zero() const {
    for (int i = 0; i < d; ++i) {
        if (coeffs[i].value != 0) {
            return false;
        }
    }
    return true;
}

Fq Fq::operator+(const Fq& other) const {
    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i] = coeffs[i] + other.coeffs[i];
    }
    return result;
}

Fq Fq::operator-(const Fq& other) const {
    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i] = coeffs[i] - other.coeffs[i];
    }
    return result;
}

Fq Fq::operator-() const {
    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i] = -coeffs[i];
    }
    return result;
}

Fq Fq::operator*(const Fq& other) const {
    // Because x^46 = -2, result coefficient i is
    //   sum over j of  a[j] * c[i - j],
    // where c[m] = b[m] for m >= 0 and c[m] = -2 * b[m + 46] for m < 0.
    // The sums are built with plain integers and reduced mod 47 once at the
    // end. All loops have a fixed length and work on small integer types, so
    // the compiler can vectorize them.

    // Sums are computed for 48 positions (a multiple of the SIMD width);
    // only the first 46 are used.
    constexpr int width = 48;

    // c shifted by d - 1 so indices start at 0: c_ext[m + 45] = c[m].
    // The zero padding at the end is read by the two unused positions.
    uint8_t c_ext[(d - 1) + width] = {};
    for (int m = 0; m < d - 1; ++m) {
        c_ext[m] = (other.coeffs[m + 1].value * (p - 2)) % p;
    }
    for (int m = 0; m < d; ++m) {
        c_ext[d - 1 + m] = other.coeffs[m].value;
    }

    // Every product is at most 46 * 46 = 2116, and 23 of them fit in a
    // uint16_t (23 * 2116 = 48668), so the 46 terms go into two accumulators.
    constexpr int half = d / 2;
    uint16_t first_half[width] = {};
    uint16_t second_half[width] = {};
    for (int j = 0; j < half; ++j) {
        const uint16_t a = coeffs[d - 1 - j].value;
        for (int i = 0; i < width; ++i) {
            first_half[i] += a * c_ext[i + j];
        }
    }
    for (int j = half; j < d; ++j) {
        const uint16_t a = coeffs[d - 1 - j].value;
        for (int i = 0; i < width; ++i) {
            second_half[i] += a * c_ext[i + j];
        }
    }

    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i].value = (uint32_t(first_half[i]) + second_half[i]) % p;
    }
    return result;
}

Fq Fq::operator*(const int other) const {
    const Fp factor(other);
    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i] = coeffs[i] * factor;
    }
    return result;
}

Fq Fq::pow(const std::array<unsigned int, 8>& exp) const {
    Fq result;
    result.coeffs[0] = Fp(1);
    Fq power = *this;
    // square-and-multiply, starting at the least significant word (the last)
    for (int idx = 7; idx >= 0; --idx) {
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
    result.coeffs[0] = Fp(1);
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

// Since x^46 = -2 and 47 = 46 + 1, x^47 = -2 * x. So the Frobenius map
// a -> a^47 only scales each coefficient: applied k times, coefficient i is
// multiplied by (-2)^(i*k mod 46). frobenius_factors[k][i] holds that factor.
constexpr std::array<std::array<uint8_t, d>, d> make_frobenius_factors() {
    std::array<uint8_t, d> neg2_pow{}; // neg2_pow[e] = (-2)^e mod 47
    int v = 1;
    for (int e = 0; e < d; ++e) {
        neg2_pow[e] = v;
        v = (v * (p - 2)) % p;
    }

    std::array<std::array<uint8_t, d>, d> table{};
    for (int k = 0; k < d; ++k) {
        for (int i = 0; i < d; ++i) {
            table[k][i] = neg2_pow[(i * k) % d];
        }
    }
    return table;
}

constexpr auto frobenius_factors = make_frobenius_factors();

Fq Fq::frobenius(int k) const {
    const std::array<uint8_t, d>& factors = frobenius_factors[k % d];
    Fq result;
    for (int i = 0; i < d; ++i) {
        result.coeffs[i].value = (coeffs[i].value * factors[i]) % p;
    }
    return result;
}

Fq Fq::invert() const {
    // Itoh-Tsujii inversion. With r = 1 + 47 + 47^2 + ... + 47^45 = (q-1)/(p-1),
    // a^r is the norm of a and lies in Fp, so a^-1 = a^(r-1) * (a^r)^-1.
    //
    // t_k = a^(1 + 47 + ... + 47^(k-1)), built with t_(m+n) = t_m * frobenius^m(t_n)
    const Fq& t1 = *this;
    const Fq t2 = t1 * t1.frobenius(1);
    const Fq t4 = t2 * t2.frobenius(2);
    const Fq t5 = t4 * t1.frobenius(4);
    const Fq t10 = t5 * t5.frobenius(5);
    const Fq t20 = t10 * t10.frobenius(10);
    const Fq t40 = t20 * t20.frobenius(20);
    const Fq t45 = t40 * t5.frobenius(40);

    // a^(r-1) = a^(47 + 47^2 + ... + 47^45)
    const Fq a_r_minus_1 = t45.frobenius(1);

    // norm = a * a^(r-1) is in Fp, so only its constant coefficient is needed:
    // a[0]*b[0] - 2 * (a[1]*b[45] + a[2]*b[44] + ... + a[45]*b[1])
    uint32_t cross = 0;
    for (int i = 1; i < d; ++i) {
        cross += coeffs[i].value * a_r_minus_1.coeffs[d - i].value;
    }
    const uint32_t norm = (coeffs[0].value * a_r_minus_1.coeffs[0].value + (p - 2) * cross) % p;

    return a_r_minus_1 * int(inv_list[norm]);
}

Fq Fq::operator/(const Fq& other) const {
    return (*this) * other.invert();
}

bool Fq::operator==(const Fq& other) const {
    return coeffs == other.coeffs;
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
