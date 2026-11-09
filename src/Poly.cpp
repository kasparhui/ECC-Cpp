#include "Poly.h"

#include <algorithm>
#include <ostream>
#include <utility>

namespace {

// Up to this size the simple schoolbook product is used. A multiplication
// in Fq costs about 20 additions, so Karatsuba pays off almost immediately;
// 3 was the fastest value measured for degrees 84, 924 and 3960.
const std::size_t karatsuba_threshold = 3;

void schoolbook(const Fq* x, const std::size_t nx, const Fq* y, const std::size_t ny, Fq* out) {
    std::fill(out, out + nx + ny - 1, Fq());
    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) {
            out[i + j] = out[i + j] + x[i] * y[j];
        }
    }
}

// Scratch space karatsuba() needs for inputs of length n
std::size_t karatsuba_scratch_size(std::size_t n) {
    std::size_t total = 0;
    while (n > karatsuba_threshold) {
        const std::size_t high = n - n / 2;
        total += 4 * high;
        n = high;
    }
    return total + 1;
}

// Product of two lists of the same length n; out gets 2n - 1 entries.
// With x = x_low + x_high * X^low and the same for y:
//   x * y = z0 + (z1 - z0 - z2) * X^low + z2 * X^(2*low)
// where z0 = x_low * y_low, z2 = x_high * y_high and
// z1 = (x_low + x_high) * (y_low + y_high): three products instead of four.
void karatsuba(const Fq* x, const Fq* y, const std::size_t n, Fq* out, Fq* scratch) {
    if (n <= karatsuba_threshold) {
        schoolbook(x, n, y, n, out);
        return;
    }
    const std::size_t low = n / 2;
    const std::size_t high = n - low; // high >= low

    karatsuba(x, y, low, out, scratch);                    // z0 in out[0 .. 2*low - 2]
    out[2 * low - 1] = Fq();
    karatsuba(x + low, y + low, high, out + 2 * low, scratch); // z2 in out[2*low ..]

    Fq* x_sum = scratch;
    Fq* y_sum = scratch + high;
    Fq* z1 = scratch + 2 * high;
    Fq* rest = z1 + 2 * high - 1;
    for (std::size_t i = 0; i < low; ++i) {
        x_sum[i] = x[i] + x[low + i];
        y_sum[i] = y[i] + y[low + i];
    }
    if (high > low) {
        x_sum[low] = x[2 * low];
        y_sum[low] = y[2 * low];
    }
    karatsuba(x_sum, y_sum, high, z1, rest);

    for (std::size_t i = 0; i < 2 * low - 1; ++i) {
        z1[i] = z1[i] - out[i];
    }
    for (std::size_t i = 0; i < 2 * high - 1; ++i) {
        z1[i] = z1[i] - out[2 * low + i];
    }
    for (std::size_t i = 0; i < 2 * high - 1; ++i) {
        out[low + i] = out[low + i] + z1[i];
    }
}

} // namespace

void multiply(const Fq* x, std::size_t nx, const Fq* y, std::size_t ny, Fq* out) {
    if (nx < ny) {
        std::swap(x, y);
        std::swap(nx, ny);
    }
    if (ny == 0) {
        return;
    }
    if (ny <= karatsuba_threshold) {
        schoolbook(x, nx, y, ny, out);
        return;
    }
    if (nx == ny) {
        std::vector<Fq> scratch(karatsuba_scratch_size(nx));
        karatsuba(x, y, nx, out, scratch.data());
        return;
    }

    // Different lengths: cut the longer list into pieces as long as the shorter one
    std::fill(out, out + nx + ny - 1, Fq());
    std::vector<Fq> scratch(karatsuba_scratch_size(ny));
    std::vector<Fq> piece_product(2 * ny - 1);
    for (std::size_t offset = 0; offset < nx; offset += ny) {
        const std::size_t piece = std::min(ny, nx - offset);
        if (piece == ny) {
            karatsuba(x + offset, y, ny, piece_product.data(), scratch.data());
        } else {
            multiply(y, ny, x + offset, piece, piece_product.data());
        }
        for (std::size_t i = 0; i < piece + ny - 1; ++i) {
            out[offset + i] = out[offset + i] + piece_product[i];
        }
    }
}

Poly::Poly(std::vector<Fq> c) : coeffs(std::move(c)) {
    trim();
}

Poly Poly::constant(const Fq& value) {
    return Poly(std::vector<Fq>{value});
}

Poly Poly::constant(const int value) {
    Fq element;
    element.coeffs[0] = Fp(value);
    return constant(element);
}

Poly Poly::x() {
    return Poly(std::vector<Fq>{Fq(), Fq::one()});
}

int Poly::degree() const {
    return static_cast<int>(coeffs.size()) - 1;
}

bool Poly::is_zero() const {
    return coeffs.empty();
}

const Fq& Poly::leading() const {
    return coeffs.back();
}

void Poly::trim() {
    while (!coeffs.empty() && coeffs.back().is_zero()) {
        coeffs.pop_back();
    }
}

Poly Poly::operator+(const Poly& other) const {
    const Poly& longer = coeffs.size() >= other.coeffs.size() ? *this : other;
    const Poly& shorter = coeffs.size() >= other.coeffs.size() ? other : *this;
    Poly result = longer;
    for (std::size_t i = 0; i < shorter.coeffs.size(); ++i) {
        result.coeffs[i] = result.coeffs[i] + shorter.coeffs[i];
    }
    result.trim();
    return result;
}

Poly Poly::operator-(const Poly& other) const {
    Poly result = *this;
    if (result.coeffs.size() < other.coeffs.size()) {
        result.coeffs.resize(other.coeffs.size());
    }
    for (std::size_t i = 0; i < other.coeffs.size(); ++i) {
        result.coeffs[i] = result.coeffs[i] - other.coeffs[i];
    }
    result.trim();
    return result;
}

Poly Poly::operator-() const {
    Poly result = *this;
    for (Fq& coefficient : result.coeffs) {
        coefficient = -coefficient;
    }
    return result;
}

Poly Poly::operator*(const Poly& other) const {
    Poly result;
    if (is_zero() || other.is_zero()) {
        return result;
    }
    result.coeffs.resize(coeffs.size() + other.coeffs.size() - 1);
    multiply(coeffs.data(), coeffs.size(), other.coeffs.data(), other.coeffs.size(), result.coeffs.data());
    return result;
}

Poly Poly::operator*(const Fq& factor) const {
    Poly result = *this;
    for (Fq& coefficient : result.coeffs) {
        coefficient = coefficient * factor;
    }
    result.trim(); // only needed if factor is zero
    return result;
}

Poly Poly::operator*(const int factor) const {
    Poly result = *this;
    for (Fq& coefficient : result.coeffs) {
        coefficient = coefficient * factor;
    }
    result.trim(); // only needed if factor is a multiple of 47
    return result;
}

bool Poly::operator==(const Poly& other) const {
    return coeffs == other.coeffs;
}

Poly Poly::monic() const {
    if (is_zero()) {
        return *this;
    }
    return *this * leading().invert();
}

DivisionResult divide(const Poly& numerator, const Poly& denominator) {
    DivisionResult result;
    const int denominator_degree = denominator.degree();
    if (numerator.degree() < denominator_degree) {
        result.remainder = numerator;
        return result;
    }

    std::vector<Fq> rest = numerator.coeffs;
    result.quotient.coeffs.resize(numerator.degree() - denominator_degree + 1);
    const Fq leading_inverse = denominator.leading().invert();

    // Remove the highest coefficient of rest, one after the other
    for (int i = numerator.degree(); i >= denominator_degree; --i) {
        const Fq factor = rest[i] * leading_inverse;
        result.quotient.coeffs[i - denominator_degree] = factor;
        if (factor.is_zero()) {
            continue;
        }
        for (int j = 0; j < denominator_degree; ++j) {
            Fq& target = rest[i - denominator_degree + j];
            target = target - factor * denominator.coeffs[j];
        }
    }

    rest.resize(denominator_degree);
    result.remainder = Poly(std::move(rest));
    return result;
}

Poly gcd(const Poly& first, const Poly& second) {
    Poly larger = first;
    Poly smaller = second;
    while (!smaller.is_zero()) {
        Poly rest = divide(larger, smaller).remainder;
        larger = std::move(smaller);
        smaller = std::move(rest);
    }
    return larger.monic();
}

std::ostream& operator<<(std::ostream& os, const Poly& poly) {
    os << '[';
    for (std::size_t i = 0; i < poly.coeffs.size(); ++i) {
        if (i > 0) {
            os << ", ";
        }
        os << poly.coeffs[i];
    }
    os << ']';
    return os;
}
