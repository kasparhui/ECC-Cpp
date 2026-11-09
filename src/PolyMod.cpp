#include "PolyMod.h"

#include <stdexcept>
#include <utility>

namespace {

// Fast (Barrett) reduction is used from this modulus degree on
const int barrett_min_degree = 16;

// Fast reduction is used if at least this many coefficients must be removed;
// below that, plain long division is cheaper
const int barrett_min_excess = 12;

// rest = rest mod divisor, and quotient = rest / divisor, working on raw
// coefficient lists without trailing zeros. divisor must not be empty.
void divide_in_place(std::vector<Fq>& rest, const std::vector<Fq>& divisor, std::vector<Fq>& quotient) {
    quotient.clear();
    if (rest.size() < divisor.size()) {
        return;
    }
    const std::size_t divisor_degree = divisor.size() - 1;
    quotient.resize(rest.size() - divisor_degree);
    const Fq leading_inverse = divisor.back().invert();
    for (std::size_t i = rest.size(); i-- > divisor_degree;) {
        const Fq factor = rest[i] * leading_inverse;
        quotient[i - divisor_degree] = factor;
        if (factor.is_zero()) {
            continue;
        }
        for (std::size_t j = 0; j < divisor_degree; ++j) {
            Fq& target = rest[i - divisor_degree + j];
            target = target - factor * divisor[j];
        }
    }
    rest.resize(divisor_degree);
    while (!rest.empty() && rest.back().is_zero()) {
        rest.pop_back();
    }
}

// target = target - factor * other, on raw coefficient lists
void subtract_product(std::vector<Fq>& target, const std::vector<Fq>& factor, const std::vector<Fq>& other) {
    if (factor.empty() || other.empty()) {
        return;
    }
    const std::size_t size = factor.size() + other.size() - 1;
    if (target.size() < size) {
        target.resize(size);
    }
    for (std::size_t i = 0; i < factor.size(); ++i) {
        if (factor[i].is_zero()) {
            continue;
        }
        for (std::size_t j = 0; j < other.size(); ++j) {
            target[i + j] = target[i + j] - factor[i] * other[j];
        }
    }
    while (!target.empty() && target.back().is_zero()) {
        target.pop_back();
    }
}

} // namespace

Modulus::Modulus(const Poly& modulus) : modulus_(modulus.monic()) {
    const int n = modulus_.degree();
    if (n < 1) {
        throw std::invalid_argument("Modulus: the modulus must have degree >= 1");
    }
    if (n >= barrett_min_degree) {
        std::vector<Fq> power(2 * n - 1);
        power.back() = Fq::one(); // x^(2n-2)
        barrett_ = divide(Poly(std::move(power)), modulus_).quotient;
    }
}

int Modulus::degree() const {
    return modulus_.degree();
}

const Poly& Modulus::poly() const {
    return modulus_;
}

Poly Modulus::reduce(const Poly& value) const {
    const int n = modulus_.degree();
    if (value.degree() < n) {
        return value;
    }
    const int excess = value.degree() - n + 1; // number of coefficients to remove
    if (barrett_.is_zero() || excess < barrett_min_excess || value.degree() > 2 * n - 2) {
        return divide(value, modulus_).remainder;
    }

    // Barrett reduction: with mu = floor(x^(2n-2) / modulus), the quotient of
    // value / modulus is exactly floor(floor(value / x^n) * mu / x^(n-2)).
    // This replaces the long division by two multiplications.
    const Fq* high = value.coeffs.data() + n;
    std::vector<Fq> product(excess + barrett_.coeffs.size() - 1);
    multiply(high, excess, barrett_.coeffs.data(), barrett_.coeffs.size(), product.data());

    // quotient = product / x^(n-2); it has exactly `excess` coefficients
    const Fq* quotient = product.data() + (n - 2);

    // remainder = value - quotient * modulus; only the lowest n coefficients matter
    std::vector<Fq> quotient_times_modulus(excess + modulus_.coeffs.size() - 1);
    multiply(quotient, excess, modulus_.coeffs.data(), modulus_.coeffs.size(), quotient_times_modulus.data());

    std::vector<Fq> rest(n);
    for (int i = 0; i < n; ++i) {
        rest[i] = value.coeffs[i] - quotient_times_modulus[i];
    }
    return Poly(std::move(rest));
}

Poly Modulus::mul(const Poly& first, const Poly& second) const {
    return reduce(first * second);
}

Poly Modulus::pow(const Poly& base, const U256& exponent) const {
    const Poly reduced_base = reduce(base);
    Poly result = reduce(Poly::constant(1));
    // square-and-multiply from the highest bit down. This way every
    // multiplication is by the original base, which is cheap when the base
    // is a short polynomial such as x.
    for (int index = bit_length(exponent) - 1; index >= 0; --index) {
        result = mul(result, result);
        if (get_bit(exponent, index)) {
            result = mul(result, reduced_base);
        }
    }
    return result;
}

InverseResult Modulus::inverse(const Poly& value) const {
    // Extended Euclidean algorithm. It keeps the invariant
    //   rest_k = cofactor_k * value   (mod modulus)
    // so when the rest becomes a constant, cofactor / constant is the inverse.
    std::vector<Fq> previous_rest = modulus_.coeffs;
    std::vector<Fq> rest = reduce(value).coeffs;
    std::vector<Fq> previous_cofactor;                  // 0
    std::vector<Fq> cofactor = {Fq::one()};             // 1
    std::vector<Fq> quotient;

    while (!rest.empty()) {
        divide_in_place(previous_rest, rest, quotient); // previous_rest becomes the new rest
        subtract_product(previous_cofactor, quotient, cofactor);
        std::swap(previous_rest, rest);
        std::swap(previous_cofactor, cofactor);
    }

    // previous_rest is now the greatest common divisor of value and modulus
    InverseResult result;
    if (previous_rest.size() != 1) {
        result.invertible = false;
        result.value = Poly(std::move(previous_rest)).monic();
        return result;
    }
    result.invertible = true;
    result.value = Poly(std::move(previous_cofactor)) * previous_rest[0].invert();
    return result;
}

Composer::Composer(const Modulus& modulus, const Poly& inner) : modulus_(modulus) {
    // Block size: about the square root of the number of coefficients
    int block = 1;
    while (block * block < modulus.degree()) {
        ++block;
    }
    powers_.reserve(block + 1);
    powers_.push_back(modulus.reduce(Poly::constant(1)));
    for (int j = 1; j <= block; ++j) {
        powers_.push_back(modulus.mul(powers_.back(), inner));
    }
}

Poly Composer::operator()(const Poly& outer) const {
    // Brent-Kung: cut outer into blocks of `block` coefficients,
    //   outer(x) = B_0(x) + B_1(x) * x^block + B_2(x) * x^(2*block) + ...
    // Each B_i(inner) is a sum of known powers of inner times constants, and
    // the blocks are combined with Horner's rule in inner^block. This needs
    // about 2 * sqrt(n) modular multiplications instead of n.
    const int block = static_cast<int>(powers_.size()) - 1;
    const int n = modulus_.degree();
    const int count = static_cast<int>(outer.coeffs.size());
    const int blocks = (count + block - 1) / block;

    Poly result;
    std::vector<Fq> sum(n);
    for (int index = blocks - 1; index >= 0; --index) {
        std::fill(sum.begin(), sum.end(), Fq());
        for (int j = 0; j < block && index * block + j < count; ++j) {
            const Fq& coefficient = outer.coeffs[index * block + j];
            if (coefficient.is_zero()) {
                continue;
            }
            const std::vector<Fq>& power = powers_[j].coeffs;
            for (std::size_t i = 0; i < power.size(); ++i) {
                sum[i] = sum[i] + coefficient * power[i];
            }
        }
        result = modulus_.mul(result, powers_[block]) + Poly(sum);
    }
    return result;
}
