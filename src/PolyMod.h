#pragma once

#include <vector>

#include "Poly.h"
#include "U256.h"

// Result of trying to invert a polynomial modulo another one
struct InverseResult {
    // True if the inverse exists
    bool invertible = false;
    // The inverse if it exists. Otherwise the greatest common divisor of the
    // value and the modulus, which is a factor of the modulus of degree >= 1.
    Poly value;
};

// Arithmetic with polynomials modulo a fixed polynomial, i.e. in the ring
// Fq[x] / (modulus). All results have a degree below the modulus degree.
class Modulus {
public:
    // The modulus must have degree >= 1. It is divided by its leading coefficient.
    explicit Modulus(const Poly& modulus);

    int degree() const;

    const Poly& poly() const;

    // Remainder of value after division by the modulus
    Poly reduce(const Poly& value) const;

    // first * second mod modulus; both must already be reduced
    Poly mul(const Poly& first, const Poly& second) const;

    // base^exponent mod modulus
    Poly pow(const Poly& base, const U256& exponent) const;

    // 1 / value mod modulus, or a factor of the modulus if that does not exist
    InverseResult inverse(const Poly& value) const;

private:
    Poly modulus_;
    // floor(x^(2n-2) / modulus) for n = degree; used for fast reduction.
    // Empty if the modulus is too small for that to pay off.
    Poly barrett_;
};

// Computes outer(inner) mod modulus for one fixed inner polynomial and any
// number of outer polynomials. The powers of inner are computed only once.
class Composer {
public:
    // inner must be reduced. The modulus must outlive the Composer.
    Composer(const Modulus& modulus, const Poly& inner);

    // outer(inner) mod modulus; outer must be reduced
    Poly operator()(const Poly& outer) const;

private:
    const Modulus& modulus_;
    std::vector<Poly> powers_; // powers_[j] = inner^j mod modulus, j = 0..block size
};
