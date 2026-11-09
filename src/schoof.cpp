#include "schoof.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "DivisionPoly.h"
#include "PolyMod.h"
#include "Torsion.h"
#include "U256.h"

namespace {

Poly curve_polynomial(const Fq& curve_a, const Fq& curve_b) {
    return Poly({curve_b, curve_a, Fq(), Fq::one()}); // x^3 + ax + b
}

double seconds_since(const std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// The Frobenius map and its square on the torsion points, as polynomials
struct FrobeniusMaps {
    TorsionPoint once;  // pi(x, y)   = (x^q,     y * curve^((q-1)/2))
    TorsionPoint twice; // pi^2(x, y) = (x^(q^2), y * curve^((q^2-1)/2))
};

// Finds t mod l with pi^2(P) + q * P = t * pi(P), by trying t = 0, 1, -1, 2, -2, ...
// May throw FactorFound.
int search_trace(const TorsionGroup& group, const Modulus& modulus, const FrobeniusMaps& frobenius, const int l) {
    // q * P only depends on q mod l; use the representative closest to zero
    int q_mod_l = remainder(field_size(), l);
    const bool negative = q_mod_l > l / 2;
    if (negative) {
        q_mod_l = l - q_mod_l;
    }
    // The point P = (x, y) itself, or -P = (x, -y)
    TorsionPoint point;
    point.is_infinity = false;
    point.x = modulus.reduce(Poly::x());
    point.y = modulus.reduce(Poly::constant(negative ? -1 : 1));

    const TorsionPoint q_times_point = group.multiply(point, q_mod_l);
    const TorsionPoint target = group.add(frobenius.twice, q_times_point);

    TorsionPoint multiple; // i * pi(P), starting with i = 0
    for (int i = 0; i <= l / 2; ++i) {
        if (group.equal(target, multiple)) {
            return i;
        }
        if (group.equal(target, group.negate(multiple))) {
            return l - i;
        }
        multiple = group.add(multiple, frobenius.once);
    }
    throw std::runtime_error("schoof: no trace found for l = " + std::to_string(l));
}

// Adds the condition "t = rest (mod prime)" to "t = trace (mod modulus)"
void combine(U256& trace, U256& modulus, const int rest, const int prime) {
    // Find k in 0 .. prime-1 with trace + modulus * k = rest (mod prime)
    const int trace_mod = remainder(trace, prime);
    const int modulus_mod = remainder(modulus, prime);
    int k = 0;
    while ((trace_mod + modulus_mod * k) % prime != rest) {
        ++k;
        if (k == prime) {
            throw std::invalid_argument("schoof: the primes must be distinct");
        }
    }
    U256 step = modulus;
    const bool fits = multiply(step, static_cast<std::uint32_t>(k)) && add(trace, step) &&
                      multiply(modulus, static_cast<std::uint32_t>(prime));
    if (!fits) {
        throw std::overflow_error("schoof: the product of the primes does not fit in 256 bits");
    }
}

} // namespace

int trace_mod_2(const Fq& curve_a, const Fq& curve_b) {
    // t is even exactly if the curve has a point of order 2, i.e. a point
    // with y = 0, i.e. if x^3 + ax + b has a root in Fq. The roots in Fq of
    // any polynomial are the common roots with x^q - x.
    const Poly curve = curve_polynomial(curve_a, curve_b);
    const Modulus modulus(curve);
    const Poly x_to_q = modulus.pow(Poly::x(), field_size());
    const Poly common = gcd(x_to_q - Poly::x(), curve);
    return common.degree() == 0 ? 1 : 0;
}

int trace_mod_l(const Fq& curve_a, const Fq& curve_b, const int l, const Poly& division_polynomial) {
    if (l < 3 || l % 2 == 0 || l == p) {
        throw std::invalid_argument("trace_mod_l: l must be an odd prime other than 47");
    }
    const auto start = std::chrono::steady_clock::now();
    const U256 q = field_size();
    Modulus modulus(division_polynomial);
    const Poly curve = curve_polynomial(curve_a, curve_b);

    // pi(x, y) = (x^q, y^q) and y^q = y * (y^2)^((q-1)/2) = y * curve^((q-1)/2)
    FrobeniusMaps frobenius;
    frobenius.once.is_infinity = false;
    frobenius.once.x = modulus.pow(Poly::x(), q);
    std::cout << "  l = " << l << ": x^q done after " << seconds_since(start) << " s" << std::endl;
    frobenius.once.y = modulus.pow(curve, half(q));
    std::cout << "  l = " << l << ": y^q done after " << seconds_since(start) << " s" << std::endl;

    // pi^2 is pi applied twice. The coefficients are in Fq and therefore not
    // changed by the q-th power, so for any polynomial g: g(x)^q = g(x^q).
    //   x^(q^2) = (x^q)^q               = pi_x(pi_x)
    //   y^(q^2) = (y * pi_y(x))^q       = y * pi_y(x) * pi_y(pi_x)
    // This is much cheaper than two more exponentiations.
    {
        const Composer compose(modulus, frobenius.once.x);
        frobenius.twice.is_infinity = false;
        frobenius.twice.x = compose(frobenius.once.x);
        frobenius.twice.y = modulus.mul(compose(frobenius.once.y), frobenius.once.y);
    }
    std::cout << "  l = " << l << ": pi^2 done after " << seconds_since(start) << " s, searching t" << std::endl;

    while (true) {
        try {
            const TorsionGroup group(curve_a, curve_b, modulus);
            return search_trace(group, modulus, frobenius, l);
        } catch (const FactorFound& found) {
            // A division failed because the division polynomial is not
            // irreducible. The factor that was found describes a subset of
            // the torsion points, for which the same equation holds, so the
            // search is repeated modulo this smaller polynomial.
            std::cout << "  l = " << l << ": continuing with a factor of degree " << found.factor.degree()
                      << " after " << seconds_since(start) << " s" << std::endl;
            modulus = Modulus(found.factor);
            frobenius.once.x = modulus.reduce(frobenius.once.x);
            frobenius.once.y = modulus.reduce(frobenius.once.y);
            frobenius.twice.x = modulus.reduce(frobenius.twice.x);
            frobenius.twice.y = modulus.reduce(frobenius.twice.y);
        }
    }
}

std::vector<int> default_schoof_primes() {
    return {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97, 101, 103, 107};
}

SchoofResult combine_remainders(const std::vector<int>& primes, const std::vector<int>& remainders) {
    if (primes.size() != remainders.size()) {
        throw std::invalid_argument("combine_remainders: one remainder per prime is needed");
    }
    SchoofResult result;
    result.primes = primes;
    result.remainders = remainders;
    result.modulus = from_small(1);
    for (std::size_t i = 0; i < primes.size(); ++i) {
        combine(result.trace_mod, result.modulus, remainders[i], primes[i]);
    }

    // Hasse: |t| <= 2 * sqrt(q), and sqrt(q) = 47^23 exactly
    U256 bound = from_small(2);
    for (int n = 0; n < d / 2; ++n) {
        multiply(bound, p);
    }
    U256 range = bound; // 4 * sqrt(q), the width of the allowed interval
    add(range, bound);
    if (compare(result.modulus, range) <= 0) {
        return result; // not enough primes to pin t down
    }

    // Exactly one of trace_mod and trace_mod - modulus lies within the bound
    U256 negative_candidate = result.modulus; // modulus - trace_mod = -(trace_mod - modulus)
    subtract(negative_candidate, result.trace_mod);
    U256 q_plus_1 = field_size();
    add(q_plus_1, from_small(1));
    result.order = q_plus_1;
    if (compare(result.trace_mod, bound) <= 0) {
        result.trace_is_negative = false;
        result.trace_magnitude = result.trace_mod;
        subtract(result.order, result.trace_magnitude);
    } else if (compare(negative_candidate, bound) <= 0) {
        result.trace_is_negative = true;
        result.trace_magnitude = negative_candidate;
        add(result.order, result.trace_magnitude);
    } else {
        throw std::runtime_error("schoof: the remainders do not fit any trace allowed by Hasse's theorem");
    }
    result.determined = true;
    return result;
}

SchoofResult schoof(const Fq& curve_a, const Fq& curve_b, const std::vector<int>& odd_primes) {
    std::vector<int> primes;
    std::vector<int> remainders;

    const auto start = std::chrono::steady_clock::now();
    primes.push_back(2);
    remainders.push_back(trace_mod_2(curve_a, curve_b));
    std::cout << "l = 2: t = " << remainders.back() << " (mod 2), took " << seconds_since(start) << " s" << std::endl;

    if (!odd_primes.empty()) {
        // The l-th division polynomial is built from those with index near l/2
        const int largest = *std::max_element(odd_primes.begin(), odd_primes.end());
        const DivisionPolynomials division_polynomials(curve_a, curve_b, largest / 2 + 2);
        std::cout << "division polynomials done after " << seconds_since(start) << " s" << std::endl;

        for (const int l : odd_primes) {
            const auto prime_start = std::chrono::steady_clock::now();
            const int rest = trace_mod_l(curve_a, curve_b, l, division_polynomials.odd(l));
            primes.push_back(l);
            remainders.push_back(rest);
            std::cout << "l = " << l << ": t = " << rest << " (mod " << l << "), took "
                      << seconds_since(prime_start) << " s (total " << seconds_since(start) << " s, "
                      << primes.size() - 1 << " of " << odd_primes.size() << " odd primes done)" << std::endl;
        }
    }
    return combine_remainders(primes, remainders);
}
