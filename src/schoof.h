#pragma once

#include <string>
#include <vector>

#include "Fq.h"
#include "Poly.h"
#include "U256.h"

// Schoof's algorithm for the curve y^2 = x^3 + ax + b over Fq.
//
// The number of points on the curve is q + 1 - t, where t is the trace of
// Frobenius. Schoof's algorithm finds t modulo small primes l, using the fact
// that the Frobenius map pi(x, y) = (x^q, y^q) satisfies
//   pi^2(P) - t * pi(P) + q * P = infinity
// for every point P, in particular for the points with l * P = infinity.
// The remainders are then combined with the Chinese remainder theorem.

// What is known about the trace t after combining the remainders
struct SchoofResult {
    std::vector<int> primes;     // the primes that were used
    std::vector<int> remainders; // remainders[i] = t mod primes[i]
    U256 modulus{};              // product of all primes
    U256 trace_mod{};            // t mod modulus, in 0 .. modulus - 1

    // By Hasse's theorem t lies between -2*sqrt(q) and 2*sqrt(q). If the
    // product of the primes is larger than that range (4*sqrt(q) = 4 * 47^23,
    // about 1.15 * 10^39), t is determined uniquely and the fields below are
    // set. Otherwise only t mod modulus is known.
    bool determined = false;
    bool trace_is_negative = false;
    U256 trace_magnitude{};      // |t|
    U256 order{};                // number of points on the curve, q + 1 - t
};

// t mod 2
int trace_mod_2(const Fq& curve_a, const Fq& curve_b);

// t mod l for an odd prime l other than 47. division_polynomial is the l-th
// division polynomial of the curve. Progress is written to std::cout.
int trace_mod_l(const Fq& curve_a, const Fq& curve_b, int l, const Poly& division_polynomial);

// The primes used by default: 3 .. 107 without 47. Together with 2 their
// product is about 5.5 * 10^40, the smallest such product above 4*sqrt(q),
// so the trace is determined uniquely.
std::vector<int> default_schoof_primes();

// Combines the remainders t mod primes[i] with the Chinese remainder theorem
// and, if the primes are sufficient, derives t and the number of points.
// The primes must be distinct.
SchoofResult combine_remainders(const std::vector<int>& primes, const std::vector<int>& remainders);

// Runs the whole algorithm: t mod 2, t mod l for every given odd prime, and
// the combination of all remainders. Progress is written to std::cout.
SchoofResult schoof(const Fq& curve_a, const Fq& curve_b, const std::vector<int>& odd_primes);
