#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "Bsgs.h"
#include "Ec.h"
#include "schoof.h"

// Finishing the point count with a search (baby-step giant-step).
//
// Schoof's algorithm gives t modulo the product M of the primes. If M is not
// large enough to pin t down, the possible orders are
//   N = q + 1 - t,   t = t0 + k * M,   |t| <= 2 * sqrt(q)
// The real order N satisfies N * P = infinity for every point P, so the right
// k can be searched for with points of the curve. That is much cheaper than
// running Schoof's algorithm for the largest primes.
//
// Everything here works on the curve defined in Ec.h / Ec.cpp.

// Square root in Fq. Returns false if value is not a square; otherwise root
// is one of the two square roots.
bool square_root(const Fq& value, Fq& root);

// A random point of the curve (never the point at infinity)
Point random_point(std::mt19937_64& generator);

// Baby-step giant-step on the curve: every n in 0 .. limit with n * G == H
Progression solve_multiples(const Point& G, const Point& H, std::uint64_t limit);

enum class SearchStatus {
    Unique,    // exactly one order is possible; result has been completed
    Ambiguous, // several orders fit all tested points; more primes are needed
    TooLarge   // too many candidates for a search; more primes are needed
};

// Tries to determine t and the order from the remainders in result. On
// Unique, result.determined is set and the trace and order are filled in.
// The order is only accepted if it is the single candidate that annihilates
// the tested points, and it is checked again on a fresh random point.
SearchStatus complete_with_search(SchoofResult& result);

// Counts the points of the curve: Schoof's algorithm for the given odd
// primes, then the search. If the search cannot decide, further primes are
// added one at a time until the order is certain.
SchoofResult count_points(const std::vector<int>& odd_primes, unsigned threads = 0);
