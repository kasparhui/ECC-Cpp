#include "OrderSearch.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "DivisionPoly.h"
#include "U256.h"

namespace {

// More candidates than this are not searched (memory and time grow with the
// square root: 10^14 candidates need 10^7 table entries, about 160 MB)
const std::uint64_t search_limit = 100000000000000ULL;

// Up to this many remaining candidates are tested one by one
const std::uint64_t explicit_limit = 100000;

// Point::operator* and Fq::pow expect the most significant word first
std::array<unsigned int, 8> to_scalar(const U256& value) {
    std::array<unsigned int, 8> scalar{};
    for (int i = 0; i < 8; ++i) {
        scalar[i] = value[7 - i];
    }
    return scalar;
}

// Constants for square roots in Fq (Tonelli-Shanks): q - 1 = 2^s * m, m odd
struct SquareRootConstants {
    int s = 0;
    U256 m{};
    U256 m_plus_1_half{};  // (m + 1) / 2
    U256 q_minus_1_half{}; // (q - 1) / 2
    Fq generator;          // z^m for a non-square z; has order 2^s
};

SquareRootConstants make_square_root_constants() {
    SquareRootConstants constants;
    U256 q_minus_1 = field_size();
    subtract(q_minus_1, from_small(1));
    constants.q_minus_1_half = half(q_minus_1);

    constants.m = q_minus_1;
    while (!get_bit(constants.m, 0)) {
        constants.m = half(constants.m);
        ++constants.s;
    }
    constants.m_plus_1_half = constants.m;
    add(constants.m_plus_1_half, from_small(1));
    constants.m_plus_1_half = half(constants.m_plus_1_half);

    // Find a non-square z: z^((q-1)/2) = -1. Half of all elements qualify.
    const Fq minus_one = -Fq::one();
    for (int k = 0;; ++k) {
        Fq z;
        z.coeffs[1] = Fp(1);
        z.coeffs[0] = Fp(k % p);
        z.coeffs[2] = Fp(k / p);
        if (z.pow(to_scalar(constants.q_minus_1_half)) == minus_one) {
            constants.generator = z.pow(to_scalar(constants.m));
            return constants;
        }
    }
}

const SquareRootConstants& square_root_constants() {
    static const SquareRootConstants constants = make_square_root_constants();
    return constants;
}

// A number identifying a point, for the lookup table. Different points can
// share a number in rare cases, so every hit is verified afterwards.
std::uint64_t point_fingerprint(const Point& point) {
    if (point.is_infinity) {
        return 0;
    }
    static_assert(sizeof(Fp) == 1, "Fp must be a single byte");
    std::uint64_t hash = 1469598103934665603ULL; // FNV-1a
    for (const Fq* coordinate : {&point.x, &point.y}) {
        for (const Fp& coefficient : coordinate->coeffs) {
            hash = (hash ^ coefficient.value) * 1099511628211ULL;
        }
    }
    return hash;
}

double seconds_since(const std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

Point scalar_multiple(const Point& point, const U256& factor) {
    return point * to_scalar(factor);
}

} // namespace

bool square_root(const Fq& value, Fq& root) {
    if (value.is_zero()) {
        root = Fq();
        return true;
    }
    const SquareRootConstants& constants = square_root_constants();
    const Fq one = Fq::one();
    // Euler's criterion: value is a square exactly if value^((q-1)/2) = 1
    if (!(value.pow(to_scalar(constants.q_minus_1_half)) == one)) {
        return false;
    }

    // Tonelli-Shanks. Invariant: root^2 = value * t, and t has order 2^(e-1)
    // or less; each round lowers e until t = 1.
    int e = constants.s;
    Fq c = constants.generator;
    Fq t = value.pow(to_scalar(constants.m));
    root = value.pow(to_scalar(constants.m_plus_1_half));
    while (!(t == one)) {
        // smallest i with t^(2^i) = 1
        int i = 0;
        for (Fq power = t; !(power == one); power = power * power) {
            ++i;
        }
        Fq b = c;
        for (int k = 0; k < e - i - 1; ++k) {
            b = b * b;
        }
        root = root * b;
        c = b * b;
        t = t * c;
        e = i;
    }
    return true;
}

Point random_point(std::mt19937_64& generator) {
    while (true) {
        Fq x;
        for (Fp& coefficient : x.coeffs) {
            coefficient = Fp(static_cast<int>(generator() % p));
        }
        // y^2 = x^3 + ax + b has a solution for about half of all x
        Fq y;
        if (square_root(x * x * x + a * x + b, y)) {
            return Point(x, y);
        }
    }
}

namespace {

// The curve group, in the form the generic search in Bsgs.h expects
struct CurveOperations {
    Point identity() const { return Point(); }
    Point add(const Point& first, const Point& second) const { return first + second; }
    Point negate(const Point& point) const { return -point; }
    Point multiply(const Point& point, const std::uint64_t factor) const {
        return scalar_multiple(point, from_uint64(factor));
    }
    bool is_identity(const Point& point) const { return point.is_infinity; }
    bool equal(const Point& first, const Point& second) const { return first == second; }
    std::uint64_t fingerprint(const Point& point) const { return point_fingerprint(point); }
};

} // namespace

Progression solve_multiples(const Point& G, const Point& H, const std::uint64_t limit) {
    return solve_multiples(G, H, limit, CurveOperations());
}

SearchStatus complete_with_search(SchoofResult& result) {
    if (result.determined) {
        return SearchStatus::Unique;
    }
    const auto start = std::chrono::steady_clock::now();

    // Hasse: |t| <= bound = 2 * sqrt(q) = 2 * 47^23
    U256 bound = from_small(2);
    for (int n = 0; n < d / 2; ++n) {
        multiply(bound, p);
    }
    const U256& modulus = result.modulus;
    const U256& t0 = result.trace_mod;

    // Candidates t = t0 + k * modulus with -below <= k <= above
    U256 numerator = bound;
    add(numerator, t0);
    const U256 below = quotient(numerator, modulus);
    U256 above{}; // only used if t0 <= bound; otherwise k must be negative
    const bool nonnegative_allowed = compare(t0, bound) <= 0;
    if (nonnegative_allowed) {
        numerator = bound;
        subtract(numerator, t0);
        above = quotient(numerator, modulus);
    }
    if (!fits_in_64_bits(below) || !fits_in_64_bits(above) || to_uint64(below) > search_limit ||
        to_uint64(above) > search_limit) {
        return SearchStatus::TooLarge;
    }
    if (!nonnegative_allowed && to_uint64(below) == 0) {
        throw std::runtime_error("order search: the remainders do not fit any trace allowed by Hasse's theorem");
    }
    // Number the candidates n = 0 .. limit, starting at the smallest t:
    //   t_n = t_smallest + n * modulus,   N_n = N_largest - n * modulus
    const std::uint64_t limit = nonnegative_allowed ? to_uint64(below) + to_uint64(above) : to_uint64(below) - 1;
    if (limit > search_limit) {
        return SearchStatus::TooLarge;
    }

    U256 largest_order = field_size(); // q + 1 - t0 + below * modulus
    add(largest_order, from_small(1));
    subtract(largest_order, t0);
    U256 shift = modulus;
    multiply64(shift, to_uint64(below));
    add(largest_order, shift);

    const auto candidate_order = [&](const std::uint64_t n) {
        U256 order = largest_order;
        U256 offset = modulus;
        multiply64(offset, n);
        subtract(order, offset);
        return order;
    };

    // N_n * P = infinity  <=>  n * (modulus * P) = N_largest * P
    std::mt19937_64 generator(20240917); // fixed seed: the run is reproducible
    Progression best;
    for (int attempt = 0; attempt < 5; ++attempt) {
        const Point point = random_point(generator);
        const Progression found = solve_multiples(scalar_multiple(point, modulus), scalar_multiple(point, largest_order), limit);
        if (found.count == 0) {
            throw std::runtime_error("order search: no possible order annihilates a point; a remainder must be wrong");
        }
        if (best.count == 0 || found.count < best.count) {
            best = found;
        }
        if (best.count <= explicit_limit) {
            break;
        }
    }
    if (best.count > explicit_limit) {
        schoof_log("search: too many orders fit the tested points");
        return SearchStatus::Ambiguous;
    }

    // Test the remaining candidates one by one against further points
    std::vector<std::uint64_t> candidates;
    for (std::uint64_t i = 0; i < best.count; ++i) {
        candidates.push_back(best.first + i * best.step);
    }
    for (int round = 0; round < 10 && candidates.size() > 1; ++round) {
        const Point point = random_point(generator);
        std::vector<std::uint64_t> kept;
        for (const std::uint64_t n : candidates) {
            if (scalar_multiple(point, candidate_order(n)).is_infinity) {
                kept.push_back(n);
            }
        }
        candidates = std::move(kept);
    }
    if (candidates.empty()) {
        throw std::runtime_error("order search: no possible order annihilates all points; a remainder must be wrong");
    }
    if (candidates.size() > 1) {
        std::ostringstream line;
        line << "search: " << candidates.size() << " orders fit all tested points";
        schoof_log(line.str());
        return SearchStatus::Ambiguous;
    }

    // Exactly one candidate is left. The real order is among the candidates
    // and annihilates every point, so this must be it. Check it once more.
    const U256 order = candidate_order(candidates[0]);
    if (!scalar_multiple(random_point(generator), order).is_infinity) {
        throw std::runtime_error("order search: the remaining order fails on a fresh point; a remainder must be wrong");
    }

    U256 q_plus_1 = field_size();
    add(q_plus_1, from_small(1));
    result.order = order;
    result.trace_is_negative = compare(order, q_plus_1) > 0;
    if (result.trace_is_negative) {
        result.trace_magnitude = order;
        subtract(result.trace_magnitude, q_plus_1);
    } else {
        result.trace_magnitude = q_plus_1;
        subtract(result.trace_magnitude, order);
    }
    result.determined = true;
    {
        std::ostringstream line;
        line << "search among " << limit + 1 << " possible orders: " << std::fixed << std::setprecision(1)
             << seconds_since(start) << " s";
        schoof_log(line.str());
    }
    return SearchStatus::Unique;
}

SchoofResult count_points(const std::vector<int>& odd_primes, const unsigned threads) {
    SchoofResult result = schoof(a, b, odd_primes, threads);

    // Normally the search settles it. If not, add the next prime and try again;
    // this ends at the latest when the primes alone determine the trace.
    while (complete_with_search(result) != SearchStatus::Unique) {
        int next = result.primes.empty() ? 3 : *std::max_element(result.primes.begin(), result.primes.end());
        const auto is_usable = [&](const int candidate) {
            if (candidate == p || candidate % 2 == 0) {
                return false;
            }
            for (int divisor = 3; divisor * divisor <= candidate; divisor += 2) {
                if (candidate % divisor == 0) {
                    return false;
                }
            }
            return std::find(result.primes.begin(), result.primes.end(), candidate) == result.primes.end();
        };
        do {
            ++next;
        } while (next < 3 || !is_usable(next));

        schoof_log("the search could not decide; adding the prime " + std::to_string(next));
        const auto start = std::chrono::steady_clock::now();
        const DivisionPolynomials division_polynomials(a, b, next / 2 + 2);
        const int rest = trace_mod_l(a, b, next, division_polynomials.odd(next));
        {
            std::ostringstream line;
            line << "l = " << next << ": t = " << rest << "   " << std::fixed << std::setprecision(1)
                 << seconds_since(start) << " s";
            schoof_log(line.str());
        }
        std::vector<int> primes = result.primes;
        std::vector<int> remainders = result.remainders;
        primes.push_back(next);
        remainders.push_back(rest);
        result = combine_remainders(primes, remainders);
    }
    return result;
}
