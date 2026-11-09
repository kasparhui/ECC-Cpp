#include "schoof.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <future>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
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

// The progress line for one finished prime, e.g. "l =  73: t =  15   65.6 s"
std::string prime_line(const int l, const int trace, const double seconds) {
    std::ostringstream line;
    line << "l = " << std::setw(3) << l << ": t = " << std::setw(3) << trace << "   " << std::fixed
         << std::setprecision(1) << std::setw(5) << seconds << " s";
    return line.str();
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

void schoof_log(const std::string& line) {
    static std::mutex mutex;
    const std::lock_guard<std::mutex> lock(mutex);
    std::cout << line << std::endl;
}

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
    const U256 q = field_size();
    Modulus modulus(division_polynomial);
    const Poly curve = curve_polynomial(curve_a, curve_b);

    // pi(x, y) = (x^q, y^q) and y^q = y * (y^2)^((q-1)/2) = y * curve^((q-1)/2).
    // The two powers do not depend on each other, so y^q runs on a second
    // thread while this one computes x^q.
    FrobeniusMaps frobenius;
    frobenius.once.is_infinity = false;
    {
        std::future<Poly> y_part = std::async(std::launch::async, [&]() { return modulus.pow(curve, half(q)); });
        frobenius.once.x = modulus.pow(Poly::x(), q);
        frobenius.once.y = y_part.get();
    }

    // pi^2 is pi applied twice. The coefficients are in Fq and therefore not
    // changed by the q-th power, so for any polynomial g: g(x)^q = g(x^q).
    //   x^(q^2) = (x^q)^q               = pi_x(pi_x)
    //   y^(q^2) = (y * pi_y(x))^q       = y * pi_y(x) * pi_y(pi_x)
    // This is much cheaper than two more exponentiations. The two
    // substitutions share the powers of pi_x and again run side by side.
    {
        const Composer compose(modulus, frobenius.once.x);
        frobenius.twice.is_infinity = false;
        std::future<Poly> y_part = std::async(std::launch::async, [&]() { return compose(frobenius.once.y); });
        frobenius.twice.x = compose(frobenius.once.x);
        frobenius.twice.y = modulus.mul(y_part.get(), frobenius.once.y);
    }

    while (true) {
        try {
            const TorsionGroup group(curve_a, curve_b, modulus);
            return search_trace(group, modulus, frobenius, l);
        } catch (const FactorFound& found) {
            // A division failed because the division polynomial is not
            // irreducible. The factor that was found describes a subset of
            // the torsion points, for which the same equation holds, so the
            // search is repeated modulo this smaller polynomial.
            modulus = Modulus(found.factor);
            frobenius.once.x = modulus.reduce(frobenius.once.x);
            frobenius.once.y = modulus.reduce(frobenius.once.y);
            frobenius.twice.x = modulus.reduce(frobenius.twice.x);
            frobenius.twice.y = modulus.reduce(frobenius.twice.y);
        }
    }
}

std::vector<int> default_schoof_primes() {
    return {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 59, 61, 67, 71, 73};
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

SchoofResult schoof(const Fq& curve_a, const Fq& curve_b, const std::vector<int>& odd_primes,
                    unsigned threads) {
    const auto start = std::chrono::steady_clock::now();
    const int t2 = trace_mod_2(curve_a, curve_b);
    schoof_log(prime_line(2, t2, seconds_since(start)));

    // Largest prime first: the large primes take by far the longest, so they
    // should start right away instead of waiting behind the small ones.
    std::vector<int> order = odd_primes;
    std::sort(order.begin(), order.end(), [](const int first, const int second) { return first > second; });
    std::vector<int> rests(order.size(), 0);

    if (!order.empty()) {
        // The l-th division polynomial is built from those with index near l/2
        const DivisionPolynomials division_polynomials(curve_a, curve_b, order.front() / 2 + 2);

        if (threads == 0) {
            threads = std::max(1u, std::thread::hardware_concurrency());
        }
        threads = std::min<unsigned>(threads, static_cast<unsigned>(order.size()));
        schoof_log("computing " + std::to_string(order.size()) + " odd primes on " + std::to_string(threads) +
                   " threads, largest first");

        std::atomic<std::size_t> next{0};      // index of the next prime nobody has taken yet
        std::atomic<std::size_t> finished{0};
        std::exception_ptr failure;
        std::mutex failure_mutex;

        // Each thread keeps taking the next prime until none are left
        const auto worker = [&]() {
            while (true) {
                const std::size_t index = next.fetch_add(1);
                if (index >= order.size()) {
                    return;
                }
                const int l = order[index];
                try {
                    const auto prime_start = std::chrono::steady_clock::now();
                    rests[index] = trace_mod_l(curve_a, curve_b, l, division_polynomials.odd(l));
                    const double seconds = seconds_since(prime_start);
                    schoof_log(prime_line(l, rests[index], seconds) + "   (" +
                               std::to_string(finished.fetch_add(1) + 1) + " of " + std::to_string(order.size()) + ")");
                } catch (...) {
                    const std::lock_guard<std::mutex> lock(failure_mutex);
                    if (!failure) {
                        failure = std::current_exception();
                    }
                    next.store(order.size()); // let the other threads stop early
                    return;
                }
            }
        };

        std::vector<std::thread> pool;
        for (unsigned i = 1; i < threads; ++i) {
            pool.emplace_back(worker);
        }
        worker(); // the calling thread works too
        for (std::thread& thread : pool) {
            thread.join();
        }
        if (failure) {
            std::rethrow_exception(failure);
        }
    }

    // Report the primes in increasing order, starting with 2
    std::vector<int> primes = {2};
    std::vector<int> remainders = {t2};
    for (std::size_t i = order.size(); i-- > 0;) {
        primes.push_back(order[i]);
        remainders.push_back(rests[i]);
    }
    return combine_remainders(primes, remainders);
}
