#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Ec.h"
#include "schoof.h"

// Runs Schoof's algorithm on the curve defined in Ec.cpp.
//
//   ./schoof            uses all default primes (3 .. 107 without 47), which
//                       is enough to determine the number of points
//   ./schoof 3 5 7 11   uses only the given odd primes
int main(int argc, char** argv) {
    std::vector<int> primes;
    for (int i = 1; i < argc; ++i) {
        primes.push_back(std::atoi(argv[i]));
    }
    if (primes.empty()) {
        primes = default_schoof_primes();
    }

    const auto start = std::chrono::steady_clock::now();
    const SchoofResult result = schoof(a, b, primes);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    std::cout << "\nRemainders of the trace t:\n";
    for (std::size_t i = 0; i < result.primes.size(); ++i) {
        std::cout << "  t = " << result.remainders[i] << " (mod " << result.primes[i] << ")\n";
    }
    std::cout << "Combined: t = " << to_decimal(result.trace_mod) << " (mod " << to_decimal(result.modulus) << ")\n";

    if (result.determined) {
        std::cout << "Trace: t = " << (result.trace_is_negative ? "-" : "") << to_decimal(result.trace_magnitude) << "\n";
        std::cout << "Number of points: " << to_decimal(result.order) << "\n";
        // The same number in the format Point::operator* expects
        std::cout << "As scalar (most significant word first): {" << std::hex << std::uppercase;
        for (int i = 7; i >= 0; --i) {
            std::cout << "0x" << std::setfill('0') << std::setw(8) << result.order[i] << (i > 0 ? ", " : "");
        }
        std::cout << std::dec << "}\n";
    } else {
        std::cout << "The product of the primes is not larger than 4*sqrt(q), so t itself is not determined yet.\n";
    }
    std::cout << "Time taken: " << seconds << " s" << std::endl;
    return 0;
}
