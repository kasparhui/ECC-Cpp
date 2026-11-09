#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Ec.h"
#include "OrderSearch.h"
#include "schoof.h"

// Counts the points of the curve defined in Ec.cpp.
//
//   ./schoof              Schoof's algorithm for the default primes (3 .. 73
//                         without 47), then a search for the exact order
//   ./schoof 3 5 7 11     Schoof's algorithm for the given odd primes only;
//                         prints the remainders of the trace
//   ./schoof -j 4 ...     use 4 threads instead of one per CPU core
int main(int argc, char** argv) {
    std::vector<int> primes;
    unsigned threads = 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-j" && i + 1 < argc) {
            threads = static_cast<unsigned>(std::atoi(argv[++i]));
        } else {
            primes.push_back(std::atoi(argv[i]));
        }
    }

    const auto start = std::chrono::steady_clock::now();
    const SchoofResult result = primes.empty() ? count_points(default_schoof_primes(), threads)
                                               : schoof(a, b, primes, threads);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    std::cout << "\nt = " << to_decimal(result.trace_mod) << " (mod " << to_decimal(result.modulus) << ")\n";

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
        std::cout << "These primes alone do not determine t. Run without arguments for the full point count.\n";
    }
    std::cout << "Time taken: " << std::fixed << std::setprecision(1) << seconds << " s" << std::endl;
    return 0;
}
