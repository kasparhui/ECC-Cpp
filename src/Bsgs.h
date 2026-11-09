#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

// All n in 0 .. limit with n * G == H. They always form an arithmetic
// progression first, first + step, ...; count is 0 if there is none.
struct Progression {
    std::uint64_t first = 0;
    std::uint64_t step = 0;
    std::uint64_t count = 0;
};

// Baby-step giant-step: finds every n in 0 .. limit with n * G == H, using
// about 2 * sqrt(limit) group additions and sqrt(limit) table entries.
//
// It is written for any group so that it can be tested on simple ones. The
// group is described by `ops`, which must provide
//   identity(), add(x, y), negate(x), multiply(x, n), is_identity(x),
//   equal(x, y) and fingerprint(x) (a 64-bit number; equal elements must
//   have equal fingerprints, different ones may rarely collide).
template <class Element, class Operations>
Progression solve_multiples(const Element& G, const Element& H, const std::uint64_t limit, const Operations& ops) {
    Progression result;
    const auto block = static_cast<std::uint64_t>(std::ceil(std::sqrt(static_cast<double>(limit) + 1.0)));

    // Baby steps: remember j * G for j = 0 .. block - 1
    std::vector<std::pair<std::uint64_t, std::uint64_t>> table; // (fingerprint, j)
    table.reserve(block);
    Element current = ops.identity(); // j * G
    std::uint64_t order = 0;          // order of G, if it turns out to be at most block
    for (std::uint64_t j = 0; j < block; ++j) {
        if (j > 0 && ops.is_identity(current)) {
            order = j;
            break;
        }
        table.emplace_back(ops.fingerprint(current), j);
        current = ops.add(current, G);
    }
    if (order == 0 && ops.is_identity(current)) {
        order = block;
    }
    std::sort(table.begin(), table.end());

    // Every j in the table with j * G == target. A fingerprint hit is only
    // accepted after comparing the real elements.
    const auto verified_hits = [&](const Element& target, std::vector<std::uint64_t>& hits) {
        const std::uint64_t key = ops.fingerprint(target);
        auto entry = std::lower_bound(table.begin(), table.end(), std::make_pair(key, std::uint64_t(0)));
        for (; entry != table.end() && entry->first == key; ++entry) {
            if (ops.equal(ops.multiply(G, entry->second), target)) {
                hits.push_back(entry->second);
            }
        }
    };

    if (order != 0) {
        // The table holds every multiple of G; the solutions repeat with the order
        std::vector<std::uint64_t> hits;
        verified_hits(H, hits);
        if (!hits.empty() && hits.front() <= limit) {
            result.first = hits.front();
            result.step = order;
            result.count = (limit - result.first) / order + 1;
        }
        return result;
    }

    // Giant steps: H - i * block * G for i = 0, 1, ... A solution n = i * block + j
    // shows up as a table hit j in step i. Because the order of G is larger
    // than block, each step contains at most one solution, and none is missed.
    const Element giant_step = ops.negate(current); // current = block * G
    Element target = H;
    std::vector<std::uint64_t> solutions;
    for (std::uint64_t i = 0; i * block <= limit && solutions.size() < 2; ++i) {
        std::vector<std::uint64_t> hits;
        verified_hits(target, hits);
        for (const std::uint64_t j : hits) {
            if (i * block + j <= limit) {
                solutions.push_back(i * block + j);
            }
        }
        target = ops.add(target, giant_step);
    }

    if (solutions.empty()) {
        return result;
    }
    result.first = solutions[0];
    if (solutions.size() == 1) {
        result.count = 1;
    } else {
        // Two neighbouring solutions differ by exactly the order of G
        result.step = solutions[1] - solutions[0];
        result.count = (limit - result.first) / result.step + 1;
    }
    return result;
}
