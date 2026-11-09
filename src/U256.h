#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>

#include "Fq.h"

// A 256-bit unsigned integer as eight 32-bit words, least significant word
// first. Only the few operations Schoof's algorithm needs are provided.
// Operations that can overflow report it instead of wrapping around silently.
using U256 = std::array<std::uint32_t, 8>;

// q = 47^46, the number of elements of Fq
inline U256 field_size() {
    U256 value{};
    value[0] = 1;
    for (int n = 0; n < d; ++n) {
        std::uint64_t carry = 0;
        for (std::uint32_t& word : value) {
            carry += static_cast<std::uint64_t>(word) * p;
            word = static_cast<std::uint32_t>(carry);
            carry >>= 32;
        }
    }
    return value;
}

inline bool get_bit(const U256& value, const int index) {
    return (value[index / 32] >> (index % 32)) & 1;
}

// Number of bits needed to write the value; 0 for the value 0
inline int bit_length(const U256& value) {
    for (int index = 255; index >= 0; --index) {
        if (get_bit(value, index)) {
            return index + 1;
        }
    }
    return 0;
}

// value / 2, rounded down
inline U256 half(const U256& value) {
    U256 result{};
    for (int i = 0; i < 8; ++i) {
        result[i] = value[i] >> 1;
        if (i + 1 < 8) {
            result[i] |= value[i + 1] << 31;
        }
    }
    return result;
}

// value mod modulus, for a small positive modulus
inline int remainder(const U256& value, const int modulus) {
    std::uint64_t rest = 0;
    for (int i = 7; i >= 0; --i) {
        rest = ((rest << 32) | value[i]) % static_cast<std::uint64_t>(modulus);
    }
    return static_cast<int>(rest);
}

// -1, 0 or 1 if first is smaller than, equal to or larger than second
inline int compare(const U256& first, const U256& second) {
    for (int i = 7; i >= 0; --i) {
        if (first[i] != second[i]) {
            return first[i] < second[i] ? -1 : 1;
        }
    }
    return 0;
}

// value += other; returns false if the result does not fit in 256 bits
inline bool add(U256& value, const U256& other) {
    std::uint64_t carry = 0;
    for (int i = 0; i < 8; ++i) {
        carry += static_cast<std::uint64_t>(value[i]) + other[i];
        value[i] = static_cast<std::uint32_t>(carry);
        carry >>= 32;
    }
    return carry == 0;
}

// value -= other; other must not be larger than value
inline void subtract(U256& value, const U256& other) {
    std::int64_t borrow = 0;
    for (int i = 0; i < 8; ++i) {
        const std::int64_t difference = static_cast<std::int64_t>(value[i]) - other[i] - borrow;
        value[i] = static_cast<std::uint32_t>(difference);
        borrow = difference < 0 ? 1 : 0;
    }
}

// value *= factor; returns false if the result does not fit in 256 bits
inline bool multiply(U256& value, const std::uint32_t factor) {
    std::uint64_t carry = 0;
    for (std::uint32_t& word : value) {
        carry += static_cast<std::uint64_t>(word) * factor;
        word = static_cast<std::uint32_t>(carry);
        carry >>= 32;
    }
    return carry == 0;
}

inline U256 from_small(const std::uint32_t small) {
    U256 value{};
    value[0] = small;
    return value;
}

// Decimal representation
inline std::string to_decimal(U256 value) {
    std::string digits;
    do {
        // value /= 10, keeping the remainder as the next digit
        std::uint64_t rest = 0;
        for (int i = 7; i >= 0; --i) {
            const std::uint64_t current = (rest << 32) | value[i];
            value[i] = static_cast<std::uint32_t>(current / 10);
            rest = current % 10;
        }
        digits.push_back(static_cast<char>('0' + rest));
    } while (compare(value, U256{}) != 0);
    std::reverse(digits.begin(), digits.end());
    return digits;
}

// value *= factor for a 64-bit factor; returns false on overflow
inline bool multiply64(U256& value, const std::uint64_t factor) {
    U256 high_part = value;
    const bool low_fits = multiply(value, static_cast<std::uint32_t>(factor));
    // high_part = value * (factor >> 32) * 2^32
    bool high_fits = multiply(high_part, static_cast<std::uint32_t>(factor >> 32));
    if (high_part[7] != 0) {
        high_fits = false;
    }
    for (int i = 7; i > 0; --i) {
        high_part[i] = high_part[i - 1];
    }
    high_part[0] = 0;
    const bool sum_fits = add(value, high_part);
    return low_fits && high_fits && sum_fits;
}

// floor(numerator / denominator); the denominator must not be zero
inline U256 quotient(const U256& numerator, const U256& denominator) {
    U256 result{};
    U256 rest{};
    // binary long division, from the highest bit down
    for (int index = 255; index >= 0; --index) {
        // rest = rest * 2 + bit
        for (int i = 7; i > 0; --i) {
            rest[i] = (rest[i] << 1) | (rest[i - 1] >> 31);
        }
        rest[0] = (rest[0] << 1) | (get_bit(numerator, index) ? 1u : 0u);
        if (compare(rest, denominator) >= 0) {
            subtract(rest, denominator);
            result[index / 32] |= 1u << (index % 32);
        }
    }
    return result;
}

// True if the value fits in 64 bits
inline bool fits_in_64_bits(const U256& value) {
    for (int i = 2; i < 8; ++i) {
        if (value[i] != 0) {
            return false;
        }
    }
    return true;
}

inline std::uint64_t to_uint64(const U256& value) {
    return (static_cast<std::uint64_t>(value[1]) << 32) | value[0];
}

inline U256 from_uint64(const std::uint64_t small) {
    U256 value{};
    value[0] = static_cast<std::uint32_t>(small);
    value[1] = static_cast<std::uint32_t>(small >> 32);
    return value;
}
