#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>

constexpr int p = 47;

// inv_list[v] is the multiplicative inverse of v mod 47 (inv_list[0] is unused)
constexpr std::array<uint8_t, p> inv_list = {
    0, 1, 24, 16, 12, 19, 8, 27, 6, 21,
    33, 30, 4, 29, 37, 22, 3, 36, 34, 5,
    40, 9, 15, 45, 2, 32, 38, 7, 42, 13,
    11, 44, 25, 10, 18, 43, 17, 14, 26,
    41, 20, 39, 28, 35, 31, 23, 46
};

// An element of the finite field Fp with p = 47. The value is always in 0..46.
// Everything is defined here in the header so the compiler can inline it.
struct Fp {
    uint8_t value = 0;

    constexpr Fp() = default;

    // Any integer, reduced into 0..46 (negative values included)
    constexpr Fp(const int v) : value(reduce(v)) {}

    // A value that is already in 0..46; not reduced
    constexpr Fp(const uint8_t v) : value(v) {}

    constexpr Fp operator+(const Fp& other) const {
        const int sum = value + other.value;
        return Fp(static_cast<uint8_t>(sum >= p ? sum - p : sum));
    }

    constexpr Fp operator-(const Fp& other) const {
        const int diff = value - other.value;
        return Fp(static_cast<uint8_t>(diff < 0 ? diff + p : diff));
    }

    constexpr Fp operator-() const {
        return Fp(static_cast<uint8_t>(value == 0 ? 0 : p - value));
    }

    constexpr Fp operator*(const Fp& other) const {
        return Fp(static_cast<uint8_t>((value * other.value) % p));
    }

    // Division by zero is not defined; it returns 0
    constexpr Fp operator/(const Fp& other) const {
        return *this * Fp(inv_list[other.value]);
    }

    constexpr Fp operator+(const int other) const { return *this + Fp(other); }
    constexpr Fp operator-(const int other) const { return *this - Fp(other); }
    constexpr Fp operator*(const int other) const { return *this * Fp(other); }
    constexpr Fp operator/(const int other) const { return *this / Fp(other); }

    constexpr bool operator==(const Fp& other) const {
        return value == other.value;
    }

private:
    static constexpr uint8_t reduce(const int v) {
        const int r = v % p;
        return static_cast<uint8_t>(r < 0 ? r + p : r);
    }
};

std::ostream& operator<<(std::ostream& os, const Fp& fp);
