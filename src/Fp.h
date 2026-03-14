#include <iostream>
#include <array>
#include <cassert>
#include <cstdint>

constexpr int p = 47;

constexpr std::array<uint8_t, p> inv_list = {
    0, 1, 24, 16, 12, 19, 8, 27, 6, 21, 
    33, 30, 4, 29, 37, 22, 3, 36, 34, 5, 
    40, 9, 15, 45, 2, 32, 38, 7, 42, 13, 
    11, 44, 25, 10, 18, 43, 17, 14, 26, 
    41, 20, 39, 28, 35, 31, 23, 46
};

constexpr std::array<std::array<uint8_t, p>, p> make_add_table() {
    std::array<std::array<uint8_t, p>, p> table{};

    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
            table[i][j] = (i + j) % p;

    return table;
}

constexpr auto add_table = make_add_table();

constexpr std::array<std::array<uint8_t, p>, p> make_sub_table() {
    std::array<std::array<uint8_t, p>, p> table{};

    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
            table[i][j] = (i - j + p) % p;

    return table;
}

constexpr auto sub_table = make_sub_table();

constexpr std::array<std::array<uint8_t, p>, p> make_mul_table() {
    std::array<std::array<uint8_t, p>, p> table{};

    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
            table[i][j] = (i * j) % p;

    return table;
}

constexpr auto mul_table = make_mul_table();

// constexpr combining div by multiplication with the precomputed inverses
constexpr std::array<std::array<uint8_t, p>, p> make_div_table() {
    std::array<std::array<uint8_t, p>, p> table{};

    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
            table[i][j] = (i * inv_list[j]) % p;

    return table;
}

constexpr auto div_table = make_div_table();

// A class representing an element in the finite field Fp with p = 47
struct Fp {
    uint8_t value;
    
    Fp();

    Fp(const int v);

    Fp(const Fp &other);

    void mod();

    Fp operator+(const Fp &other) const;

    Fp operator+(const int other) const;

    Fp operator-(const Fp& other) const;

    Fp operator-(const int other) const;

    Fp operator-() const;

    Fp operator*(const Fp& other) const;

    Fp operator*(const int other) const;

    Fp operator/(const Fp& other) const;

    Fp operator/(const int other) const;

    bool operator==(const Fp& other) const;

};
std::ostream& operator<<(std::ostream& os, const Fp& fp);

