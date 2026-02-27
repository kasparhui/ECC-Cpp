#include "Fp.h"

constexpr int d = 46; // Degree of the field extension
constexpr std::array<unsigned int, 8> q_sub_2 = {0xB66A744B, 0xFCDBBC9B, 0x96BBF58B, 0x7AE788BB,
                                            0xFF979CA4, 0x4C7FA532, 0xAEEBC56D, 0xE1751A5F};; // q - 2 for exponentiation as a vector of 8 32-bit words

// A class representing an element in the finite field Fq with q = 47**46 with irreducible polynomial x^46 + 2
struct Fq {
    std::array<Fp, 46> coeffs; // Coefficients of the polynomial representation

    Fq();

    Fq(const std::array<Fp, 46>& c);

    Fq(const Fq &other);

    Fq(Fq&& other) noexcept;

    Fq& operator=(Fq&& other) noexcept;

    void mod();

    Fq operator+(const Fq &other) const;

    Fq operator-(const Fq& other) const;

    Fq operator-() const;

    Fq operator*(const Fq& other) const;

    Fq operator*(const int other) const;

    Fq pow(const std::array<unsigned int, 8>& exp) const;

    Fq pow(const unsigned int exp) const;

    Fq invert() const;

    Fq operator/(const Fq& other) const;

    bool operator==(const Fq& other) const;
};

std::ostream& operator<<(std::ostream& os, const Fq& fq);