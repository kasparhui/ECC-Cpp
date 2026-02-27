#include "Fp.h"

extern const int d; // Degree of the field extension
extern const std::array<unsigned int, 8> q_sub_2; // q - 2 for exponentiation as a vector of 8 32-bit words

// A class representing an element in the finite field Fq with q = 47**46 with irreducible polynomial x^46 + 2
struct Fq {
    std::array<Fp, 46> coeffs; // Coefficients of the polynomial representation

    Fq();

    Fq(const std::array<Fp, 46>& c);

    Fq(const Fq &other);

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