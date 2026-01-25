#include "Fp.h"

extern const int d; // Degree of the field extension
extern const std::vector<unsigned int> q_sub_2; // q - 2 for exponentiation as a vector of 8 32-bit words

// A class representing an element in the finite field Fq with q = 47**46 with irreducible polynomial x^46 + 2
struct Fq {
    std::vector<Fp> coeffs; // Coefficients of the polynomial representation

    Fq(const std::vector<Fp>& c);

    Fq(const Fq &other);

    Fq operator+(const Fq &other) const;

    Fq operator-(const Fq& other) const;

    Fq operator-() const;

    Fq operator*(const Fq& other) const;

    Fq operator*(const int other) const;

    Fq pow(const std::vector<unsigned int>& exp) const;

    Fq invert() const;

    Fq operator/(const Fq& other) const;
};

std::ostream& operator<<(std::ostream& os, const Fq& fq);