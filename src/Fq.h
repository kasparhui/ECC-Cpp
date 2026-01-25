#include "Fp.h"

extern const int d; // Degree of the field extension

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

    Fq operator/(const Fq& other) const;
};

std::ostream& operator<<(std::ostream& os, const Fq& fq);