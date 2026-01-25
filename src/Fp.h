#include <iostream>
#include <vector>
#include <cassert>

extern const std::vector<int> inv_list;
extern const int p;
// A class representing an element in the finite field Fp with p = 47
struct Fp {
    int value;
    
    Fp();

    Fp(const int v);

    Fp(const Fp &other);

    Fp operator+(const Fp &other) const;

    Fp operator+(const int other) const;

    Fp operator-(const Fp& other) const;

    Fp operator-(const int other) const;

    Fp operator-() const;

    Fp operator*(const Fp& other) const;

    Fp operator*(const int other) const;

    Fp operator/(const Fp& other) const;

    Fp operator/(const int other) const;

};
std::ostream& operator<<(std::ostream& os, const Fp& fp);

