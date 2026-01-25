#include <iostream>
#include <vector>

extern const std::vector<int> inv_list;

struct Fp {
    int value;
    
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

