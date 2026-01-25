#include <iostream>
#include "Fq.h"

int main() {
    std::vector<Fp> a_coeffs(46, Fp(1));
    Fq a(a_coeffs);
    std::vector<Fp> b_coeffs(46, Fp(2));
    Fq b(b_coeffs);
    Fq c = a + b;
    std::cout << "a + b = " << c << std::endl;
    return 0;
}