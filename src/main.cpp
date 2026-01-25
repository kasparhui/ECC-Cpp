#include <iostream>
#include "Fq.h"

int main() {
    std::vector<Fp> a_coeffs(46, Fp(1));
    Fq a(a_coeffs);
    std::vector<Fp> b_coeffs(46, Fp(2));
    Fq b(b_coeffs);
    Fq c = a / b;
    Fq d = a * b;
    Fq e = c * d;
    std::cout << "e: "<< e << std::endl;
    std::cout << "a * a: "<< a * a << std::endl;
    return 0;
}