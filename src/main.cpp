#include <iostream>
#include "Fp.h"

int main() {
    Fp a(35);
    Fp b(43);
    Fp c = a  / b;
    std::cout << "5 / 10 mod 47 = " << c<< std::endl;
    Fp d = c * b;
    std::cout << " (5 / 10) * 10 mod 47 = " << d << std::endl;
    return 0;
}