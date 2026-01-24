#include <iostream>
#include "Fp.h"

int main() {
    Fp a(5);
    Fp b(10);
    Fp c = a + b;
    std::cout << "5 + 10 mod 47 = " << c.value << std::endl;
    return 0;
}