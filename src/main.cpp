#include <iostream>
#include "Ec.h"

int main() {
    Fq x = Fq({
        Fp(15), Fp(46), Fp(9), Fp(33), Fp(3), Fp(31), Fp(14), Fp(8),
        Fp(19), Fp(46), Fp(21), Fp(30), Fp(3), Fp(45), Fp(20), Fp(8),
        Fp(22), Fp(45), Fp(45), Fp(33), Fp(37), Fp(20), Fp(31), Fp(24),
        Fp(9), Fp(37), Fp(36), Fp(8), Fp(38), Fp(44), Fp(11), Fp(43),
        Fp(34), Fp(43), Fp(25), Fp(22), Fp(0), Fp(13), Fp(32), Fp(23),
        Fp(5), Fp(9), Fp(17), Fp(7), Fp(4), Fp(7)
    });

    Fq y = Fq({
        Fp(14), Fp(19), Fp(30), Fp(37), Fp(45), Fp(16), Fp(38), Fp(28),
        Fp(35), Fp(24), Fp(40), Fp(39), Fp(6), Fp(25), Fp(40), Fp(10),
        Fp(39), Fp(36), Fp(17), Fp(41), Fp(28), Fp(10), Fp(28), Fp(17),
        Fp(43), Fp(42), Fp(4), Fp(22), Fp(5), Fp(15), Fp(37), Fp(23),
        Fp(35), Fp(7), Fp(23), Fp(43), Fp(34), Fp(42), Fp(45), Fp(13),
        Fp(1), Fp(19), Fp(8), Fp(25), Fp(13), Fp(41)
    });


    Point P(x, y);
    std::vector<unsigned int> scalar = {
        0xB66A744B, 0xFCDBBC9B, 0x96BBF58B, 0x7AE788BD,
        0x51D52358, 0x2D2B640A, 0x72A64822, 0x726C5FAB
    };
    std::cout << "P on curve: " << P.is_on_curve() << std::endl;
    Point R = P * scalar;
    std::cout << "Point scalar: " << R << std::endl;
    std::cout << "R on curve: " << R.is_on_curve() << std::endl;
    return 0;
}