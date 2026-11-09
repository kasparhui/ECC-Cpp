#include "DivisionPoly.h"

#include <algorithm>
#include <stdexcept>

DivisionPolynomials::DivisionPolynomials(const Fq& curve_a, const Fq& curve_b, int max_index) {
    max_index = std::max(max_index, 4);
    const Fq zero;
    const Fq one = Fq::one();

    const Poly curve({curve_b, curve_a, zero, one}); // x^3 + ax + b
    curve_squared_ = curve * curve;

    table_.reserve(max_index + 1);
    // psi_0 = 0, psi_1 = 1, psi_2 = 2y
    table_.push_back(Poly());
    table_.push_back(Poly::constant(1));
    table_.push_back(Poly::constant(2));
    // psi_3 = 3x^4 + 6a x^2 + 12b x - a^2
    table_.push_back(Poly({-(curve_a * curve_a), curve_b * 12, curve_a * 6, zero, one * 3}));
    // psi_4 = 4y (x^6 + 5a x^4 + 20b x^3 - 5a^2 x^2 - 4ab x - 8b^2 - a^3)
    table_.push_back(Poly({-(curve_a * curve_a * curve_a * 4 + curve_b * curve_b * 32),
                           -(curve_a * curve_b * 16),
                           -(curve_a * curve_a * 20),
                           curve_b * 80,
                           curve_a * 20,
                           zero,
                           one * 4}));

    for (int n = 5; n <= max_index; ++n) {
        table_.push_back(n % 2 == 1 ? odd_step(n / 2) : even_step(n / 2));
    }
}

int DivisionPolynomials::max_index() const {
    return static_cast<int>(table_.size()) - 1;
}

const Poly& DivisionPolynomials::stored(const int n) const {
    return table_.at(n);
}

Poly DivisionPolynomials::odd(const int n) const {
    if (n < 1 || n % 2 == 0) {
        throw std::invalid_argument("DivisionPolynomials::odd: n must be odd and positive");
    }
    if (n <= max_index()) {
        return table_[n];
    }
    const int m = n / 2;
    if (m + 2 > max_index()) {
        throw std::invalid_argument("DivisionPolynomials::odd: n is too large for the stored table");
    }
    return odd_step(m);
}

Poly DivisionPolynomials::odd_step(const int m) const {
    // psi_(2m+1) = psi_(m+2) * psi_m^3 - psi_(m-1) * psi_(m+1)^3
    // One of the two products contains y^4 = (x^3 + ax + b)^2: the first one
    // if m is even, the second one if m is odd.
    const Poly& f_m = table_[m];
    const Poly& f_next = table_[m + 1];
    Poly first = table_[m + 2] * f_m * (f_m * f_m);
    Poly second = table_[m - 1] * f_next * (f_next * f_next);
    if (m % 2 == 0) {
        first = first * curve_squared_;
    } else {
        second = second * curve_squared_;
    }
    return first - second;
}

Poly DivisionPolynomials::even_step(const int m) const {
    // psi_(2m) = psi_m * (psi_(m+2) * psi_(m-1)^2 - psi_(m-2) * psi_(m+1)^2) / (2y)
    // The powers of y always cancel to a single y, so for the x-parts:
    //   f_(2m) = f_m * (f_(m+2) * f_(m-1)^2 - f_(m-2) * f_(m+1)^2) / 2
    const Poly& f_previous = table_[m - 1];
    const Poly& f_next = table_[m + 1];
    const Poly difference = table_[m + 2] * (f_previous * f_previous) - table_[m - 2] * (f_next * f_next);
    return (table_[m] * difference) * int(inv_list[2]);
}
