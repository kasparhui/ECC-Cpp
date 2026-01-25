#include "Ec.h"

const Fq a = Fq({
    Fp(45), Fp(38), Fp(38), Fp(37), Fp(38), Fp(12), Fp(25), Fp(27),
    Fp(45), Fp(34), Fp(29), Fp(4), Fp(27), Fp(7), Fp(36), Fp(18),
    Fp(0), Fp(19), Fp(23), Fp(46), Fp(10), Fp(16), Fp(7), Fp(39),
    Fp(6), Fp(22), Fp(7), Fp(20), Fp(2), Fp(26), Fp(12), Fp(37),
    Fp(38), Fp(29), Fp(36), Fp(37), Fp(10), Fp(15), Fp(15), Fp(10),
    Fp(7), Fp(41), Fp(25), Fp(22), Fp(17), Fp(3)
});
const Fq b = Fq({
    Fp(45), Fp(17), Fp(27), Fp(27), Fp(38), Fp(26), Fp(46), Fp(39),
    Fp(10), Fp(31), Fp(24), Fp(17), Fp(6), Fp(40), Fp(7),
    Fp(24), Fp(34), Fp(41), Fp(4), Fp(20), Fp(10), Fp(38), Fp(33),
    Fp(36), Fp(38), Fp(20), Fp(13), Fp(27), Fp(0), Fp(29), Fp(28),
    Fp(29), Fp(41), Fp(32), Fp(46), Fp(26), Fp(41), Fp(40), Fp(38),
    Fp(26), Fp(46), Fp(22), Fp(5), Fp(40), Fp(29), Fp(10)
});

Point::Point() : is_infinity(true) {}

Point::Point(const Fq& x_coord, const Fq& y_coord) : x(x_coord), y(y_coord), is_infinity(false) {}

Point Point::operator+(const Point& other) const {
    if (is_infinity) return other;
    if (other.is_infinity) return *this;
    if(x == other.x){
        if(y == -other.y){
            return Point(); // Point at infinity
        } else {
            // Implement point doubling logic here
            // Placeholder logic for point doubling
            Fq slope = (x * x * 3 + a) / (y * 2);
            Fq x3 = slope * slope - x*2;
            Fq y3 = slope * (x - x3) - y;
            return Point(x3, y3);
        }
    }
    else {
        // Implement point addition logic here
        // Placeholder logic for point addition
        Fq slope = (other.y - y) / (other.x - x);
        Fq x3 = slope * slope - x - other.x;
        Fq y3 = slope * (x - x3) - y;
        return Point(x3, y3);
    }

    
    return Point();
}

Point Point::operator-() const {
    if (is_infinity) return *this;

    // Implement point negation logic here
    return Point(-x, y);
}

Point Point::operator-(const Point& other) const {
    return *this + (-other);
}

Point Point::operator*(int scalar) const {
    if (is_infinity) return *this;

    // double-and-add algorithm for scalar multiplication
    Point result; // Point at infinity
    Point addend = *this;

    while (scalar > 0) {
        if (scalar & 1) {
            result = result + addend;
        }
        addend = addend + addend;
        scalar >>= 1;
    }

    return result;
}

Point Point::operator*(const std::vector<unsigned int>& scalar) const {
    if (is_infinity) return *this;
    assert(!scalar.empty());
    // double-and-add algorithm for scalar multiplication with big integer
    Point result; // Point at infinity
    Point addend = *this;

    for (int i = scalar.size()-1; i >= 0; --i) {
        unsigned int part = scalar[i];
        for (int j = 0; j < 32; ++j) { // Assuming 32 bits per part
            if (part & 1) {
                result = result + addend;
            }
            addend = addend + addend;
            part >>= 1;
        }
    }

    return result;
}

bool Point::operator==(const Point& other) const {
    if (is_infinity && other.is_infinity) return true;
    if (is_infinity || other.is_infinity) return false;
    return (x == other.x) && (y == other.y);
}

bool Point::is_on_curve() const {
    if (is_infinity) return true;

    // Check if the point satisfies the elliptic curve equation: y^2 = x^3 + ax + b
    Fq left = y * y;
    Fq right = (x * x * x) + (a * x) + b;
    return left == right;
}

std::ostream& operator<<(std::ostream& os, const Point& point) {
    if (point.is_infinity) {
        os << "Point at Infinity";
    } else {
        os << "Point(" << point.x << ", " << point.y << ")";
    }
    return os;
}