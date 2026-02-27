#include "Fq.h"
extern const Fq a; // Coefficient 'a' of the elliptic curve equation
extern const Fq b; // Coefficient 'b' of the elliptic curve equation
struct Point {
    Fq x;
    Fq y;
    bool is_infinity;

    Point();

    Point(const Fq& x_coord, const Fq& y_coord);

    Point(const Point &other);

    Point(Point&& other) noexcept;

    Point& operator=(Point&& other) noexcept;

    Point operator+(const Point& other) const;

    Point operator-() const;
    
    Point operator-(const Point& other) const;

    Point operator*(int scalar) const;

    Point operator*(const std::array<unsigned int, 8>& scalar) const;

    bool operator==(const Point& other) const;

    bool is_on_curve() const;
};

std::ostream& operator<<(std::ostream& os, const Point& point);