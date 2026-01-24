struct Fp {
    int value;
    
    Fp(const int v);

    Fp(const Fp &other);

    Fp operator+(const Fp &other) const;

    Fp operator+(const int other) const;

    Fp operator-(const Fp& other) const;

    Fp operator-(const int other) const;

    Fp operator*(const Fp& other) const;

    Fp operator*(const int other) const;

    Fp operator/(const Fp& other) const;

    Fp operator/(const int other) const;
};

