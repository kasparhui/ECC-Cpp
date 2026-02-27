# ECC-Cpp
Elliptic curve cryptographie in c++
Implements Point multiplication for the field $F_{47^{46}}$  
In Main you have a point, which you can multiply by some scalar.
The scalar is an array of lenght 8 with long long unsigned int, where the scalar is the number in hex, when you concatenate the elements. With this, scalars up to $2^{256}-1$ are possible.
