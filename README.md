# ECC-Cpp
Elliptic curve cryptographie in c++
Implements point multiplication for the field $F_{47^{46}}$  
In main you have a point, which you can multiply by some scalar.
The scalar is an array of lenght 8 with unsigned int, where the scalar is the number in hex, which you get, when you concatenate the elements. With this, scalars up to $2^{256}-1$ are possible.
This implementation uses irreduceble polynomial x^46 + 2. 
