#include "Fp.h"

#include <ostream>

std::ostream& operator<<(std::ostream& os, const Fp& fp) {
    os << int(fp.value);
    return os;
}
