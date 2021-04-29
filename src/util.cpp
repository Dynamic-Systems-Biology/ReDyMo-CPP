#include "util.hpp"

bool operator==(const constitutive_origin_t &a, const constitutive_origin_t &b)
{
    return (a.base == b.base);
}
