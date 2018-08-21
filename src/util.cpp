#include "util.hpp"

bool operator==(const constitutive_origin_t &a, const constitutive_origin_t &b)
{
    return (b.base == b.base);
}