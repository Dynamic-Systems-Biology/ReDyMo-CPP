#ifndef __UTIL_HPP__
#define __UTIL_HPP__

typedef unsigned int uint;

typedef struct
{
    uint start;
    uint end;

} transcription_region_t;

typedef struct
{
    uint base;

} constitutive_origin_t;

bool operator==(const constitutive_origin_t &a, const constitutive_origin_t &b);
#endif