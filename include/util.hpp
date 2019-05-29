#ifndef __UTIL_HPP__
#define __UTIL_HPP__

typedef unsigned int uint;

typedef struct
{
    int start;
    int end;

} transcription_region_t;

typedef struct
{
    int base;

} constitutive_origin_t;

bool operator==(const constitutive_origin_t &a, const constitutive_origin_t &b);
#endif