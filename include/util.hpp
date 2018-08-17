#ifndef __UTIL_HPP__
#define __UTIL_HPP__

typedef struct
{
    unsigned int start;
    unsigned int end;

} transcription_region_t;

typedef struct
{
    unsigned int base;

    bool operator==(const constitutive_origin_t &other)
    {
        return (base == other.base);
    }
} constitutive_origin_t;

#endif