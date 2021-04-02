#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <chrono>
#include <string>

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

// S-Phase Checkpoint times for duration calculation
typedef struct
{
    std::chrono::steady_clock::time_point start_create;
    std::chrono::steady_clock::time_point end_create;
    std::chrono::steady_clock::time_point start_sim;
    std::chrono::steady_clock::time_point end_sim;
    std::chrono::steady_clock::time_point start_save;
    std::chrono::steady_clock::time_point end_save;
} s_phase_checkpoints_t;

/*
 *
 *
 *
 */
const char *getCLErrorString(int error);

#endif
