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

/*! Compresses a C++ string into a byte array using zstandard.
 * @param std::string src The source string to be compressed
 * @param size_t A reference to where the size of the compressed contents will
 * be written.
 * @param int comp_level Level of compression passed to zstd. Defaults to 1.
 * @returns void * The pointer to where the compressed result is
 * written.
 * @note This function will reallocate the dest pointer.
 */
void *compress_cpp_string(std::string src, size_t &comp_size,
                          int comp_level = 1) throw();

/*
 *
 *
 *
 */
const char *getCLErrorString(int error);

#endif