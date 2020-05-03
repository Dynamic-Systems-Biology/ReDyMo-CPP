#include "util.hpp"
#include <zstd.h>

bool operator==(const constitutive_origin_t &a, const constitutive_origin_t &b)
{
    return (a.base == b.base);
}

void *compress_cpp_string(std::string src, size_t &comp_size,
                          int comp_level) throw()
{
    const size_t est_size = ZSTD_compressBound(src.size());
    void *dest            = malloc(est_size);
    comp_size = ZSTD_compress(dest, est_size, (void *)src.c_str(), src.size(),
                              comp_level);

    if (ZSTD_isError(comp_size))
    {
        char *error_msg;
        sprintf(error_msg, "Failed to compress string. Error: %s",
                ZSTD_getErrorName(comp_size));
        throw error_msg;
    }

    return dest;
}