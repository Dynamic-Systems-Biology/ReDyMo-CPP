#include "fork.cuh"
#include <stdio.h>

typedef struct _transcription_region_t
{
    int start;
    int end;
} transcription_region_t;

__device__ unsigned int int_rand(unsigned int *state)
{
    printf("line: Random %d\n", __LINE__);
    unsigned int x = *state;
    // 32-bit XOR Shift RNG
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;

    *state = x;
    return x;
}

__device__ double uniform_rand(unsigned int *state)
{
    printf("line: uniform rand %d\n", __LINE__);
    return ((double)int_rand(state)) / UINT_MAX;
}

__device__ bool collided(int transcription_period,
                         int transcription_regions_size,
                         const transcription_region_t *transcription_regions, int at,
                         int direction, int time)
{
    printf("line: check collision %d\n", __LINE__);
    for (int i = 0; i < transcription_regions_size; i++)
    {
        printf("line: for each transc. reg. %d\n", __LINE__);
        transcription_region_t region = transcription_regions[i];

        int t_dir   = 1;
        int t_start = region.start;
        int t_end   = region.end;

        int replisome_position_within_region = at - t_start;

        if (region.end < region.start)
        {
            printf("line: resersed region %d\n", __LINE__);
            t_dir   = -1;
            t_start = region.end;
            t_end   = region.start;
        }

        // If inside transcription region
        if (!(t_start <= at && at <= t_end)) continue;

        printf("line: inside transc region %d\n", __LINE__);
        replisome_position_within_region = t_end - at;
        int RNAP_position =
            t_dir == 1 ? time % transcription_period
                       : transcription_period - (time % transcription_period);
        int replisome_position =
            t_dir == 1
                ? replisome_position_within_region % transcription_period
                : transcription_period -
                      (replisome_position_within_region % transcription_period);

        return replisome_position == RNAP_position && direction != t_dir;
    }
    return false;
}
__device__ int get_global_id()
{
    return blockIdx.x * blockDim.x + threadIdx.x;
}
/**
 *
 */
__global__ void fork(uint seed)
{
    uint seed_boggled;
    seed_boggled = seed ^ get_global_id();
    printf("[GPU] Hi there from %d, here is a rand %u\n", get_global_id(), int_rand(&seed_boggled));
}
