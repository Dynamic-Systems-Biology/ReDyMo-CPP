#ifndef __FORK_CUH__
#define __FORK_CUH__

#include "util.hpp"

__global__ void cuda_fork(
    // End Time
    int *end_time,
    // Number of already replicated bases
    int *replicated,
    // Number of free forks
    int *free_forks,
    // Start Location and Direction Arrays
    int *start_locations, int *start_directions,
    // Replication-transcription collision counts
    int *rt_collisions,
    // Replication timestamps for base pairs
    unsigned int *replication_times,
    // Repliction initiation probability landscape
    const float *probability_landscape,
    // Genome chromosome boundaries
    const int *chromosome_boundaries,
    // Transcription period
    const int transcription_period,
    // Transcription region count
    const int transcription_regions_size,
    // Transcription regions
    const transcription_region_t *transcription_regions,
    // Timeout
    const unsigned int timeout,
    // Total genome base pairs
    const int genome_size,
    // Chromosome count
    const int chromosome_count,
    // Number of max forks
    const int max_forks,
    // Initial seed for random
    const unsigned long seed,
    // Colldown after detach before next attempt
    const int fork_cooldown);

#endif
