#include "fork.cuh"
#include <stdio.h>

__device__ unsigned int int_rand(unsigned int *state)
{
    unsigned int x = *state;
    // 32-bit XOR Shift RNG
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;

    *state = x;
    return x;
}

__device__ float uniform_rand(unsigned int *state)
{
    return ((float)int_rand(state)) / UINT_MAX;
}

__device__ bool collided(int transcription_period,
                         int transcription_regions_size,
                         const transcription_region_t *transcription_regions,
                         int at, int direction, int time)
{
    if (!transcription_period) return false;
    for (int i = 0; i < transcription_regions_size; i++)
    {
        transcription_region_t region = transcription_regions[i];

        int t_dir   = 1;
        int t_start = region.start;
        int t_end   = region.end;

        int replisome_position_within_region = at - t_start;

        if (region.end < region.start)
        {
            t_dir   = -1;
            t_start = region.end;
            t_end   = region.start;
        }

        // If inside transcription region
        if (!(t_start <= at && at <= t_end)) continue;

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
__device__ int get_global_id() { return blockIdx.x * blockDim.x + threadIdx.x; }

// TODO: da pra tirar a thread gerente. Cada forquilha pode sortear sua propria origem de replicacao. Da pra usar tbm uma thread por chromosomo

__global__ void cuda_fork(
    // End Time
    int *end_time,
    // Number of already replicated bases
    int *replicated,
    // Number of free forks
    int *free_forks,
    // Number of workers that are still in the main loop
    int *workers_running,
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
    const int fork_cooldown)
{
    // Replication fork ID, make manager thread -1
    int fork_id = get_global_id() - 1;

    // Local time
    int time = 0;

    // Fork cooldown
    int cooldown = 0;

    int fired_origins = 0;

    // Fork location
    int at = -1;

    // Fork direction
    int direction = 0;

    // Next chromosome boundary
    int boundary = -1;

    // Replicated
    int local_replicated = 0;

    // If fork is free
    bool free = true;

    // Initial RNG state
    unsigned int state = (unsigned int)(fork_id + 2) * seed;

    if (fork_id < 0 && true)
    {
        // Genome
        for (int i = 0; i < genome_size; i++)
            printf("%d ", replication_times[i]);
        printf("\n\n");

        // Probs
        for (int i = 0; i < genome_size; i++)
            printf("%.1f ", probability_landscape[i]);
        printf("\n\n");

        // chrm boundaries
        for (int i = 0; i < chromosome_count + 1; i++)
            printf("%d ", chromosome_boundaries[i]);
        printf("\n\n");

        // Transcr reg
        for (int i = 0; i < transcription_regions_size; i++)
            printf("(%d, %d) ", transcription_regions[i].start,
                   transcription_regions[i].end);
        printf("\n\n");

        printf("free_forks: %d\n\n", *free_forks);
        printf("end_time: %d\n\n", *end_time);
        printf("replicated: %d\n\n", *replicated);
    }

    __syncthreads();
    if (fork_id >= 0) atomicAdd(workers_running, 1);

    // TODO: see https://core.ac.uk/download/pdf/77011478.pdf
    //  Do until entire genome is replicated
    while (((fork_id < 0 && *workers_running > 0) ||
            (fork_id >= 0 && time < timeout)) &&
           (*replicated) < genome_size)
    {
        // if (*workers_running == max_forks) __syncthreads();
        // printf("[%d] n_workers_running %d\n", fork_id, *workers_running);

        if ((fork_id == 0 || fork_id == 99) && time % 50 == 0)
            printf("[%d] time %d/%d  replicated %d/%d\n", fork_id, time,
                   timeout, *replicated, genome_size);
        if (fork_id >= 0 || true)
        {
            time++;
            if (fork_id < 0 && time > 10000)
            {
                for (int i = 0; i < max_forks; i++)
                    printf("%d", start_locations[i]);
                printf("\n\n");
                break;
            }
        }

        // Try to attach forks to genome (fork manager)
        if (fork_id < 0 && (*free_forks) > 1)
        {
            printf("free forks %d\n", *free_forks);
            const int prev_free = *free_forks;
            int free_cnt        = prev_free;

            for (int attempt = 0; attempt < prev_free && free_cnt > 1;
                 attempt++)
            {
                unsigned int location = int_rand(&state) % (genome_size - 1);

                // Check if replicated and check probability
                if (!replication_times[location] &&
                    uniform_rand(&state) < probability_landscape[location])
                {
                    int i = 0;

                    // Spawn forks

                    // Fork to -1
                    for (bool attached = false; i < max_forks && !attached; i++)
                    {
                        if (start_locations[i] == -1)
                        {
                            // TODO: the actual base replicated is the next one
                            start_locations[i]  = location + 1;
                            start_directions[i] = -1;

                            attached = true;

                            break;
                        }
                    }

                    // Fork to +1
                    for (bool attached = false; i < max_forks && !attached; i++)
                    {
                        if (start_locations[i] == -1)
                        {
                            start_locations[i]  = location;
                            start_directions[i] = 1;

                            attached = true;

                            break;
                        }
                    }

                    atomicSub(free_forks, 2);
                    fired_origins += 2;
                    free_cnt -= 2;
                }
            }
        }
        if (fork_id < 0)
        {
            for (int i = 0; i < max_forks; i++)
                if (start_locations[i] != -1) printf("%d", start_locations[i]);
        }

        // If not attached, check for start_location
        if (fork_id >= 0 && !cooldown && start_locations[fork_id] > -1 &&
            at < 0)
        {
            at        = start_locations[fork_id];
            free      = false;
            direction = start_directions[fork_id];

            // Find relevant chromosome boundary
            for (int c = 0; c < chromosome_count + 1; c++)
            {
                if (chromosome_boundaries[c] > at)
                {
                    if (direction < 0)
                        boundary = chromosome_boundaries[c - 1];
                    else
                        boundary = chromosome_boundaries[c];
                }
            }
        }

        // If attached, loop copy until machinery release
        if (!cooldown && at >= 0)
        {
            int next = at + direction;

            bool on_boundary = at < 0 ? next == boundary : at == boundary;

            bool collision =
                collided(transcription_period, transcription_regions_size,
                         transcription_regions, at, direction, time);

            if (collision) atomicAdd(rt_collisions, 1);
            // If has not reached a boundary, collided with transcription
            // machinery or hit an already replicated base, continue
            if (next != boundary + direction && !collision &&
                !atomicCAS(&replication_times[next], 0, time))
            {
                local_replicated++;
                at = next;
            }
            // Set cooldown and detach fork
            else
            {
                atomicAdd(replicated, local_replicated);
                local_replicated = 0;
                cooldown         = fork_cooldown;
                at               = -1;
                float pct_rep    = (float)(*replicated) / genome_size;
                printf("Replicated: %d / %d %f%\n", *replicated, genome_size,
                       pct_rep * 100);
            }
        }

        // Fork Grace period before attaching again
        if (fork_id >= 0 && at < 0)
        {
            if (cooldown)
                cooldown--;
            else if (!free)
            {
                atomicAdd(free_forks, 1);
                start_locations[fork_id] = -1;
                free                     = true;
            }
        }
    }

    // TODO: there are 31 (yeah very odd) workers that are locking the manager
    printf("======================= thread %d out!!  // ================\n", fork_id);
    // thread inside the main loop
    __syncthreads();
    if (fork_id >= 0) atomicSub(workers_running, 1);
    if (!free) atomicAdd(free_forks, 1);

    // TODO: maybe use the largest time with cmpexch
    (*end_time) = time;
    if (fork_id < 0 && time >= timeout) printf("Timeout\n");
    if (fork_id < 0) printf("Origins Fired %d \n", fired_origins);
    if (fork_id < 0) printf("Bases Replicated %d \n", *replicated);
    if (fork_id < 0) printf("Workers at the end %d \n", *workers_running);
}
