#include "cuda_s_phase.hpp"
#include "fork.cuh"
#include "util.hpp"
#include <iostream>

#ifdef CUDA_ENABLED

CUDASPhase::CUDASPhase(Configuration &configuration,
                       std::shared_ptr<DataProvider> data,
                       unsigned long long seed)
    : data(data)
{
    auto args = configuration.arguments();

    origins_range        = args.constitutive;
    n_resources          = args.resources;
    replication_speed    = args.speed;
    timeout              = args.timeout;
    transcription_period = args.period;
    has_dormant          = args.dormant;
    organism             = args.organism;
    name                 = args.name;
    output_folder        = args.output;

    checkpoint_times.start_create = std::chrono::steady_clock::now();

    std::vector<std::shared_ptr<Chromosome>> chromosomes;

    auto codes = data->get_codes();
    for (auto code = codes.begin(); code != codes.end(); code++)
    {
        chromosomes.push_back(std::make_shared<Chromosome>(*code, data));
    }

    genome = std::make_shared<Genome>(chromosomes, seed);

    fork_manager =
        std::make_shared<ForkManager>(n_resources, genome, replication_speed);

    checkpoint_times.end_create = std::chrono::steady_clock::now();
}

CUDASPhase::~CUDASPhase() {}

void CUDASPhase::simulate(int sim_number)
{

    checkpoint_times.start_sim = std::chrono::steady_clock::now();

    int alpha                     = 1;
    int time                      = 0;
    int constitutive_origins      = (int)genome->n_constitutive_origins();
    int n_collisions              = 0;
    bool use_constitutive_origins = origins_range > 0;

    std::cout << "[INFO] Preparing CUDA simulation " << sim_number << std::endl;

    std::cout << "[INFO] Loading data to GPU" << std::endl;

    //////////////////////////
    // Single value buffers //
    //////////////////////////
    int *d_end_time, *d_replicated, *d_free_forks;
    cudaMalloc(&d_end_time, sizeof(int));
    cudaMalloc(&d_replicated, sizeof(int));
    cudaMalloc(&d_free_forks, sizeof(int));

    cudaMemset(d_end_time, 0, sizeof(int));
    cudaMemset(d_replicated, 0, sizeof(int));
    cudaMemcpy(d_free_forks, &n_resources, sizeof(int), cudaMemcpyHostToDevice);

    ////////////////////////////////
    // Start locations/directions //
    ////////////////////////////////
    int *d_start_locations, *d_start_directions;
    cudaMalloc(&d_start_locations, n_resources * sizeof(int));
    cudaMalloc(&d_start_directions, n_resources * sizeof(int));

    cudaMemset(d_start_locations, -1, n_resources * sizeof(int));
    cudaMemset(d_start_directions, 0, n_resources * sizeof(int));

    /////////////////////////////
    // Collision count buffers //
    /////////////////////////////
    int *d_rt_collisions;
    cudaMalloc(&d_rt_collisions, genome->size() * sizeof(int));

    cudaMemset(d_rt_collisions, 0, genome->size() * sizeof(int));

    ///////////////////////////////////////////
    // Replication timestamps for base pairs //
    ///////////////////////////////////////////
    unsigned int *d_replication_times;
    cudaMalloc(&d_replication_times, genome->size() * sizeof(unsigned int));

    cudaMemset(d_replication_times, 0, genome->size() * sizeof(unsigned int));

    ///////////////////////////////////////////
    // Probability landscape for full genome //
    ///////////////////////////////////////////
    float *d_probability_landscape;
    cudaMalloc(&d_probability_landscape, genome->size() * sizeof(float));

    std::vector<float> probability_landscape;
    for (auto chrm : genome->chromosomes)
    {
        std::vector<double> chrm_probabilities =
            data->get_probability_landscape(chrm->get_code());
        probability_landscape.insert(probability_landscape.end(),
                                     chrm_probabilities.begin(),
                                     chrm_probabilities.end());
    }
    cudaMemcpy(d_probability_landscape, probability_landscape.data(),
               probability_landscape.size() * sizeof(float),
               cudaMemcpyHostToDevice);

    /////////////////////////////////////////////////////////
    // Boundaries of chromosomes (for all the flat arrays) //
    /////////////////////////////////////////////////////////
    int *d_chromosome_boundaries;
    int boundary_count = genome->chromosomes.size() + 1;

    cudaMalloc(&d_chromosome_boundaries, boundary_count * sizeof(int));

    std::vector<int> boundaries(boundary_count, 0);

    // Generate boundaries list
    for (long int c = 0, bd = -1; c < boundary_count; c++)
    {
        boundaries[c] = bd;
        // TODO: if a bug appears, check this logic
        if (c < genome->chromosomes.size())
            bd += genome->chromosomes[c]->size();
    }
    cudaMemcpy(d_chromosome_boundaries, boundaries.data(),
               boundary_count * sizeof(int), cudaMemcpyHostToDevice);

    ///////////////////////////
    // Transcription regions //
    ///////////////////////////
    transcription_region_t *d_transcription_regions;
    std::vector<transcription_region_t> transcription_regions_v;

    for (int c = 0; c < genome->chromosomes.size(); c++)
    {
        auto regions = genome->chromosomes[c]->get_transcription_regions();

        for (int r = 0; r < regions->size(); r++)
        {
            transcription_regions_v.push_back((*regions)[r]);
        }
    }

    cudaMalloc(&d_transcription_regions,
               transcription_regions_v.size() * sizeof(transcription_region_t));

    cudaMemcpy(d_transcription_regions, transcription_regions_v.data(),
               transcription_regions_v.size() * sizeof(transcription_region_t),
               cudaMemcpyHostToDevice);

    int *d_workers_running;
    cudaMalloc(&d_workers_running, sizeof(int));

    cudaMemset(d_workers_running, 0, sizeof(int));

    ///////////////////
    // Create kernel //
    ///////////////////
    cudaDeviceSynchronize();
    std::cout << "[INFO] Launching GPU simulation " << sim_number << std::endl;
    // Add an extra thread for management
    cuda_fork<<<1, n_resources + 1>>>(
        d_end_time, d_replicated, d_free_forks, d_workers_running,
        d_start_locations, d_start_directions, d_rt_collisions,
        d_replication_times, d_probability_landscape, d_chromosome_boundaries,
        transcription_period, transcription_regions_v.size(),
        d_transcription_regions, timeout, genome->size(),
        genome->chromosomes.size(), n_resources, genome->seed, 0);

    // Wait kernel end
    cudaDeviceSynchronize();
    // TODO: add cudaFree or use cudaDeviceReset()
}

#endif
