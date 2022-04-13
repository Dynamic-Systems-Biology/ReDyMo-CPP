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
    std::cout << "Will run CUDA simulation " << sim_number << std::endl;
    fork<<<1, 100>>>(1234);
    cudaDeviceSynchronize();

    checkpoint_times.start_sim = std::chrono::steady_clock::now();

    int alpha                     = 1;
    int time                      = 0;
    int constitutive_origins      = (int)genome->n_constitutive_origins();
    int n_collisions              = 0;
    bool use_constitutive_origins = origins_range > 0;

    std::cout << "[INFO] Starting CUDA simulation " << sim_number << std::endl;

    std::cout << "[INFO] Loading data to GPU" << std::endl;

    //////////////////////////
    // Single value buffers //
    //////////////////////////
    int *end_time, *replicated, *free_forks;
    cudaMalloc(&end_time, sizeof(int));
    cudaMalloc(&replicated, sizeof(int));
    cudaMalloc(&free_forks, sizeof(int));

    cudaMemset(end_time, 0, sizeof(int));
    cudaMemset(replicated, 0, sizeof(int));
    cudaMemset(free_forks, n_resources, n_resources * sizeof(int));

    ////////////////////////////////
    // Start locations/directions //
    ////////////////////////////////
    int *start_locations, *start_directions;
    cudaMalloc(&start_locations, n_resources * sizeof(int));
    cudaMalloc(&start_directions, n_resources * sizeof(int));

    cudaMemset(start_locations, -1, n_resources * sizeof(int));
    cudaMemset(start_directions, 0, n_resources * sizeof(int));

    /////////////////////////////
    // Collision count buffers //
    /////////////////////////////
    int *rt_collisions;
    cudaMalloc(&rt_collisions, genome->size() * sizeof(int));

    cudaMemset(rt_collisions, 0, genome->size() * sizeof(int));

    ///////////////////////////////////////////
    // Replication timestamps for base pairs //
    ///////////////////////////////////////////
    int *replication_times;
    cudaMalloc(&replication_times, genome->size() * sizeof(unsigned int));

    cudaMemset(replication_times, 0, genome->size() * sizeof(int));

    ///////////////////////////////////////////
    // Probability landscape for full genome //
    ///////////////////////////////////////////
    float *probability_landscape;
    cudaMalloc(&probability_landscape, genome->size() * sizeof(float));

    // TODO: make a flat probability landsacpe for entire genome
    // std::vector<float> probabilities = data->get_probability_landscape();
    // cudaMemcpy(replication_times, probabilities.data(), sizeof(int), );

    cudaMemset(probability_landscape, 0.f, genome->size() * sizeof(float));

    /////////////////////////////////////////////////////////
    // Boundaries of chromosomes (for all the flat arrays) //
    /////////////////////////////////////////////////////////
    int *chromosome_boundaries;
    cudaMalloc(&chromosome_boundaries, genome->size() * sizeof(int));

    std::vector<int> boundaries(genome->chromosomes.size() + 1, 0);

    // Generate boundaries list
    for (long unsigned int c = 0, bd = -1; c < genome->chromosomes.size() + 1;
         c++)
    {
        boundaries[c] = bd;
        // TODO: if a bug appears, check this logic
        if (c < genome->chromosomes.size())
            bd += genome->chromosomes[c]->size();
    }
    cudaMemcpy(replication_times, boundaries.data(),
               boundaries.size() * sizeof(int), cudaMemcpyHostToDevice);

    ///////////////////////////
    // Transcription regions //
    ///////////////////////////
    transcription_region_t *transcription_regions;
    std::vector<transcription_region_t> transcription_regions_v;

    for (int c = 0; c < genome->chromosomes.size(); c++)
    {
        auto regions = genome->chromosomes[c]->get_transcription_regions();

        for (int r = 0; r < regions->size(); r++)
        {
            transcription_regions_v.push_back((*regions)[r]);
        }
    }

    cudaMalloc(&transcription_regions,
               transcription_regions_v.size() * sizeof(transcription_region_t));

    cudaMemcpy(replication_times, boundaries.data(),
               transcription_regions_v.size() * sizeof(transcription_region_t),
               cudaMemcpyHostToDevice);

    ///////////////////
    // Create kernel //
    ///////////////////
    // Add an extra thread for management
    cuda_fork<<<1, n_resources + 1>>>();
    // TODO: add cudaFree or use cudaDeviceReset()
}

#endif
