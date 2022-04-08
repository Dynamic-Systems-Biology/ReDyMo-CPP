#include "fork.cuh"
#include "cuda_s_phase.hpp"
#include <iostream>

#ifdef CUDA_ENABLED

// CUDASPhase::CUDASPhase(int origins_range, int n_resources,
//                      int replication_speed, int timeout,
//                      int transcription_period, bool has_dormant,
//                      std::shared_ptr<DataProvider> data, std::string organism,
//                      std::string name, std::string output_folder,
//                      unsigned long long seed)
//     : origins_range(origins_range), n_resources(n_resources),
//       replication_speed(replication_speed), timeout(timeout),
//       transcription_period(transcription_period), has_dormant(has_dormant),
//       data(data), organism(organism), name(name), output_folder(output_folder),
// {
//     checkpoint_times.start_create = std::chrono::steady_clock::now();

//     std::vector<std::shared_ptr<Chromosome>> chromosomes;

//     auto codes = data->get_codes();
//     for (auto code = codes.begin(); code != codes.end(); code++)
//     {
//         chromosomes.push_back(std::make_shared<Chromosome>(*code, data));
//     }

//     genome = std::make_shared<Genome>(chromosomes, seed);

//     fork_manager =
//         std::make_shared<ForkManager>(n_resources, genome, replication_speed);

//     checkpoint_times.end_create = std::chrono::steady_clock::now();
// }

CUDASPhase::~CUDASPhase() {
}

void CUDASPhase::simulate(int sim_number){
    std::cout << "Will run CUDA simulation " << sim_number << std::endl;
    fork<<<1,100>>>(1234);
    cudaDeviceSynchronize();
}

#endif
