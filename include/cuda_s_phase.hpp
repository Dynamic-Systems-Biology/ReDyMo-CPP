#ifndef __CUDA_S_PHASE__
#define __CUDA_S_PHASE__

#ifdef CUDA_ENABLED

#include "chromosome.hpp"
#include "data_manager.hpp"
#include "fork_manager.hpp"
#include "genome.hpp"
#include "util.hpp"
#include "s_phase.hpp"
#include <vector>

/*! This class represents the whole Synthsis(S) phase of the cell cycle.
 *  It is a wrapper for all the other classes and serves as a main.
 */
class CUDASPhase
{
  private:
    int origins_range, n_resources, replication_speed, timeout,
        transcription_period;
    bool has_dormant;

    simulation_stats stats;
    s_phase_checkpoints_t checkpoint_times;

    std::shared_ptr<DataProvider> data;
    std::vector<std::shared_ptr<Chromosome>> *chromosomes;
    std::shared_ptr<Genome> genome;
    std::shared_ptr<ForkManager> fork_manager;

    std::string organism;
    std::string name;
    std::string output_folder;

  public:
    CUDASPhase(Configuration &configuration, std::shared_ptr<DataProvider> data,
               unsigned long long seed);
    ~CUDASPhase();
    void simulate(int sim_number);

    // void output(int sim_number, int time, int iod,
    //             std::shared_ptr<Genome> genome);
    // void semantic_compression_output(int sim_number, int time, int iod,
    //                                  std::shared_ptr<Genome> genome,
    //                                  std::string path);

    // const s_phase_checkpoints_t getTimes() const;
};

#endif
#endif
