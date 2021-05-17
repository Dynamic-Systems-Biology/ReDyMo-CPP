#ifndef __GPU_S_PHASE__
#define __GPU_S_PHASE__

#define GPU_ENABLED

#ifdef GPU_ENABLED

#include "chromosome.hpp"
#include "data_manager.hpp"
#include "fork_manager.hpp"
#include "genome.hpp"
#include "opencl.hpp"
#include "util.hpp"
#include <vector>

/*! This class represents the whole Synthsis(S) phase of the cell cycle.
 *  It is a wrapper for all the other classes and serves as a main.
 */
class GPUSPhase
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

    // OpenCL Kernel Source Code
    std::string kernel_source;

    // Compiled OpenCL Kernel Source Code
    cl::Program kernel_program;

    // OpenCL Context
    cl::Context &clContext;

    void initialize_opencl();

  public:
    GPUSPhase(cl::Context &clContext, int origins_range, int n_resources,
              int replication_speed, int timeout, int transcription_period,
              bool has_dormant, std::shared_ptr<DataProvider> data,
              std::string organism, std::string name,
              std::string output_folder = "output", int seed = 0);
    ~GPUSPhase();
    void simulate(int sim_number);

    void output(int sim_number, int time, int iod,
                std::shared_ptr<Genome> genome);
    void semantic_compression_output(int sim_number, int time, int iod,
                                     std::shared_ptr<Genome> genome,
                                     std::string path);
    void zstd_compression_output(int sim_number, int time, int iod,
                                 std::shared_ptr<Genome> genome,
                                 std::string path);

    const s_phase_checkpoints_t getTimes() const;
};

#endif
#endif
