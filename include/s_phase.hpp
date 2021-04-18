#ifndef __S_PHASE__
#define __S_PHASE__

#include "chromosome.hpp"
#include "configuration.hpp"
#include "data_manager.hpp"
#include "fork_manager.hpp"
#include "genome.hpp"
#include "util.hpp"
#include <vector>

/** This struct stores simulation statistics for use in the evolution simulator
 * and data output.
 */
typedef struct
{
    unsigned int collisions = 0;
    unsigned int time       = 0;
} simulation_stats;

/*! This class represents the whole Synthsis(S) phase of the cell cycle.
 *  It is a wrapper for all the other classes and serves as a main.
 */
class SPhase
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
    std::string output_folder;
    std::string name;

    void initialize(int origins_range, int n_resources, int replication_speed,
                    int timeout, int transcription_period, bool has_dormant,
                    std::shared_ptr<DataProvider> data, std::string organism,
                    std::string name, std::string output_folder);

  public:
    SPhase(int origins_range, int n_resources, int replication_speed,
           int timeout, int transcription_period, bool has_dormant,
           std::shared_ptr<DataProvider> data, std::string organism,
           std::string name, std::string output_folder = "output",
           unsigned long long seed = 0);
    SPhase(Configuration &configuration, std::shared_ptr<DataProvider> data,
           unsigned long long seed = 0);
    ~SPhase();

    simulation_stats get_stats();

    void output(int sim_number, int time, int iod,
                std::shared_ptr<Genome> genome);
    void semantic_compression_output(int sim_number, int time, int iod,
                                     std::shared_ptr<Genome> genome,
                                     std::string path);
    void zstd_compression_output(int sim_number, int time, int iod,
                                 std::shared_ptr<Genome> genome,
                                 std::string path);
    void simulate(int sim_number);

    void reset();

    const s_phase_checkpoints_t getTimes() const;
};

#endif
