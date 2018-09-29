#ifndef __S_PHASE__
#define __S_PHASE__

#include "chromosome.hpp"
#include "data_manager.hpp"
#include "fork_manager.hpp"
#include "genome.hpp"
#include "util.hpp"
#include <vector>

/*! This class represents the whole Synthsis(S) phase of the cell cycle.
 *  It is a wrapper for all the other classes and serves as a main.
 */
class SPhase
{
  public:
    int origins_range, n_resources, replication_speed, timeout,
        transcription_period, time;
    bool has_dormant;
    std::shared_ptr<DataManager> data;
    std::vector<std::shared_ptr<Chromosome>> *chromosomes;
    std::shared_ptr<Genome> genome;
    std::shared_ptr<ForkManager> fork_manager;
    std::string organism;
    bool ended;
    int sim_number;

  public:
    SPhase(int origins_range, int n_resources, int replication_speed,
           int timeout, int transcription_period, bool has_dormant,
           std::shared_ptr<DataManager> data, std::string organism,
           int sim_number);
    ~SPhase();
    static void output(std::shared_ptr<SPhase> s_phase);
    void simulate();
};

#endif
