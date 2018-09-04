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
  private:
    int origins_range, n_resources, replication_speed, timeout,
        transcription_period;
    bool has_dormant;
    DataManager *data;
    std::vector<Chromosome *> *chromosomes;
    Genome *genome;
    ForkManager *fork_manager;
    std::string organism;

  public:
    SPhase(int origins_range, int n_resources, int replication_speed,
           int timeout, int transcription_period, bool has_dormant,
           DataManager *data, std::string organism);
    ~SPhase();
    void output(int sim_number, int time, int iod, Genome *genome);
    void simulate(int sim_number);
};

#endif
