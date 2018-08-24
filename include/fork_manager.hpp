#ifndef __FORK_MANAGER__
#define __FORK_MANAGER__

#include "replication_fork.hpp"
#include "util.hpp"
#include <map>
#include <vector>

class ForkManager
{
  private:
    uint n_forks, n_free_forks;
    std::vector<ReplicationFork *> replication_forks;

  public:
    ForkManager(uint size, Genome &genome, uint speed);

    /*! This function checks if there is any fork (replication) colliding with
     * any RNAP (transcription) and handles the collision by rainsing the
     * collision counter, changing the activation probability landscape around
     * the base affected and unattaching the affected fork.
     * @param time The simulation time when it was checked.
     * @param period The period of the RNAP carousel.
     * @param has_dormant Assigns if this chromosome has dormant origins.
     * @return uint The of occured collisions.
     */
    uint check_replication_transcription_conflicts(uint time, uint period,
                                                   bool has_dormant);

    /*! This function moves all forks that are attached and prepares forks just
     * unattached that were not treated and makes the available as free forks.
     * @param time The time in the simulation when the forks were advanced.
     */
    void advance_attached_forks(uint time);

    /*! This function attaches available forks to a given genomic location. If
     * there are not enough forks, the just one or none is attached.
     * @param genomic_location The location where the fork will be attached.
     * @param time The simulation time when the attachment was done.
     * @see GenomicLocation
     */
    void attach_forks(GenomicLocation &location, uint time);
};

#endif