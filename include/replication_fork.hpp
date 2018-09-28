#ifndef __REPLICATION_FORK_HPP__
#define __REPLICATION_FORK_HPP__

#include "genome.hpp"

/*! This class represents a replication fork.
 */
class ReplicationFork
{
  private:
    std::shared_ptr<Genome> genome;
    std::shared_ptr<Chromosome> chromosome;
    uint speed;
    int base, direction;
    bool just_detached;

  public:
    /*! The constructor.
     * @param Genome genome The Genome to which this fork belongs to.
     * @param int speed The speed, in bases per step, at which the
     * fork replicates.
     */
    ReplicationFork(std::shared_ptr<Genome> genome, uint speed);

    /*! This function assigns the fork to a given base and chromosome
     * (genomic location) and replicates this base right away.
     * @param GenomicLocation genomic_location This is the position to which the
     * fork will attach.
     * @param int direction The direction to which the fork should move
     * (left or right) in a Chromosome.
     * @param int time The time when this attachment occurs in the simulation.
     * @throw invalid_argument if fork is already attached
     */
    void attach(GenomicLocation &gen_loc, int direction, uint time);

    /*! direction getter.*/
    int get_direction();

    /*! base getter.*/
    int get_base();

    /*! chromosome getter.*/
    std::shared_ptr<Chromosome> get_chromosome();

    /*! Unbinds the fork from the position where it was.
     * @param problem If the detachent is caused by a problem in replication.
     * Defaults to false.
     */
    void detach(bool problem = false);

    /*! Advances the fork proportionally to its speed and replicates the bases
     * in the path.
     * @param int time The time of the simulation when these base replications
     * happen.
     * @return True if the replication went well.
     */
    bool advance(uint time);

    /*! This function queries the attachment status of the fork.
     * @return True if the fork is attached to some base in any chromosome.
     */
    bool is_attached();

    /*! This function queries if the fork was recently detached and was not
     * treated yet.
     * @return bool True if the fork has not received post detachement
     * treatment.
     */
    bool get_just_detached();

    /*!Setter for the just_detached property.
     */
    void set_just_detached(bool new_value);
};

#endif