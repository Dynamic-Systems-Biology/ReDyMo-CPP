/*! File Genome.hpp
    Contains the Genome class.
*/
#ifndef __GENOME_HPP__
#define __GENOME_HPP__

#include "chromosome.hpp"
#include "genomic_location.hpp"
#include <memory>
#include <vector>

/*! This class represents a Genome.
 * It stores a set of Chromosomes and has methods to manipulate, replicate and
 * verify a Genome.
 */
class Genome
{
  public:
    std::vector<std::shared_ptr<Chromosome>> chromosomes;

  public:
    /*! Empty constructor */
    Genome();

    /*! Constructor */
    Genome(std::vector<std::shared_ptr<Chromosome>> &chromosomes);

    /*! Initializer to fill a previously empty object*/
    void initialize(std::vector<std::shared_ptr<Chromosome>> &chromosomes);

    /*! The combined length of the Genome
     * @return uint The total number of bases in the genomes.
     */
    uint size();

    /*! This function chooses a random base from a random Chromosome.
     * It does this using random integers from an uniform distribution.
     * @return A GenomicLocation object referencing the random base selected.
     * @see GenomicLocation
     */
    std::shared_ptr<GenomicLocation> random_genomic_location();

    /*! This function chooses a random UNREPLICATED base from a random
     *Chromosome. It assumes that the genome is not completely replicated yet.
     * It does this using random integers from an uniform distribution.
     * @return A GenomicLocation object referencing the random base selected.
     * @see GenomicLocation
     */
    std::shared_ptr<GenomicLocation> random_unreplicated_genomic_location();

    /*! Checks if the Genome is entirely replicated.
     * It checks if all Chromosomes are completely replicated.
     * @return True if all bases of all Chromosomes have been replicated.
     * @see Chromosome
     */
    bool is_replicated();

    /*! This function calculates the average inter-origin distance across all
     * Chromosomes in the Genome.
     * @return The average inter-origin distance measured in number of bases.
     */
    float average_interorigin_distance();

    /*! Retrieve the number of constitutive origins in the whole genome.
     * @return The number of constitutive origins in the whole genome.
     */
    uint n_constitutive_origins();
};

#endif