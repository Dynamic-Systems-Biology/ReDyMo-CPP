/*! File Chromosome.hpp
 *  Contains the Chromosome class.
 */
#ifndef __CHROMOSOME_HPP__
#define __CHROMOSOME_HPP__

#include "data_provider.hpp"
#include "util.hpp"
#include <memory>
#include <string>
#include <vector>

// Interval between bases in the output
#define CHRM_OUTPUT_STEP 1

class GenomicLocation;

/*! The Chromosome class stores relevant data like length, number of bases
 * replicated, transcription regins and has methods to query and modify the
 * Chromosome.
 */
class Chromosome
{
    friend class GenomicLocation;
    friend class TestDataManager;

  private:
    std::string code;
    uint length;
    uint n_replicated_bases;
    uint n_fired_origins;
    std::vector<int> strand;

    std::vector<double> probability_landscape;
    const std::vector<transcription_region_t> &transcription_regions;

  public:
    std::vector<constitutive_origin_t> fired_constitutive_origins;
    const std::vector<constitutive_origin_t> &constitutive_origins;

  public:
    /*! @var int number_of_replicated_bases
     * Stores the number os bases already replicated in that Chromosome.
     */

    /*! @var int number_of_origins
     * Number of replication origins available for the Chromosome.
     */

    /*! @var vector strand
     * Stores the time when each base was replicated.
     */

    /*! The constructor for a Chromosome object.
     * @param code The id of the Chromosome.
     * @param length The length of the chromosome
     * @param probability_landscape The probability of each base to be an \
     * activation point
     * @param transcription_regions List of the transcription regions of the
     * Chromosome.
     */

    Chromosome(std::string code, DataProvider &provider);

    /*! Query the length of the Chromosome.
     * @return The length of the Chromosome.
     */
    uint size();

    /*! Print method for better visualization.
     * @return A string representation of the chromosome state.
     */
    std::string to_string();

    /*! Queries if a given base is already replicated.
     * @param base The index of a base to check.
     * @return True if the given base was replicated.
     */
    bool base_is_replicated(uint base);

    /*! Query the activation probability of a base.
     * @param int base The index of a base to check. Note that it starts at 0.
     * @return The probability of a dormant origin to attach to the given
     * base based on the probability_landscape.
     * @see probability_landscape
     */
    double activation_probability(uint base);

    /*! This method changes the probability landscape around the location of a
     * head-to-head collision. It sets the probability landscape with a
     * Gaussian function centered on the collision location.
     * @param int base The index of the base around which the
     * probability_landscape will be changed. Note that it starts at 0.
     */
    void set_dormant_activation_probability(uint base);

    /*! This function replicates the genome inside a given interval of bases.
     * This sets all the bases in the interval as replicated, and increasing
     * the number of replicated bases.
     * @see number_of_replicated_bases.
     * @param int start The index of the first base to replicate.
     * @param int end The index of the last base to replicate(exclusive).
     * @param double time The simulation iteration (time) when the replication
     * occurrs.
     * @return true if the replication did not overlap an already
     * replicated area, nor included bases outside the Chromosome.
     */
    bool replicate(int start, int end, int time);

    /*! Checks if the entire Chromosome is replicated.
     * @return True if all bases have been replicated.
     * @see base_is_replicated
     * @see replicate
     */
    bool is_replicated();

    /*! Query the id of the Chromosome.
     * @return The code of the Chromosome.
     */
    std::string get_code();

    /*! Retrieve the number of constitutive origins of this chromosome.
     * @return The number of constitutive origins.
     */
    uint n_constitutive_origins();

    /*! Getter for the number of replicated bases
     * @return uint The number of already replicated bases.
     */
    uint get_n_replicated_bases();

    uint get_n_fired_origins();

    void add_fired_origin();

    const std::vector<transcription_region_t> &
    get_transcription_regions() const;

    bool operator==(Chromosome &other);

    // Strand accessor
    int operator[](int index);
};

#endif