/*! File Chromosome.hpp
    Contains the Chromosome class.
*/
#ifndef __CHROMOSOME_HPP__
#define __CHROMOSOME_HPP__

#include "util.hpp"
#include <string>
#include <vector>

/*! The Chromosome class stores relevant data like length, number of bases
 * replicated, transcription regins and has methods to query and modify the
 * Chromosome.
 */
class Chromosome
{
  private:
    unsigned int code;
    unsigned int length;
    unsigned int n_replicated_bases;
    unsigned int n_fired_origins;
    std::vector<int> strand;
    std::vector<float> probability_landscape;
    std::vector<int> fired_constitutive_origins;
    std::vector<transcription_region_t> transcription_regions;
    std::vector<constitutive_origin_t> constitutive_origins;

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

    Chromosome(unsigned int code, unsigned int length, std::vector<float> &probability_landscape,
               std::vector<transcription_region_t> &transcription_regions,
               std::vector<constitutive_origin_t> &constitutive_origins);
    
    Chromosome();

    void initialize(unsigned int code, unsigned int length, std::vector<float> &probability_landscape,
               std::vector<transcription_region_t> &transcription_regions,
               std::vector<constitutive_origin_t> &constitutive_origins);

    /*! Query the length of the Chromosome.
     * @return The length of the Chromosome.
     */
    unsigned int size();

    /*! Print method for better visualization.
     * @return A string representation of the chromosome state.
     */
    std::string to_string();

    /*! Queries if a given base is already replicated.
     * @param base The index of a base to check.
     * @return True if the given base was replicated.
     */
    bool base_is_replicated(unsigned int base);

    /*! Query the activation probability of a base.
     * @param int base The index of a base to check. Note that it starts at 0.
     * @return The probability of a dormant origin to attach to the given
     * base based on the probability_landscape.
     * @see probability_landscape
     */
    float activation_probability(unsigned int base);

    /*! This method changes the probability landscape around the location of a
     * head-to-head collision. It sets the probability landscape with a
     * Gaussian function centered on the collision location.
     * @param int base The index of the base around which the
     * probability_landscape will be changed. Note that it starts at 0.
     */
    void set_dormant_activation_probability(unsigned int base);

    /*! This function replicates the genome inside a given interval of bases.
     * This sets all the bases in the interval as replicated, and increasing
     * the number of replicated bases.
     * @see number_of_replicated_bases.
     * @param int start The index of the first base to replicate.
     * @param int end The index of the last base to replicate(exclusive).
     * @param float time The simulation iteration (time) when the replication
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
    unsigned int get_code();

    /*! Retrieve the number of constitutive origins of this chromosome.
     * @return The number of constitutive origins.
     */
    unsigned int n_constitutive_origins();
};

#endif
