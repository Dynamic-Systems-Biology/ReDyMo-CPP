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
    int code;
    int size;
    int n_replicated_bases;
    int n_origins;
    std::vector<float> strand;
    std::vector<float> prob_landscape;
    std::vector<bool> fired_constitutive_origins;
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

    Chromosome(int code, int length, std::vector<float> probability_landscape,
               std::vector<transcription_region_t> transcription_regions,
               std::vector<constitutive_origin_t> constitutive_origins);

    ~Chromosome();

    /*! Query the length of the Chromosome.
     * @return The length of the Chromosome.
     */
    int length();

    /*! Print method for better visualization.
     * @return A string representation of the chromosome state.
     */
    std::string to_string();

    /*! Queries if a given base is already replicated.
     * @param base The index of a base to check.
     * @return True if the given base was replicated.
     */
    bool base_is_replicated(int base);

    /*! Query the activation probability of a base.
     * @param int base The index of a base to check. Note that it starts at 0.
     * @return The probability of a dormant origin to attach to the given
     * base based on the probability_landscape.
     * @see probability_landscape
     */
    float activation_probability();

    /*! This method changes the probability landscape around the location of a
     * head-to-head collision. It sets the probability landscape with a
     * Gaussian function centered on the collision location.
     * @param int base The index of the base around which the
     * probability_landscape will be changed. Note that it starts at 0.
     */
    int set_dormant_activation_probability();

    /*! This function replicates the genome inside a given interval of bases.
     * This sets all the bases in the interval as replicated, and increasing
     * the number of replicated bases.
     * @see number_of_replicated_bases.
     * @param int start The index of the first base to replicate.
     * @param int end The index of the last base to replicate.
     * @param float time The simulation iteration (time) when the replication
     * occurrs.
     * @return true if was a normal (not at the very end) transcription.
     */
    int replicate(uint start, uint end, float time);

    /*! Checks if the entire Chromosome is replicated.
     * @return True if all bases have been replicated.
     * @see base_is_replicated
     * @see replicate
     */
    bool is_replicated();

    /*! Query the id of the Chromosome.
     * @return The code of the Chromosome.
     */
    int get_code();

    /*! Retrieve the number of constitutive origins of this chromosome.
     * @return The number of constitutive origins.
     */
    int get_n_constitutive_origins();
};

#endif