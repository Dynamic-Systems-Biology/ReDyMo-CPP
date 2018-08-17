#ifndef __GENOMIC_LOCATION_HPP__
#define __GENOMIC_LOCATION_HPP__

#include "chromosome.hpp"
#include "util.hpp"

/*! This class represents a base from a chromosome called a Genomic Location.
 * It stores a set of Chromosomes and has methods to check and query a base
 * from a Genome.
 */
class GenomicLocation
{
  private:
    unsigned int base;
    Chromosome chromosome;

  public:
    /*! The constructor.
     * @param int base The index of the base that this location represents.
     * @param Chromosome chromosome The Chromosome which contains the given
     * base.
     * @throw invalid_argument If base doesn't belong to chromosome.
     * @see Chromosome
     */
    GenomicLocation(unsigned int base, Chromosome &chromosome);

    /*! Queries if the particular base has been replicated.
     * @return True if the base is replicated.
     */
    bool is_replicated();

    /*! Tests the probability of the base to be activated.
     * @param Boolean use_constitutive_origins True if uses this type of origin.
     * @param int origins_range Considered range around a constitutive origin.
     * @return True if the base will be activated.
     */
    bool will_activate(bool use_constitutive_origin,
                       unsigned int origins_range);

    /*! Retrieve a constitutive origin located in this range.
     * @param int origins_range Considered range around a constitutive origin.
     * @return True if the base will be activate, False otherwise.
     * @throw range_error if no origin is found within the range.
     */
    constitutive_origin_t *get_constitutive_origin(unsigned int origins_range);

    /*! Update the list of fired constitutive origins with an fired origin.
     * @param constitutive_origin_t location of the constitutive origin.
     * @return True if list was successfully updated and False otherwise.
     * @throw invalid_argument if origin doesn't exist in list of origins.
     */
    bool put_fired_constitutive_origin(constitutive_origin_t &origin);
};

#endif