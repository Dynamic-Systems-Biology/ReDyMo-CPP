#ifndef __DATA_PROVIDER_HPP__
#define __DATA_PROVIDER_HPP__

#include "chromosome.hpp"
#include "util.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/*! The DataProvider class is an abstract class for passing data to a chromosome
 */
class DataProvider
{
  public:
    virtual const std::vector<std::string> &get_codes() = 0;
    virtual int get_length(std::string code)            = 0;
    virtual const std::vector<double> &
    get_probability_landscape(std::string code) = 0;
    virtual const std::vector<transcription_region_t> &
    get_transcription_regions(std::string code) = 0;
    virtual const std::shared_ptr<std::vector<constitutive_origin_t>>
    get_constitutive_origins(std::string code) = 0;
};

#endif