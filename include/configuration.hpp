/*! File Evolution.hpp
 *  Contains the Evolution class.
 */
#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

#include <algorithm>
#include <functional>
#include <getopt.h>
#include <ryml.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#define PUSH_ULL(field)                                                        \
    {                                                                          \
        ([]() -> std::string {                                                 \
            std::string path = #field;                                         \
            return path.substr(path.rfind(".") + 1, path.size());              \
        })(),                                                                  \
            [](cl_configuration_data &data, std::string val,                   \
               ryml::NodeRef &base) -> void { data.field = std::stoull(val); } \
    }
#define PUSH_D(field)                                                          \
    {                                                                          \
        ([]() -> std::string {                                                 \
            std::string path = #field;                                         \
            return path.substr(path.rfind(".") + 1, path.size());              \
        })(),                                                                  \
            [](cl_configuration_data &data, std::string val,                   \
               ryml::NodeRef &base) -> void { data.field = std::stod(val); }   \
    }
#define PUSH_STR(field)                                                        \
    {                                                                          \
        ([]() -> std::string {                                                 \
            std::string path = #field;                                         \
            return path.substr(path.rfind(".") + 1, path.size());              \
        })(),                                                                  \
            [](cl_configuration_data &data, std::string val,                   \
               ryml::NodeRef &base) -> void { data.field = std::string(val); } \
    }
#define PUSH_BOOL(field)                                                       \
    {                                                                          \
        ([]() -> std::string {                                                 \
            std::string path = #field;                                         \
            return path.substr(path.rfind(".") + 1, path.size());              \
        })(),                                                                  \
            [](cl_configuration_data &data, std::string val,                   \
               ryml::NodeRef &base) {                                          \
                data.field = !std::string(val).compare("true");                \
            }                                                                  \
    }

#define PUSH_FUNCS(field, functions)                                           \
    {                                                                          \
        ([]() -> std::string {                                                 \
            std::string path = #field;                                         \
            return path.substr(path.rfind(".") + 1, path.size());              \
        })(),                                                                  \
            [](cl_configuration_data &data, std::string val,                   \
               ryml::NodeRef &base) {                                          \
                std::string path = #field;                                     \
                std::string k =                                                \
                    path.substr(path.rfind(".") + 1, path.size()).c_str();     \
                conf_function_map map = functions;                             \
                for (ryml::NodeRef c : base.children())                        \
                {                                                              \
                    if (!c.key().compare(k.c_str(), k.size()))                 \
                    {                                                          \
                        read_conf_yml(c, data, map);                           \
                        break;                                                 \
                    }                                                          \
                }                                                              \
            }                                                                  \
    }

typedef struct
{
    std::string optimizer;
    double weight;
    std::string gene;
} cl_fitness_option;

typedef struct
{
    unsigned long long population  = 0;
    unsigned long long generations = 0;
    unsigned long long survivors   = 0;

    struct
    {
        struct
        {
            double add = 0;
            double del = 0;
            struct
            {
                double prob = 0;
                double std  = 0;
            } change_mean;

            struct
            {
                double prob = 0;
                double std  = 0;
                double max  = 0;
            } change_std;
        } probability_landscape;

        struct
        {
            struct
            {
                double prob = 0;
                double std  = 0;
            } move;

            struct
            {
                double prob = 0;
            } swap;
        } genes;
    } mutations;

    struct
    {
        double min_sphase   = 0;
        double match_mfaseq = 0;

        double max_coll_all = 0;
        double min_coll_all = 0;

        struct
        {
            double weight = 0;
            std::string gene;
        } max_coll;

        struct
        {
            double weight = 0;
            std::string gene;
        } min_coll;
    } fitness;
} cl_evolution_data;

typedef struct
{
    std::string mode = "basic";

    unsigned long long cells     = 0;
    std::string organism         = "";
    unsigned long long resources = 0;
    unsigned long long speed     = 1;
    unsigned long long timeout   = 0;
    bool dormant                 = false;

    unsigned long long seed         = 0;
    std::string name                = "no_name";
    unsigned long long period       = 0;
    unsigned long long constitutive = 0;
    std::string data_dir            = "../data";
    double probability              = 0;
    std::string output              = "output";
    unsigned long long threads      = 8;
    bool gpu                        = false;

    // Other modes data
    cl_evolution_data evolution;
} cl_configuration_data;

bool operator==(const cl_evolution_data &a, const cl_evolution_data &b);

bool operator==(const cl_configuration_data &a, const cl_configuration_data &b);

typedef std::unordered_map<
    std::string,
    std::function<void(cl_configuration_data &, std::string, ryml::NodeRef &)>>
    conf_function_map;

void read_conf_yml(
    ryml::NodeRef &base, cl_configuration_data &arguments,
    conf_function_map &function_map,
    std::function<void(std::string)> on_unknown = [](std::string argument) {
        throw std::invalid_argument(
            "Unknown parameter in configuration file: " + argument);
    });

/*! This class represents a running configuration for the simulations.
 *
 * It contains parses and stores all configuration options for passing
 * to simulations classes.
 */
class Configuration
{
  private:
    // Allow the test class to set args directly
    friend class EvolutionTest;

    cl_configuration_data args;

    cl_configuration_data configure_cmd_options(int argc, char *argv[]);

    cl_configuration_data
    read_configuration_file(std::string filename,
                            cl_configuration_data &arguments);

  public:
    Configuration(int argc, char *argv[]);

    cl_configuration_data arguments();
};

#endif
