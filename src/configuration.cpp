#include "configuration.hpp"

#include <c4/yml/std/string.hpp>
#include <fstream>
#include <random>

conf_function_map cl_evolution_mutations_genes_functions = {
    PUSH_FUNCS(evolution.mutations.genes.move,
               (conf_function_map{
                   PUSH_D(evolution.mutations.genes.move.prob),
                   PUSH_D(evolution.mutations.genes.move.std),
               })),
    PUSH_FUNCS(evolution.mutations.genes.swap,
               (conf_function_map{
                   PUSH_D(evolution.mutations.genes.swap.prob),
               })),
};

conf_function_map cl_evolution_mutations_probability_landscape_functions = {
    PUSH_D(evolution.mutations.probability_landscape.add),
    PUSH_D(evolution.mutations.probability_landscape.del),
    PUSH_FUNCS(
        evolution.mutations.probability_landscape.change_mean,
        (conf_function_map{
            PUSH_D(evolution.mutations.probability_landscape.change_mean.prob),
            PUSH_D(evolution.mutations.probability_landscape.change_mean.std),
        })),
    PUSH_FUNCS(
        evolution.mutations.probability_landscape.change_std,
        (conf_function_map{
            PUSH_D(evolution.mutations.probability_landscape.change_std.prob),
            PUSH_D(evolution.mutations.probability_landscape.change_std.std),
            PUSH_D(evolution.mutations.probability_landscape.change_std.max),
        })),
};

conf_function_map cl_evolution_functions = {
    PUSH_ULL(evolution.population), PUSH_ULL(evolution.generations),
    PUSH_ULL(evolution.survivors),
    PUSH_FUNCS(
        evolution.mutations,
        (conf_function_map{
            PUSH_FUNCS(evolution.mutations.probability_landscape,
                       cl_evolution_mutations_probability_landscape_functions),
            PUSH_FUNCS(evolution.mutations.genes,
                       cl_evolution_mutations_genes_functions),
        })),
    PUSH_FUNCS(evolution.fitness,
               (conf_function_map{
                   PUSH_D(evolution.fitness.min_sphase),
                   PUSH_D(evolution.fitness.match_mfaseq),
                   PUSH_D(evolution.fitness.max_coll_all),
                   PUSH_D(evolution.fitness.min_coll_all),
                   PUSH_FUNCS(evolution.fitness.max_coll,
                              (conf_function_map{
                                  PUSH_D(evolution.fitness.max_coll.weight),
                                  PUSH_STR(evolution.fitness.max_coll.gene),
                              })),

                   PUSH_FUNCS(evolution.fitness.min_coll,
                              (conf_function_map{
                                  PUSH_D(evolution.fitness.max_coll.weight),
                                  PUSH_STR(evolution.fitness.max_coll.gene),
                              }))}))};

conf_function_map cl_configuration_functions = {
    PUSH_ULL(cells),        PUSH_STR(organism),
    PUSH_ULL(resources),    PUSH_ULL(speed),
    PUSH_ULL(timeout),      PUSH_BOOL(dormant),
    PUSH_STR(name),         PUSH_ULL(period),
    PUSH_ULL(constitutive), PUSH_STR(data_dir),
    PUSH_D(probability),    PUSH_STR(output),
    PUSH_ULL(threads),      PUSH_FUNCS(evolution, cl_evolution_functions)};

void read_conf_yml(ryml::NodeRef &base, cl_configuration_data &arguments,
                   conf_function_map &function_map,
                   std::function<void(std::string)> on_unknown)
{
    for (auto ref : base.children())
    {
        auto key = std::string();
        auto val = std::string();

        c4::from_chars(ref.key(), &key);
        if (ref.has_val()) c4::from_chars(ref.val(), &val);

        // Call callback on unknown
        if (function_map.find(key) == function_map.end())
            on_unknown(key);
        else // Else run function
            function_map[key](arguments, val, base);
    }
}

/**
 * Class definitions
 */

Configuration::Configuration(int argc, char *argv[])
{
    args = configure_cmd_options(argc, argv);
}

cl_configuration_data Configuration::arguments() { return args; }

cl_configuration_data Configuration::configure_cmd_options(int argc,
                                                           char *argv[])
{
    int c;
    cl_configuration_data arguments;

    int dormant = -1;
    int summary = 0;

    std::string config;

    while (1)
    {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"gpu", no_argument, 0, 'g'},
            {"cells", required_argument, 0, 'c'},
            {"organism", required_argument, 0, 'o'},
            {"resources", required_argument, 0, 'r'},
            {"speed", required_argument, 0, 's'},
            {"dormant", no_argument, &dormant, 1},
            {"summary", no_argument, &summary, 1},

            {"seed", required_argument, 0, 'x'},
            {"name", required_argument, 0, 'n'},
            {"config", required_argument, 0, 'C'},
            {"timeout", required_argument, 0, 'T'},
            {"period", required_argument, 0, 'P'},
            {"constitutive", required_argument, 0, 'k'},
            {"data-dir", required_argument, 0, 'd'},
            {"probability", required_argument, 0, 'p'},
            {"output", required_argument, 0, 'O'},
            {"threads", required_argument, 0, 't'},
            {NULL, 0, NULL, 0}};

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv,
                        "h:g:c:o:r:s:T:DP:n:C:d:p:O:t:x:", long_options,
                        &option_index);

        /* Detect the end of the options. */
        if (c == -1) break;

        switch (c)
        {
        case 0: break;
        case 'h':
            // TODO: add help
            throw std::invalid_argument("Usage: ");
            break;
        case 'g':
            std::cout << "GPU Processing is broken for now. Ignoring option..."
                      << std::endl;
            break;
        case 'c': arguments.cells = std::stoull(optarg); break;
        case 'o': arguments.organism = std::string(optarg); break;
        case 'r': arguments.resources = std::stoull(optarg); break;
        case 's': arguments.speed = std::stoull(optarg); break;
        case 'T': arguments.timeout = std::stoull(optarg); break;

        case 'x': arguments.seed = std::stoull(optarg); break;
        case 'n': arguments.name = std::string(optarg); break;
        case 'C': config = std::string(optarg); break;
        case 'P': arguments.period = std::stoull(optarg); break;
        case 'k': arguments.constitutive = std::stoull(optarg); break;
        case 'd': arguments.data_dir = std::string(optarg); break;
        case 'p': arguments.probability = std::stod(optarg); break;
        case 'O': arguments.output = std::string(optarg); break;
        case 't': arguments.threads = std::stoull(optarg); break;

        case '?':
            /* getopt_long already printed an error message. */
            throw std::invalid_argument("Check arguments");
            break;

        default: abort();
        }
    }

    if (config.length() > 0) read_configuration_file(config, arguments);

    if (dormant >= 0) arguments.dormant = !!dormant;

    if (!arguments.cells)
    {
        throw std::invalid_argument("Argument \"cells\" (c) is mandatory!");
    }

    if (!arguments.organism.length())
    {
        throw std::invalid_argument("Argument \"organism\" (o) is mandatory!");
    }

    if (!arguments.resources)
    {
        throw std::invalid_argument("Argument \"resources\" (r) is mandatory!");
    }

    if (!arguments.timeout)
    {
        throw std::invalid_argument("Argument \"timeout\" (T) is mandatory!");
    }

    // Set a GLOBAL seed based on a random device (possibly true random number)
    if (!arguments.seed)
    {
        std::random_device rd;
        arguments.seed = rd();
    }

    if (summary)
    {
        std::cout << std::endl
                  << "Parameter summary "
                     "================================================"
                  << std::endl
                  << std::flush;
        std::cout << "Number of cells         : " << arguments.cells
                  << std::endl
                  << std::flush;
        std::cout << "Organism                : " << arguments.organism
                  << std::endl
                  << std::flush;
        std::cout << "Number of forks         : " << arguments.resources
                  << std::endl
                  << std::flush;
        std::cout << "Steps per iteration     : " << arguments.speed
                  << std::endl
                  << std::flush;
        std::cout << "Max Iterations          : " << arguments.timeout
                  << std::endl
                  << std::flush;
        std::cout << "Use dormant origins     : "
                  << (arguments.dormant ? "Yes" : "No") << std::endl
                  << std::flush;
        std::cout << "Transcription period    : " << arguments.period
                  << std::endl
                  << std::flush;
        std::cout << "Use constitutive origins: " << arguments.constitutive
                  << std::endl
                  << std::flush;
        std::cout << "Data directory          : " << arguments.data_dir
                  << std::endl
                  << std::flush;
        if (arguments.probability != 0)
            std::cout << "Uniform probability     : " << arguments.probability
                      << std::endl
                      << std::flush;
        std::cout << "Output directory        : " << arguments.output
                  << std::endl
                  << std::flush;
        std::cout << "Thread count            : " << arguments.threads
                  << std::endl
                  << std::flush;
        std::cout << "Seed for the RNG        : "
                  << arguments.seed << std::endl
                  << std::flush;
    }

    return arguments;
}

cl_configuration_data
Configuration::read_configuration_file(std::string filename,
                                       cl_configuration_data &arguments)
{
    /** field setters
     *
     */

    std::ifstream config_file(filename);
    std::string str((std::istreambuf_iterator<char>(config_file)),
                    std::istreambuf_iterator<char>());

    c4::substr buffer = c4::to_substr(str);

    ryml::Tree tree = ryml::parse(buffer);

    ryml::NodeRef simulation = tree["simulation"];
    ryml::NodeRef parameters = tree["parameters"];

    std::string mode = std::string();
    c4::from_chars(simulation.val(), &mode);

    arguments.mode = mode;
    // root configuration handler
    read_conf_yml(parameters, arguments, cl_configuration_functions);

    return arguments;
}

bool operator==(const cl_evolution_data &a, const cl_evolution_data &b)
{
    return a.generations == b.generations;
}

bool operator==(const cl_configuration_data &a, const cl_configuration_data &b)
{
    return a.mode == b.mode && a.cells == b.cells && a.organism == b.organism &&
           a.resources == b.resources && a.speed == b.speed &&
           a.timeout == b.timeout && a.dormant == b.dormant &&
           a.seed == b.seed && a.name == b.name && a.period == b.period &&
           a.constitutive == b.constitutive && a.data_dir == b.data_dir &&
           a.probability == b.probability && a.output == b.output &&
           a.threads == b.threads && a.evolution == b.evolution;
}
