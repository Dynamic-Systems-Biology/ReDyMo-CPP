#include "s_phase.hpp"
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>

int get_cmd_option(char **input, int size, const std::string &option)
{
    for (int i = 0; i < size; i++)
        if (option == input[i]) return ++i;
    return -1;
}

bool cmd_option_exists(char **begin, char **end, const std::string &option)
{
    return std::find(begin, end, option) != end;
}

void print_arg_opt(std::string opt, std::string value)
{
    std::cout << std::left << std::setw(26) << opt << std::right
              << std::setw(40) << value << std::endl;
}

std::vector<std::string> configure_cmd_options(int argc, char *argv[])
{
    std::vector<std::string> arg_options = {
        "--cells",  "--organism",     "--resources",
        "--speed",  "--timeout",      "--dormant",
        "--period", "--constitutive", "--data-folder"};

    std::vector<std::string> optional_args = {"--period", "--constitutive",
                                              "--data-folder"};

    std::vector<std::string> arg_texts = {"Number of cells: ",
                                          "Organism: ",
                                          "Number of forks: ",
                                          "Steps per iteration: ",
                                          "Max Iterations: ",
                                          "Use dormant origins: ",
                                          "Transcription period: ",
                                          "Use constitutive origins: ",
                                          "Data path: "};

    std::vector<std::string> arg_values(arg_options.size());

    for (auto option : arg_options)
    {
        bool is_optional = std::find(optional_args.begin(), optional_args.end(),
                                     option) != optional_args.end();

        if (!cmd_option_exists(argv, argv + argc, option) && !is_optional)
        {
            std::cout << "Incorrect arguments."
                      << " Required argument " << option << " was not given."
                      << " For usage see README.md" << std::endl;
            exit(0);
        }
    }

    std::cout
        << std::endl
        << "Parameter summary ================================================"
        << std::endl;

    arg_values[6] = "0";
    arg_values[7] = "0";
    // Add default location for backward compatibility
    arg_values[8] = "../data";

    for (int i = 0; i < arg_options.size(); i++)
    {

        if (cmd_option_exists(argv, argv + argc, arg_options[i]))
        {
            int index     = get_cmd_option(argv, argc, arg_options[i]);
            arg_values[i] = argv[index];
        }

        print_arg_opt(arg_texts[i], arg_values[i]);
    }

    return arg_values;
}

int main(int argc, char *argv[])
{
    std::vector<std::string> arg_values = configure_cmd_options(argc, argv);

    std::string data_dir = arg_values[8];

    std::shared_ptr<DataManager> data = std::make_shared<DataManager>(
        data_dir + "/database.sqlite", data_dir + "/MFA-Seq_TBrucei_TREU927/");

    srand(time(NULL));

    omp_set_num_threads(40);

    int n_cells = std::stoi(arg_values[0]);
    bool has_dormant;
    std::istringstream(arg_values[5]) >> std::boolalpha >> has_dormant;

#pragma omp parallel for
    for (int i = 0; i < n_cells; i++)
    {
        SPhase s_phase(std::stoi(arg_values[7]), std::stoi(arg_values[2]),
                       std::stoi(arg_values[3]), std::stoi(arg_values[4]),
                       std::stoi(arg_values[6]), has_dormant, data,
                       arg_values[1]);
        s_phase.simulate(i);
    }
}