# ReDyMo-CPP

![Build and Test](https://github.com/msreis/ReDyMo-CPP/actions/workflows/build-test-coverage.yaml/badge.svg)
[![codecov](https://codecov.io/gh/Dynamic-Systems-Biology-Group/ReDyMo-CPP/branch/master/graph/badge.svg?token=TtcbL2EvLU)](https://codecov.io/gh/Dynamic-Systems-Biology-Group/ReDyMo-CPP)

This is the C++ implementation of the ReDyMo model.
ReDyMo is a Python-coded, stochastic dynamic model simulator that reproduces the DNA replication process of cellular organisms belonging to the kinetoplastida group. Initially, we focused on Trypanosoma brucei strain TREU927.

## Requirements

In order to compile and run this model, one must have installed:

- CMake 3.9 or greater
- gcc compatible with the C++14 standard
- doxygen (Optional) for generating documentation

## Using docker

Before using this method, first make sure to have docker installed, the **right permissions** to use it, and internet connection to download the necessary docker images.

1. Run ReDyMo-CPP image in a container

```bash
docker run -it --rm -v <local_output_folder>:/opt/redymo/output --name <container name> brunobbs/redymo:tbrucei
```

Assuming all went successfully, the model runs inside a docker container and get the outputs in a folder locally.
This command runs the simulator with example params.

**IMPORTANT** The local_output_folder **must** have write permission to other users.

2. Run the image with custom parameters

```bash
docker run -it --rm -v <local output folder>:/opt/redymo/output --name <container name> brunobbs/redymo:tbrucei ./simulator <arguments to the model>

e.g.
docker run -it --rm -v "$(pwd)"/output:/opt/redymo/output --name redymo-runtime brunobbs/redymo:tbrucei ./simulator --cells 100 --organism 'Trypanosoma brucei brucei TREU927' --resources 10000 --speed 1 --period 100 --timeout 100000000 --dormant true --data-dir ./data
```

**IMPORTANT** in docker the _data_ folder is not in the default location, so it must be specified like in the example above. The _data_ folder is located at `/opt/redymo/data` and can be referenced as just `./data`.

## Compiling

After cloning this repository, in order to compile, one needs to create a new directory to separate compilation output files and source files.

1. Go to the ReDyMo-CPP folder

```bash
cd ReDyMo-CPP
```

2. Create a new folder to compile inside

```bash
mkdir build
```

3. Go to the **build** folder

```bash
cd build
```

4. Initialize CMake in the **build**

```bash
cmake .. [-D<OPTIONS> ...]
```

- if your **build** folder not inside the ReDyMo-CPP folder, the command is:

```bash
cmake <path to ReDyMo-CPP> [-D<OPTIONS> ...]
```

For a list of options see [Build Options](#Build-Options).

5. Compile ReDyMo-CPP

```bash
make
```

6. Make sure everything worked correctly by running the tests

```bash
make test
```

### Build Options

The custom options used by this project, aside form the CMake existing ones
are listed below. To use with the `cmake` command write: `-DOPTION=value`,
for example `-DCOVERAGE=OFF`.

|              Option |  Type   | Description                                                                                                                                             |
| ------------------: | :-----: | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
|     `BUILD_TESTING` | ON\|OFF | Whether to build the executables of unit tests.                                                                                                         |
| `CMAKE_BUILD_TYPE` | STRING  | Defines the compilation flags depending on the value. If no valid value is found, defaults to Release. See [Build Types](#Build-Types) for valid values. |
|          `COVERAGE` | ON\|OFF | Whether to configure the test coverage reports. **WARNING** This option overwrites the compilation flags.                                               |
|        `BUILD_GPGPU` | ON\|OFF | Whether to build the GPU code that depends on OpenCL.                                                                                                   |

#### Build Types

The available values for the `CMAKE_BUILD_TYPE` option are listed below. For
the exact flags set please refer to `CMakeLists.txt`.

| Build type    | Description                                                                                                                                     |
| -------------:| ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `Release`     | Default build type. Has warnings disabled and optimization level 3 (`-O3`).                                                                       |
| `Debug`       | Has all warnings enabled, debug flag enabled, profiling flag enabled and no optimizations.                                                      |
| `Coverage` | Has all warnings enabled, debug flag enabled and no optimizations. Enables coverage, if not already enabled.                                                                            |
| `Performance` | Disables coverage. Has all warnings enabled, optimization set to fast (`-Ofast`), native architecture optimization, and rename-registers enabled. |

### Run tests and coverage

To check whether the tests are passing use the command, after already configured CMake:

```bash
make test
```

To check the test coverage, use:

```bash
make coverage
<browser> build/coverage/index.html
```

## Parameters

In this version of ReDyMo, most parameters are mandatory and are listed below:

- **--cells** <number_of_cells>: Number of independent simulations to be made. number_of_cells is a positive integer.

- **--dormant** <dormant_flag>: Flag that, when present activates ('true') the firing of dormant origins and when its absent, disables ('false') it. It is noteworthy that the dormant origin firing does not work when constitutive origins are used (parameter --constitutive).

- **--organism** <'organism_name'>: Name of the parasite species, as saved in the database. 'organism_name' is a string (in space-separated names, use single quotation marks).

- **--resources** <number_of_forks>: Number of available forks for the replication process. number_of_forks is a positive integer.

- **--speed** <speed_value>: Velocity of each replication fork (in number of nucleotides per iteration). speed_value is a positive integer.

- **--timeout** <timeout_value>: Maximum allowed number of iterations of a simulation; if this value is reached, then a simulation is ended even if DNA replication is not completed yet.

The three optional parameters are:

- **--constitutive** <range>: When this parameter is provided, a DNA replication must use **only** the set of constitutive origins within the database instead of the probability landscape. `range` is a positive integer, and specifies the range of nucleotides around each constitutive origin that can initiate replication.

- **--period** <period_value>: Period (in number of simulation iterations) between two consecutive activations (i.e. RNAP binding) of a transcription region. period_value is a positive integer. If this parameter is not set, then the simulation is carried out without transcription.

- **--data-dir** <data_directory>: The directory containing the MFA-Seq_TBrucei_TREU927 folder for the organism and the database file. The database file must be named **database.sqlite**.

## Running the simulation

To run the program, the syntax of the main simulator program is the following one:

```bash
./simulator --cells number_of_cells --dormant dormant_flag --organism 'organism_name' --resources number_of_forks --speed speed_value --timeout timeout_value [--constitutive range] [--period period_value] [--data-dir directory_with_data]
```

The command above must be executed within the project main directory. For example, to run a simulation of 7 cells of _T. brucei TREU927_, with 10 forks, replisome speed of 65 bp/iteration, transcription period of 150 iterations between two transcription initiations, a timeout of one million iterations and with dormant origin firing, one must type at the project main directory:

```bash
./simulator --cells 7 --organism 'Trypanosoma brucei brucei TREU927' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true
```

The simulation results will be stored into a directory named _output/True_10_150/_, in which "output" is the outer directory name and the inner directory name of composed of the concatenation of the used parameter values for dormant origin firing, resources and period.

Another example: if one wants to simulate 30 cells of _T. brucei TREU927_, with 50 forks, replisome speed of 1 bp/iteration, no transcription, no dormant origin firing, using constitutive origins with a firing initiation range of 200 Kb, and the same timeout of the previous example:

```bash
./simulator --cells 30 --organism 'Trypanosoma brucei brucei TREU927' --resources 50 --speed 1 --timeout 1000000 --dormant false --constitutive 200000
```

In this case, the simulation results will be stored into the directory _output/False_50_0_.

## Generating documentation

This project has auto-generating DoxyGen documentation. In order to generate the documentation run inside the ReDyMo-CPP root folder:

```bash
doxygen
```

and then open the documentation index with a web browser:

```bash
<browser> docs/html/index.html
```

## Outputs

### Simulator

The output folder contains the outputs o each run in separate folders which names contain information about the parameters of the **run**.
If a simulation is run with the same parameters, the folder will be **overwritten**.
Inside each run folder are located each individual cell output which comprise of a **cell.txt** and the last state of every chromosome.
The **cell.txt** file contains, respectively, values for **number of forks**, **steps per iteration**, **total iterations**, **average inter origin distance**

### cell_output_aggregator.py

The cell_output_agregator.py script outputs a file with name equal to the **run** and txt format with the following attributes: **number of forks**, **steps per iteration**, **average number of iteration between all cells of the run**, **the standard deviation of iterations**, **average inter origin distance**, **standard deviation of IOD**, **number of cells in the run**

### collision_distance_median.py

Calculates **median**, **mean** and _std deviation_ of the distances between replication origins and conflicts with RNAP by processing the output files of each simulation.

### collision_distance_median_parallel.py

Does the same as collision_distance_median.py but in a parallel fashion.

## License

This program is distributed under the GPL v3 license, see the **LICENSE** file.

## Bug report and contact

If you have any bug report and/or want to contact for other subjects (e.g., to collaborate in this project), please do not hesitate to contact us!

Please, address your message to:

msreis at butantan dot gov dot br.
