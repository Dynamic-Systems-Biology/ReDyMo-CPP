# ReDyMo-CPP
This is the C++ implementation of the ReDyMo model.
ReDyMo is a Python-coded, stochastic dynamic model simulator that reproduces the DNA replication process of cellular organisms belonging to the kinetoplastida group. Initially, we focused on Trypanosoma brucei strain TREU927.

## Requirements
In order to compile and run this model, one must have installed:

- CMake 3.0 or greater
- gcc compatible with the C++14 standard
- doxygen (Optional) for generating documentation

## Compiling
After cloning this repository, in order to compile, one needs to create a new directory to separate compilation output files and source files.
1. Go to the ReDyMo-CPP folder
```
cd ReDyMo-CPP
```
2. Create a new folder to compile inside
```
mkdir build
```
3. Go to the __build__ folder
```
cd build
```
4. Initialize CMake in the __build__
```
cmake ..
```
- if your __build__ folder not inside the ReDyMo-CPP folder, the command is:
```
cmake <path to ReDyMo-CPP>
```
5. Compile ReDyMo-CPP
```
make
```
6. Make shure everything worked correctly by running the tests
```
make test
```


## Parameters
In this version of ReDyMo, most parameters are mandatory and are listed below:

- --cells number_of_cells: Number of independent simulations to be made. number_of_cells is a positive integer.

- --dormant dormant_flag: Flag that either activates ('true') or disables ('false') the firing of dormant origins. dormant_value is a Boolean flag. It is noteworthy that the dormant origing firing does not work when constitutive origins are used (parameter --constitutive).

- --organism 'organism_name': Name of the parasite species, as saved in the database. 'organism_name' is a string (in space-separated names, use single quotation marks).

- --resources number_of_forks: Number of available forks for the replication process. number_of_forks is a positive integer.

- --speed speed_value: Velocity of each replication fork (in number of nucleotides per iteration). speed_value is a positive integer.

- --timeout timeout_value: Maximum allowed number of iterations of a simulation; if this value is reached, then a simulation is ended even if DNA replication is not completed yet.

The two optional parameters are:

- --constitutive range: When this parameter is provided, a DNA replication must use the set of constitutive origins within the database instead of the probability landscape. range is a positive integer, and specifies the range of nucleotides around each constitutive origin that can initiate replication.

- --period period_value: Period (in number of simulation iterations) between two consecutive activations (i.e. RNAP binding) of a transcription region. period_value is a positive integer. If this parameter is not set, then the simulation is carried out without transcription.

## Running the simulation

To run the program, the syntax of the main simulator program is the following one:
```
$ ./simulator --cells number_of_cells --dormant dormant_flag --organism 'organism_name' --resources number_of_forks --speed speed_value --timeout timeout_value [--constitutive range] [--period period_value]
```

The command above must be executed within the project main directory. For example, to run a simulation of 7 cells of *T. brucei TREU927*, with 10 forks, replisome speed of 65 bp/iteration, transcription period of 150 iterations between two transcription initiations, a timeout of one million iterations and with dormant origin firing, one must type at the project main directory:
```
$ ./simulator --cells 7 --organism 'Trypanosoma brucei brucei TREU927' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true
```
The simulation results will be stored into a directory named *output/True_10_150/*, in which "output" is the outer directory name and the inner directory name of composed of the concatenation of the used parameter values for dormant origin firing, resources and period.

Another example: if one wants to simulate 30 cells of *T. brucei TREU927*, with 50 forks, replisome speed of 1 bp/iteration, no transcription, no dormant origin firing, using constitutive origins with a firing initiation range of 200 Kb, and the same timeout of the previous example:
```
$ ./simulator --cells 30 --organism 'Trypanosoma brucei brucei TREU927' --resources 50 --speed 1 --timeout 1000000 --dormant false --constitutive 200000
```

In this case, the simulation results will be stored into the directory *output/False_50_0*.

## Generating documentation
This project has auto-generating DoxyGen documentation. In order to generate the documentation run inside the ReDyMo-CPP root folder:
```
doxygen
```
and then open the documentation index with a web browser:
```
<browser> docs/html/index.html
```

## License
This program is distributed under the GPL v3 license, see the __LICENSE__ file.

## Bug report and contact

If you have any bug report and/or want to contact for other subjects (e.g., to collaborate in this project), please do not hesitate to contact us!

Please, address your message to:

msreis at butantan dot gov dot br.