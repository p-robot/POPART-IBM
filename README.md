Individual-based model for the HPTN071 (PopART) trial
=====================================================

Description
-----------

PopART-IBM is an individual-based model for simulating HIV epidemics in high-prevalence settings, as used in the HPTN 071 (PopART) trial.  Using the default parameters, the model simulates approximately 50,000 individuals over several decades following demographic patterns from UN's Population Division.

Documentation
-------------

- Model description: A full description of the model is described in [Pickles et al., 2020](https://www.medrxiv.org/content/10.1101/2020.08.24.20181180v1).  
- Code documentation: Detailed documentation of the model's C code is online [here](https://p-robot.github.io/POPART-IBM/) and included with this source code.
- Data dictionary of output files: A data dictionary describing all the output files from the model is [here](doc/output_files/output_file_dictionary_overview.md).
- Data dictionary of input parameters: A data dictionary describing the parameters in the model is described [here](doc/parameters/parameters.md).
 
Compilation
-----------

For Mac and Unix-type systems, PopART-IBM requires a C compiler (such as gcc) and the [GSL](https://www.gnu.org/software/gsl/) libraries installed:

```bash
cd POPART-IBM/src
make all
```

GSL can be downloaded from [here](ftp://ftp.gnu.org/gnu/gsl/).  


For Windows systems, please see [this](./doc/running_popartibm_on_windows.md) walkthrough.  

Usage
-----

```bash
cd POPART-IBM/src
./popart-simul.exe <inputdir> <nruns>
```
 
 where
 
* `inputdir`: Directory where input parameter files ("param_processed*.csv") are located
* `nruns` : number of simulation runs in parameter files (num. of lines in parameter files to read in)

A basic [example](examples/example_101.py) illustrates how the parameter input files can be set up (steps 1 and 2) as is expected by the model.  

**Notes**

* Additional command-line arguments are described in [main.c](src/main.c).  
* The model will write all output files to the directory `inputdir/Output` (additional command-line arguments can adjust this).  
* The output files written will depending upon which macros are set to 1 within the file [constants.h](src/constants.h) (those beginning `WRITE_*`).  

### Folder structure

```
src/                 # Model C code
   popart_ibm/       # Helper Python code for handling input parameters
doc/                 # General documentation in Markdown files
docs/                # Auto-generated code documentation using doxygen
examples/            # Basic script for running the model
tests/               # Testing files
python/              # Help Python scripts for manipulating markdown documentation
doxygen-awesome-css/ # Git submodule of CSS file for styling doxygen
```

Testing
-------

Tests are written [`pytest`](https://docs.pytest.org/en/stable/) using Python version 3.6 or later, and run in the following manner:

```bash
python3 -m pytest
```

Some tests take a long time to run and so are not run by default, they can be invoked using the ` --runslow` option in `pytest`:

```
python3 -m pytest --runslow
```

It is recommended that tests are run under a Python virtual environment.  The following will set up a Python virtual environment, install required modules, and run the tests: 

```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r tests/requirements.txt
python3 -m pytest
deactivate
```

Publications
------------

The model has been published in several journals and conferences:

- [Probert et al., 2019; CROI](https://www.croiconference.org/abstract/quantifying-transmissions-age-groups-simulations-hptn071-popart-trial/)
- [Probert et al., 2022; CROI](https://www.croiconference.org/abstract/hiv-1-dynamics-following-universal-testing-and-treatment-within-hptn-071-popart/)
- [Thomas et al., 2021](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X%2821%2900034-6/fulltext)
- [Pickles et al., 2022](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009301)
- [Probert et al., 2022](https://www.thelancet.com/journals/lanhiv/article/PIIS2352-3018(22)00259-4/fulltext)
- [Hall et al., 2024](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00220-3/fulltext)

Contributing
------------

Contributions are most welcome.  Please see the documentation on [contributing](CONTRIBUTING.md) for further information.  The C code is documented using formatting for [Doxygen](https://www.doxygen.nl/).  Please contact the core team or raise an issue before making a pull request.  
