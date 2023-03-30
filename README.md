# Science Data Challenge 2: Neutral Hydrogen

This repository contains a collection of scripts used to generate mock SKA-Mid-observed neutral hydrogen (HI) datacubes for the second SKA Science Data Challenge ([SDC2](https://sdc2.astronomers.skatelescope.org/)).

### Prerequisites

#### Environment

The Python environment has been exported to `SDC2_full_dependencies.yml` which lists the Python dependencies.

Also in use is a version of fitsio that has been patched to allow arbitrarily large file sizes to be written directly to disk, and to allow subsets of data to be written to the file on disk (this functionality is currently unavailable in AstroPy). The patched version can be downloaded from [here](https://drive.google.com/drive/folders/15h0hE-cnqvS6xpX90qtX_Ji1wzC65V9R?usp=sharing). 

To set up the environment for running the pipeline, first create the Python environment via conda

`conda env create --file SDC2.yml`

Activate this environment

`conda activate SDC2`

Then install fitsio by navigating to the patched source code directory e.g.

`cd /home/software/fitsio`

before installing via

`python setup.py install`

#### Source input catalogues

The pipeline uses catalogue files produced using [T-RECS](https://github.com/abonaldi/TRECS). The catalogues produced for use during SDC2 can be downloaded from [here](https://drive.google.com/drive/folders/15h0hE-cnqvS6xpX90qtX_Ji1wzC65V9R?usp=sharing). For use with these scripts, place the SDC2_catalogues directory inside the directory that contains this repository.

### Basic usage

`python run_SDC2_pipeline.py`


### Detailed usage

#### How to run the SDC2 simulation pipeline

An end-to-end pipeline, `run_SDC2_pipeline.py`, can be used to simulate the data products produced for SDC2. 

Four data product versions were produced for SDC2:
'dev' (FoV 30 arcmin): a 'development' datacube used by SDC2 teams for pipeline development
'ldev' (FoV 60 arcmin): the same as 'dev' but covering a larger field of view
'eval' (FoV 30 arcmin): an 'evaluation' data cube used by SDC2 teams for validating subnmission file formats
'full' (FoV 273 arcmin): the full Challenge datacube

The default version that will be run is 'dev'. This can be changed by passing a version name as an argument to the main script, e.g. `python run_SDC2_pipeline.py eval`.

There are two main stages to the pipeline: `skymodel` takes a catalogue of sources to produce image models of the sky;  `observe` takes the image models to produce the sky as observed by the SKA-Mid telescope. Each stage uses several modules, following the recipe below.

Recipe for producing an SKA-Mid-observed HI sky

1. Run `skymodel_HI.py` over frequency range 950-1150 MHz
This module takes as input a T-RECS catalogue containing a list of simulated sources and their HI and continuum properties. The module converts the morphological HI emission properties of each source into a 'postage stamp' image cube. Each cube is then placed in the full HI emission field at its catalogued position. This module also outputs a new catalogue file that lists the subset of TRECS sources used in the simulation with their HI properties.

2. Run `skymodel_continuum.py` over frequency range 950-1150 MHz
This module uses the same T-RECS catalogue to produce a corresponding field of continuum emission at the same frequency as the HI field. 

3. Run `skymodel_continuum.py` over frequency range 1200-1400 MHz
This module uses the same T-RECS catalogue but produces a continuum field at a higher frequency. 

4. Run `observe.py`
This module takes as input the three skymodel outputs and uses them to create an SKA-Mid-observed model of the HI sky.  The lower frequency continuum cube is used to simulate imperfect continuum subtraction. The higher frequency continuum field is used to create HI absorption signatures. The module outputs the final HI image cube and a corresponding continuum image cube.
 
The `inis` directory contains the initialisation files that are used for each step. Each ini file is available in either a 'dev', ldev', 'eval', or 'full' variation, representing the 'development', large development', evalutaion', and 'full Challenge' datasets produced for SDC2.

At the end of the run, the SKA-Mid-observed HI and continuum image cubes can be found along with the output HI source catalogue in `out/products`.

#### Simulation description

A detailed description of the simulations is available in Section 3 of the SDC2 paper [here](https://arxiv.org/abs/2303.07943).




