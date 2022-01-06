# Science Data Challenge 2: Neutral Hydrogen

This repository contains a collection of scripts designed to generate mock neutral hydrogen (HI) datacubes for the second SKA science data challenge (SDC2).

Credit: Ian Harrison (original scripts developed for the [simuCLASS pipeline](https://bitbucket.org/itrharrison/simuclass/src/master/)).

### Prerequisites

The Python environment has been exported to `SDC2.yml` which lists the Python dependencies.

Also in use is fitsio. The fitsio version has been patched to allow arbitrarily large file sizes to be written directly to disk, and also to allow subsets of data to be written to the file on disk (this functionality is currenetly unavailable in AstroPy). 

To set up the environment for running the pipeline, first create the Python environment via conda

`conda env create --file SDC2.yml`

Activate this environment

`conda activate SDC2`

Then install fitsio by navigating to the patched source code directory e.g.

`cd /home/software/fitsio`

before installing via

`python setup.py install`





### Usage

`python run_pipeline.py inis/my_ini_file.ini` 


