import configparser
import os
import pdb
import sys

import numpy as np

config = configparser.ConfigParser()
config.read(sys.argv[1])

if not os.path.exists(
    config.get("pipeline", "base_dir")
    + config.get("pipeline", "data_path")
    + config.get("pipeline", "project_name")
):
    os.mkdir(
        config.get("pipeline", "base_dir")
        + config.get("pipeline", "data_path")
        + config.get("pipeline", "project_name")
    )

if config.getboolean("pipeline", "doskymodel"):
    from skymodel.skymodel import runSkyModel

    runSkyModel(config)

# put diagnostics option in here
if config.getboolean("pipeline", "dodiagnostics"):
    from skymodel.diagnostics import make_plots

    make_plots(config)

if config.getboolean("pipeline", "dotelescope"):
    from simulatedata.simdata import runSimulateData

    runTelescope(config)
