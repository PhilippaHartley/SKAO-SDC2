import configparser
import os
import sys

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
    if config.getboolean("skymodel", "docontinuum"):
        from skymodel.skymodel_continuum import runSkyModel
        runSkyModel(config)
    if config.getboolean("skymodel", "doHI"):
        from skymodel.skymodel_HI import runSkyModel
        runSkyModel(config)

