import pdb
import numpy as np
import sys
import os
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

if not os.path.exists(config.get('pipeline', 'base_dir')+ config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')):
  os.mkdir(config.get('pipeline', 'base_dir')+ config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name'))

if config.getboolean('pipeline', 'doskymodel'):
  from skymodel.skymodel import runSkyModel
  runSkyModel(config)

# put diagnostics option in here  
if config.getboolean('pipeline', 'dodiagnostics'):
  from skymodel.diagnostics import make_plots
  make_plots(config)



if config.getboolean('pipeline', 'dosimdata'):
  from simulatedata.simdata import runSimulateData
  runSimulateData(config)#, sys.argv[1])

if config.getboolean('pipeline', 'doimagedata'):
  from imager.imager import runNWImager, runCASAClean, runWSClean
  if config.get('imager', 'type') == 'casaclean':
    runCASAClean(config, sys.argv[1])
  elif config.get('imager', 'type') == 'nwimager':
    runNWImager(config, sys.argv[1])
  elif config.get('imager', 'type') == 'wsclean':
    runWSClean(config)
    if config.getboolean('imager', 'dopostagestamps'):
      from thumbnailer.thumbnailer import makeThumbnails
      makeThumbnails(config)
  else:
    'You picked an unknown imager!'
