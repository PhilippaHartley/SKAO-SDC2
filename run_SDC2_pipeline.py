import configparser
import os
import sys
import time

doHI = True
docontinuum = True
doobserve = True

tstart = time.time()

cver = 'dev'

if doHI:
    config = configparser.ConfigParser()
    config.read('inis/skymodel/HI/SDC2_HI_'+cver+'.ini')
    from skymodel.skymodel_HI import runSkyModel
    runSkyModel(config)
if docontinuum:
    config = configparser.ConfigParser()
    config.read('inis/skymodel/continuum/SDC2_continuum_'+cver+'_f1.ini')
    from skymodel.skymodel_continuum import runSkyModel
    runSkyModel(config)
    config = configparser.ConfigParser()
    config.read('inis/skymodel/continuum/SDC2_continuum_'+cver+'_f2.ini')
    runSkyModel(config)
if doobserve:
    config = configparser.ConfigParser()
    config.read('inis/observe/SDC2_observe_'+cver+'.ini')
    from observe.observe import run_observe
    run_observe(config)

tend = time.time()
print("pipeline finished in {0} seconds.".format(tend - tstart))


