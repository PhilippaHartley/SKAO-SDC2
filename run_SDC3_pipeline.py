import configparser
import os
import sys
import time

doHI = False
docontinuum = True
doobserve = False

tstart = time.time()
#print (sys.argv)
#if len(sys.argv)>1:
#   cver = sys.argv[1]
#   if not cver in ['dev', 'ldev', 'eval', 'full']:
#       print ("please select from: dev, ldev, eval, full")
#       exit()
#else:       
#    cver = 'dev'




if doHI:
    config = configparser.ConfigParser()
    config.read('inis/skymodel/HI/SDC2_HI_'+cver+'.ini')
    from skymodel.skymodel_HI import runSkyModel
    runSkyModel(config)
if docontinuum:
    config = configparser.ConfigParser()
    config.read('inis/SDC3/SDC3_continuum_v2.ini')
    from skymodel.skymodel_continuum import runSkyModel
    runSkyModel(config)
#    config = configparser.ConfigParser()
#    config.read('inis/skymodel/continuum/SDC2_continuum_'+cver+'_f2.ini')
#    runSkyModel(config)
if doobserve:
    config = configparser.ConfigParser()
    config.read('inis/observe/SDC2_observe_'+cver+'.ini')
    from observe.observe import run_observe
    run_observe(config)

tend = time.time()
print("pipeline finished in {0} seconds.".format(tend - tstart))


