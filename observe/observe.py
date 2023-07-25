import configparser
import os
import sys



def run_observe(config):
    if not os.path.exists('out/observe'): 
        os.system('mkdir out/observe')
    from observe.func_fcont import run_fcont
    run_fcont(config)

    from observe.func_fhia import run_fhia
    run_fhia(config)

    from observe.proc_cubs import run_cubs
    run_cubs(config)
    from observe.func_fcontIm import run_fcontIm
    run_fcontIm(config)




if __name__ == '__main__':
    run_observe(config)



