import subprocess, os
import numpy as np
import math
import numpy.matlib
import multiprocessing as mp

cnresult_list = []
def cnlog_result(result):
    cnresult_list.append(result)

def run_fcont(config):
    #set up
    miriad_path = config.get("pipeline", "miriad_path")
    cver = config.get("pipeline", "cver")
    fcont_seed = int(config.getfloat("pipeline", "fcont_seed"))
    frames =  config.get("pipeline", "frames")
    bins =  config.get("pipeline", "bins")
    np.random.seed(fcont_seed)
    indir = config.get("pipeline", "indir")
    outdir = config.get("pipeline", "outdir")
    dfreq = 0.050
    sdfreq = '%.6e' % dfreq
    ds = 2.8
    db = 7.
    dbr = db/3600.*np.pi/180.
    sdbr = '%.6e' % dbr
    pool = mp.Pool()
    # convert from Jy/pix to Jy/bm, db arcs FWHM and ds arcs pixels
    px2bm = np.pi/(4.*math.log(2))*(db/ds)**2
    spx2bm = '%.6f' % px2bm

    if os.path.exists(outdir+'cube_frncont_'+cver): 
        print ('**** message from pipeline: '+outdir+'cube_frncont_'+cver+' already exists')
        print ('**** message from pipeline: not running func_fcont')
        return
    if not os.path.exists(indir+'sky_cont_'+cver): 
        print('converting fits files to miriad format')
        convert(indir, cver,miriad_path) 


    text = miriad_path+'maths exp="(<'+indir+'sky_cont_'+cver+'>)*'+spx2bm+'" out="'+outdir+'cube_fcont"'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in='+outdir+'cube_fcont/bunit value="JY/BEAM "'
    subprocess.run(text, shell='True')

    # resample continuum cube more finely
    text = miriad_path+'regrid in="'+outdir+'cube_fcont" axes=3 desc="0.95,1,0.01,21" tol=0.001 out="'+outdir+'cube_frcont_'+cver+'"'
    subprocess.run(text, shell='True')
    # make second imperfect version of the resampled continuum cube
    sigma = 0.001    
    pool = mp.Pool()
    for mfreq in range(950,1160,10):
        pool.apply_async(make_ncont, args=(mfreq,sigma,cver, outdir, frames, bins,miriad_path), callback = cnlog_result)
    pool.close()
    pool.join()
    print (cnresult_list)

    text = miriad_path+'imcat in="*_rnc" out="'+outdir+'cube_frncont_'+cver+'"'
    subprocess.run(text, shell='True')
    text = '/bin/rm -rf *_rnc '+outdir+'cube_fcont' 
    subprocess.run(text, shell='True')

def convert(indir, cver,miriad_path):
    text = miriad_path+'fits op="xyin" in="'+indir+'sky_continuum_f1_'+cver+'.fits" out="'+indir+'sky_cont_f1_'+cver+'"'
    subprocess.run(text, shell='True')
    text = miriad_path+'fits op="xyin" in="'+indir+'sky_continuum_f2_'+cver+'.fits" out="'+indir+'sky_cont_f2_'+cver+'"'
    subprocess.run(text, shell='True')
    text = miriad_path+'imcat in='+indir+'sky_cont_f1_'+cver+','+indir+'sky_cont_f2_'+cver+' out='+indir+'sky_cont_'+cver
    subprocess.run(text, shell='True')
    text = miriad_path+'fits in='+indir+'sky_continuum_f1_'+cver+'_z.fits op="xyin" out='+indir+'sky_cont_'+cver+'_z'
    subprocess.run(text, shell='True')





def make_ncont(mfreq,sigma,cver, outdir, frames, bins, miriad_path):
    ichan = 1 + (mfreq - 950)/10
    schan = str(ichan)
    smfreq = '%d' % mfreq
    sm4freq = smfreq.zfill(4)
    sifreq = '%.6f' % (mfreq/1000.)
    text = miriad_path+'imsub in="'+outdir+'cube_frcont_'+cver+'" out="'+sm4freq+'_rc" region="images('+schan+')"'
    subprocess.run(text, shell='True')
#    text = miriad_path+'imframe in="'+sm4freq+'_rc" out="'+sm4freq+'_rcf" frame="-2944,2943,-2944,2943"'
#    text = miriad_path+'imframe in="'+sm4freq+'_rc" out="'+sm4freq+'_rcf" frame="-656,655,-656,655"'
    text = miriad_path+'imframe in="'+sm4freq+'_rc" out="'+sm4freq+'_rcf" frame="'+frames+'"'
    subprocess.run(text, shell='True')
    text = miriad_path+'imbin in="'+sm4freq+'_rcf" out="'+sm4freq+'_rcb" bin="'+bins+'"'
    subprocess.run(text, shell='True')
    text = miriad_path+'imgen in="'+sm4freq+'_rcb" out="'+sm4freq+'_rcbn" object="noise" spar=1. seed="'+smfreq+'"'
    subprocess.run(text, shell='True')
    text = miriad_path+'regrid in="'+sm4freq+'_rcbn" out="'+sm4freq+'_rcbnx" tin='+sm4freq+'_rc""'
    subprocess.run(text, shell='True')
    ssigma = '%.6e' % sigma
    text = miriad_path+'maths exp="(<'+sm4freq+'_rc>*(1.+'+ssigma+'*<'+sm4freq+'_rcbnx>))" out="'+sm4freq+'_rnc"'
    subprocess.run(text, shell='True')
    # clean up intermediate products
    text = '/bin/rm -rf '+sm4freq+'_rc '+sm4freq+'_rcf '+sm4freq+'_rcb '+sm4freq+'_rcbn '+sm4freq+'_rcbnx' 
    subprocess.run(text, shell='True')
    text = 'Finished NCont Frequency: '+sm4freq
    return text


if __name__ == '__main__':
    run_fcont(config)


