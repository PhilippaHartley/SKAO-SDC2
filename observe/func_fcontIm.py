import subprocess, os
import numpy as np
import math
import numpy.matlib
import multiprocessing as mp

result_list = []
def log_result(result):
    result_list.append(result)

def run_fcontIm(config):
    miriad_path = config.get("pipeline", "miriad_path")
    cver = config.get("pipeline", "cver")
    outdir = config.get("pipeline", "outdir")
    prdir = outdir+'cont_'+cver+'_fpr/'
    prndir = outdir+'cont_'+cver+'_fprn/'
    if os.path.exists(outdir+'cont_'+cver+'.fits'):
        print ('**** message from pipeline: '+outdir+'cube_'+cver+'.fits already exists')
        print ('**** message from pipeline: not running func_fcontIm_')
        exit()
    if not os.path.exists(prdir):
        text = '/bin/mkdir '+prdir 
        subprocess.run(text, shell='True')
    if not os.path.exists(prndir):
        text = '/bin/mkdir '+prndir 
        subprocess.run(text, shell='True')

    # Pixel increments in frequency (GHz) and spatially (arcsec), npix and ref pix
    dfreq = 0.010
    ds = 2.8
    db = 7.
    npix = int(config.get("pipeline", "npix"))
    crp12 = int(config.get("pipeline", "crp12"))
    
    # Starting frequency of continuum (sub-)cube
    freq1 = 0.95
    freq2 = 1.15
    nfreq = 21
    rfreq_int = int(config.get("pipeline", "rfreq_int2"))

    pool = mp.Pool()
    for ichan in range(1,22,1):
        pool.apply_async(cont_proc, args=(ichan,freq1,dfreq,nfreq,rfreq_int,ds,db,npix,crp12,prdir,prndir,cver,outdir,miriad_path), callback = log_result)

    pool.close()
    pool.join()
    print (result_list)
    dbrad = db*np.pi/(3600.*180.)
    sdbrad = '%.6e' % dbrad
    text = miriad_path+'imcat in="'+prdir+'*.m" out="'+prdir+'cube_c"'
    subprocess.run(text, shell='True')
    text = miriad_path+'imcat in="'+prndir+'*.m" out="'+prndir+'cube_cn"'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="'+prndir+'cube_cn/bmaj" value='+sdbrad+''
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="'+prndir+'cube_cn/bmin" value='+sdbrad+''
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="'+prndir+'cube_cn/bpa" value=0.'
    subprocess.run(text, shell='True')

    text = miriad_path+'fits op=xyout in="'+prndir+'cube_cn" out="'+outdir+'cont_'+cver+'.fits"'
    subprocess.run(text, shell='True')

def cont_proc(ichan,freq1,dfreq,nfreq,rfreq_int,ds,db,npix,crp12,prdir,prndir,cver,outdir,miriad_path):
    sdfreq = '%.6e' % dfreq
    sichan = '%d' % ichan
    mfreq = (freq1 + (ichan-1)*dfreq)*1000.
    smfreq = '%d' % mfreq
    sm4freq = smfreq.zfill(4)
    ifreq = (mfreq/1000.)
    sifreq = '%.6f' % (mfreq/1000.)
    # make modified version of freq to diversify random seeds between runs
    rfreq = ifreq+rfreq_int
    srfreq = '%.6f' % rfreq

    # make a PB for the channel at hand
    text = '/bin/cp -r telescope/PB_I_14 pb'+sm4freq+'.b'
    subprocess.run(text, shell='True')
    cdel = 5.646e-4*1.4/ifreq
    scdel = '%.6e' % cdel
    text = miriad_path+'puthd in=pb'+sm4freq+'.b/cdelt1 value=-'+scdel
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=pb'+sm4freq+'.b/cdelt2 value='+scdel
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=pb'+sm4freq+'.b/crval3 value='+sifreq
    subprocess.run(text, shell='True')

    # read input sky model slice
    text = miriad_path+'imsub in="'+outdir+'cube_frncont_'+cver+'" region="images('+sichan+')" out="'+sm4freq+'_c"'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in='+sm4freq+'_c/crval3 value='+sifreq
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in='+sm4freq+'_c/cdelt3 value='+sdfreq
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in='+sm4freq+'_c/crpix3 value=1'
    subprocess.run(text, shell='True')

    # simulate visibilities at this frequency
    text = miriad_path+'uvgen source="telescope/fonaxis.src" ant="telescope/b2v1ska1mid.neu.197.ant" baseunit=-3.33564 corr="10,1,0,10." time="20JAN1'+srfreq+'" freq="'+sifreq+',0" radec="12.,-30." harange="-4.,4.,0.01667" stokes="rr" lat=-31. pbfwhm=3300. center="0,-0" systemp=30. jyperk=20. out="c'+sm4freq+'.uv" telescop="altaz,0"'
    subprocess.run(text, shell='True')
    # make dirty beam and noise images
    snpix = '%i' % npix
    sds = '%.6f' % ds
    sdb = '%.6f' % db
    scrp12 = '%.6f' % crp12
    text = miriad_path+'invert vis="c'+sm4freq+'.uv" map="c'+sm4freq+'.m" beam="c'+sm4freq+'.b" line="channel,1,1,10,1" imsize='+snpix+' cell='+sds+' options="mfs,radial,double" sup="64,64,0.5" fwhm='+sdb+''
    subprocess.run(text, shell='True')
    text = miriad_path+'velsw in="c'+sm4freq+'.b" axis="FREQ"'
    subprocess.run(text, shell='True')
    text = miriad_path+'velsw in="c'+sm4freq+'.m" axis="FREQ"'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=c'+sm4freq+'.m/cdelt3 value='+sdfreq+''
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=c'+sm4freq+'.m/ctype3 value="FREQ    "'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=c'+sm4freq+'.m/naxis value=3'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=c'+sm4freq+'.m/crpix1 value='+scrp12+''
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in=c'+sm4freq+'.m/crpix2 value='+scrp12+''
    subprocess.run(text, shell='True')
    # the input noise images have an RMS that is extracted by channel
    # the desired noise levels have RMS that varies linearly with frequency (950 to 1350 MHz) as 13.6 to 17.1 microJy/Beam (for inc of 115 kHz)
    # the desired noise levels have RMS that varies linearly with frequency (950 to 1350 MHz) (for other inc divide by sqrt(bandwidth ratio))
    # also take account of the confusion noise following Condon et al 2012
    text = miriad_path+'histo in="c'+sm4freq+'.m" > histc'+sm4freq+'.txt'
    subprocess.run(text, shell='True')
    with open('histc'+sm4freq+'.txt') as f:
        content = f.readlines()
                
    irms = float(content[5][26:37])
    trms = (13.6e-6+3.5e-6*(ifreq-0.95)/(1.35-0.95))/(dfreq*1.e9/115.e3)**0.5
    crms = 1.2e-6*(ifreq/3.02)**-0.7*(db/8.)**3.33
    orms = (trms**2 + crms**2)**0.5
    norm = trms/irms
    snorm = '%.8f' % norm
    text = miriad_path+'maths exp="<c'+sm4freq+'.m>*'+snorm+'" out="n'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # expand the PB to full image size
    text = miriad_path+'regrid in="pb'+sm4freq+'.b" tin="c'+sm4freq+'.b" axes="1,2" out="pbr'+sm4freq+'.b"'
    subprocess.run(text, shell='True')
    # modify dirty beam for mosaic situation (multiply by PB)
    text = miriad_path+'maths exp="<c'+sm4freq+'.b>*<pbr'+sm4freq+'.b>" out="mos'+sm4freq+'.b"'
    subprocess.run(text, shell='True')
    norm = 4.0*orms
    sn = '%.8f' % norm

    # extract peaks (> 4 effective RMS including confusion) from smooth (5 arcsec) sky image and reset blanks to zero
    text = miriad_path+'maths exp="<'+sm4freq+'_c>-'+sn+'" mask="<'+sm4freq+'_c>.gt.'+sn+'" out="pm'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    text = miriad_path+'imblr in="pm'+sm4freq+'.m" out="p'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # generate residuals
    text = miriad_path+'maths exp="<'+sm4freq+'_c>-<p'+sm4freq+'.m>" out="r'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # deconvolve beam from sky residuals and update header for "Jy/pixel" units
    text = miriad_path+'convol map="r'+sm4freq+'.m" fwhm=5.0 options="divide" sigma=1.e-1 out="rd'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="rd'+sm4freq+'.m/bmaj" value=0.0'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="rd'+sm4freq+'.m/bmin" value=0.0'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="rd'+sm4freq+'.m/bpa" value=0.0'
    subprocess.run(text, shell='True')
    text = miriad_path+'puthd in="rd'+sm4freq+'.m/bunit" value="JY/PIXEL"'
    subprocess.run(text, shell='True')
    # convolve residuals with mosaic dirty beam
    text = miriad_path+'convol map="rd'+sm4freq+'.m" beam="mos'+sm4freq+'.b" out="rdc'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # add back the smooth peaks
    text = miriad_path+'imcomb in="p'+sm4freq+'.m,rdc'+sm4freq+'.m" options="nonorm" out="'+prdir+'pr'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # add in the noise image
    text = miriad_path+'imcomb in="'+prdir+'pr'+sm4freq+'.m,n'+sm4freq+'.m" options="nonormalise,relax,fqaver" out="'+prndir+'prn'+sm4freq+'.m"'
    subprocess.run(text, shell='True')
    # clean up intermediate products
    text = '/bin/rm -rf '+'pb'+sm4freq+'.b c'+sm4freq+'.uv c'+sm4freq+'.b c'+sm4freq+'.m  c.log n'+sm4freq+'.m pbr'+sm4freq+'.b '+sm4freq+'_c  '
    subprocess.run(text, shell='True')
    text = '/bin/rm -rf '+'mos'+sm4freq+'.b pm'+sm4freq+'.m p'+sm4freq+'.m r'+sm4freq+'.m rd'+sm4freq+'.m rdc'+sm4freq+'.m histc'+sm4freq+'.txt'
    subprocess.run(text, shell='True')

    text = 'Finished channel: '+sichan
    return text

if __name__ == '__main__':
    run_fcontIm(config)

