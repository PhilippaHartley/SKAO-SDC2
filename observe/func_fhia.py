import subprocess
import numpy as np
import math
import numpy.matlib
import multiprocessing as mp
import os

result_list = []
def log_result(result):
    result_list.append(result)

def run_fhia(config):
    pool = mp.Pool()
    miriad_path = config.get("pipeline", "miriad_path")
    fhia_seed = int(config.getfloat("pipeline", "fhia_seed"))
    np.random.seed(fhia_seed)
    indir = config.get("pipeline", "indir")
    outdir = config.get("pipeline", "outdir")
    cver = config.get("pipeline", "cver")
    prdir = outdir+cver+'_fpr/'
    prndir = outdir+cver+'_fprn/'
    taudir = outdir+cver+'_tau/'

    if not os.path.exists(indir+'sky_HI_'+cver):       
        convert(indir, cver,miriad_path)

    if not os.path.exists(prdir):
        text = '/bin/mkdir '+prdir 
        subprocess.run(text, shell='True')
    if not os.path.exists(prndir):
        text = '/bin/mkdir '+prndir 
        subprocess.run(text, shell='True')
    if not os.path.exists(taudir):
        text = '/bin/mkdir '+taudir 
        subprocess.run(text, shell='True')
 

    # Read file with assumed RFI modulation of SEFD
    frfi = 'telescope/mktpdat.txt'
    mkdata = np.loadtxt(frfi)
    mkfrq = np.array([row[0] for row in mkdata])
    mkrat = np.array([row[3] for row in mkdata])
    mkfrqr = mkfrq[::-1]*1.e-3
    mkratr = mkrat[::-1]
    
    # Pixel increments in frequency (GHz) and spatially (arcsec), npix and ref pix
    dfreq = 0.000030
    ds = 2.8
    db = 7.
    npix = int(config.get("pipeline", "npix"))
    crp12 = int(config.get("pipeline", "crp12"))
    
    # Starting frequency of HI (sub-)cube
    freq1 = 0.95
    freq2 = 1.15
    nfreq = 6668
    rfreq_int = int(config.get("pipeline", "rfreq_int1"))


    lfreq = np.linspace(freq1,freq2,nfreq)
    lrat = np.interp(lfreq,mkfrqr,mkratr)

    # Nominal PB at 1.4 GHz is called PB_I_14 and is constructed as below
    # maths exp="(<E_HH14f>+<E_VV14f>)/2." out=PB_I_14
    # puthd in=PB_I_14/ctype1 value='RA---SIN'
    # puthd in=PB_I_14/ctype2 value='DEC--SIN'
    # puthd in=PB_I_14/crval1 value=0.
    # puthd in=PB_I_14/crval2 value=-0.5235988
    # puthd in=PB_I_14/cdelt1 value=-5.646e-4
    # puthd in=PB_I_14/cdelt2 value=5.646e-4

    pool = mp.Pool()
    for ichan in range(1,6669,1):#(1,6669,1):
        lrati = lrat[ichan - 1]
        pool.apply_async(chan_proc, args=(ichan,freq1,dfreq,nfreq,rfreq_int,ds,db,npix,crp12,lrati,prdir,prndir,taudir,cver,indir, outdir, miriad_path), callback = log_result)

    pool.close()
    pool.join()
    print ('**** message from pipeline:', result_list)


def convert(indir, cver,miriad_path):


    text = miriad_path+'fits op=xyin in='+indir+'sky_HI_'+cver+'.fits out='+indir+'sky_HI_'+cver
    subprocess.run(text, shell='True')
    text = miriad_path+'velsw in='+indir+'sky_HI_'+cver+' axis=FREQ'
    subprocess.run(text, shell='True')

def chan_proc(ichan,freq1,dfreq,nfreq,rfreq_int,ds,db,npix,crp12,lrati,prdir,prndir,taudir,cver,indir, outdir,miriad_path):
    # HI rest freq and basic cosmology parms
    hi0freq = 1.42040575
    H0 = 67.3
    Om = 0.32
    ckms = 3.e5
    Ms = 1.99e33
    amu = 1.66e-24
    Mpc = 3.09e24
    sdfreq = '%.6e' % dfreq
    dsrad = ds*np.pi/(3600.*180.)
    px2bm = np.pi/(4.*math.log(2))*(db/ds)**2
    Hz2ch = dfreq*1.e9
    Cfact = 49.8*Hz2ch*Ms/(px2bm*amu*dsrad**2*Mpc**2)
    Cfactp = 49.8*Hz2ch*Ms/(amu*dsrad**2*Mpc**2)
    # set levels at which RFI will give some/continuous simulated flagging
    lratmin = 1.1
    lratmax = 2.
    # set uvmax (in kilolambda) where low level RFI extends to, uvmax scales non-linearly with lrati
    # first guess is 15m at 21cm
    uvmaxf = (15./0.21)/1000.

    invchan = nfreq + 1 - ichan
    ifreq = freq1 + (ichan-1)*dfreq
    wlcm = ckms*1.e5/(ifreq*1.e9)
    ifreqn = freq1 + (ichan-1-1)*dfreq
    ivel = ckms*(hi0freq-ifreq)/ifreq
    iveln = ckms*(hi0freq-ifreqn)/ifreqn
    dvel = iveln-ivel
    sdvel = '%.6f' % dvel
    redz = (hi0freq/ifreq - 1.)
    sredz = '%.6f' % redz
    Dlum = (1.+redz)*(redz*ckms/H0)/(1.+1.718*Om*redz+0.315*Om**0.5*redz**2)**0.5
    Dang = Dlum/(1.+redz)**2
    sifreq = '%.6f' % ifreq
    # make modified version of freq to diversify random seeds between runs
    rfreq = ifreq+rfreq_int
    srfreq = '%.6f' % rfreq
    sichan = str(ichan).zfill(4)
    sinvchan = str(invchan)
    # first check whether outfput file exists, in which case skip this channel
    outfile = prndir+'prn'+sichan+'.m'
    #outfile = taudir+'hia'+sichan+'.m'
    if os.path.exists(outfile):
        text = 'Previously finished channel: '+sichan
        return text
        
    else:
        # make a PB for the channel at hand
        text = '/bin/cp -r telescope/PB_I_14 pb'+sichan+'.b'
        subprocess.run(text, shell='True')
        cdel = 5.646e-4*1.4/ifreq
        scdel = '%.6e' % cdel
        text = miriad_path+'puthd in=pb'+sichan+'.b/cdelt1 value=-'+scdel
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=pb'+sichan+'.b/cdelt2 value='+scdel
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=pb'+sichan+'.b/crval3 value='+sifreq
        subprocess.run(text, shell='True')
        # simulate visibilities at this frequency (timestamp is used as random number seed so it is modified for each new frequency)
        text = miriad_path+'uvgen source="telescope/fonaxis.src" ant="telescope/b2v1ska1mid.neu.197.ant" baseunit=-3.33564 corr="1,1,0,0.030" time="20JAN1'+srfreq+'" freq="'+sifreq+',0" radec="12.,-30." harange="-4.,4.,0.01667" stokes="rr" lat=-31. pbfwhm=3300. center="0,-0" systemp=30. jyperk=20. out="c'+sichan+'.uv" telescop="altaz,0"'
        subprocess.run(text, shell='True')
        # simulate RFI flagging at this frequency using lrati, the SEFD multiplier for this channel
        # the uv range to be flagged extends from 0 to values larger than uvmaxf scaling nonlinearly with lrati
        if (lrati < 1.) :
            lrati = 1.
            
        uvr1 = 0.
        uvr2 = uvmaxf*10**((lrati-1.)**0.333)
        # when multipler is larger than lratmax, then RFI is continuous
        if (lrati > lratmax) :
            ha1 = -4.
            ha2 = 4.
            # when multipler is in range lratmin to lratmax then a fraction of the ha range is impacted 
        elif (lrati > lratmin) and (lrati <= lratmax):
            delha = ((lrati - lratmin)/(lratmax - lratmin))*8.
            ha1 = -4. + (8.-delha)*np.matlib.rand(1)
            ha2 = ha1 + delha
            # when multipler is less than lratmin then none of the ha range is impacted 
        else:
            ha1 = -4.
            ha2 = -3.999

        sha1 = '%.6f' % ha1
        sha2 = '%.6f' % ha2
        suvr1 = '%.6f' % uvr1
        suvr2 = '%.6f' % uvr2
        text = miriad_path+'uvdflag vis="c'+sichan+'.uv" select="ha('+sha1+','+sha2+')" uvrange="'+suvr1+','+suvr2+'"'
        subprocess.run(text, shell='True')
        # make dirty beam and noise images
        snpix = '%i' % npix
        sds = '%.6f' % ds
        sdb = '%.6f' % db
        dba = db/(2.)**0.5
        dbar = (dba/3600.)*np.pi/180.
        sdbar = '%.6e' % dbar
        crp12m = crp12 - 1
        scrp12 = '%.6f' % crp12
        scrp12m = '%.6f' % crp12m
        text = miriad_path+'invert vis="c'+sichan+'.uv" map="c'+sichan+'.m" beam="c'+sichan+'.b" line="channel,1,1,1,1" imsize='+snpix+' cell='+sds+' options="radial,double" sup="64,64,0.5" fwhm='+sdb
        subprocess.run(text, shell='True')
        text = miriad_path+'velsw in="c'+sichan+'.b" axis="FREQ"'
        subprocess.run(text, shell='True')
        text = miriad_path+'velsw in="c'+sichan+'.m" axis="FREQ"'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=c'+sichan+'.m/cdelt3 value='+sdfreq
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=c'+sichan+'.m/ctype3 value="FREQ    "'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=c'+sichan+'.m/naxis value=3'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=c'+sichan+'.m/crpix1 value='+scrp12
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=c'+sichan+'.m/crpix2 value='+scrp12
        subprocess.run(text, shell='True')
        # the input noise images have an RMS that is extracted by channel
        # the desired noise levels have RMS that varies linearly with frequency (950 to 1350 MHz) as 13.6 to 17.1 microJy/Beam (at 115 kHz sampling)
        # the nominal noise level is scaled with the MeerKAT total power spectrum to account for RFI
        text = miriad_path+'histo in="c'+sichan+'.m" > histr'+sichan+'.txt'
        subprocess.run(text, shell='True')
        with open('histr'+sichan+'.txt') as f:
            content = f.readlines()
                
        irms = float(content[5][26:37])
        orms = lrati*(13.6e-6+3.5e-6*(ifreq-0.95)/(1.35-0.95))*(115.e3/1.e9/dfreq)**0.5
        norm = orms/irms
        snorm = '%.8f' % norm
        text = miriad_path+'maths exp="<c'+sichan+'.m>*'+snorm+'" out="n'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # expand the PB to full image size
        text = miriad_path+'regrid in="pb'+sichan+'.b" tin="c'+sichan+'.b" axes="1,2" out="pbr'+sichan+'.b"'
        subprocess.run(text, shell='True')
        # modify dirty beam for mosaic situation (multiply by PB)
        text = miriad_path+'maths exp="<c'+sichan+'.b>*<pbr'+sichan+'.b>" out="mos'+sichan+'.b"'
        subprocess.run(text, shell='True')
        norm = 3.0*orms
        sn = '%.8f' % norm
        bnorm = 5.0*orms
        sbn = '%.8f' % bnorm
        # pull out single channel from input cube
        text = miriad_path+'imsub in="'+indir+'sky_HI_'+cver+'" out="h'+sichan+'.m" region="images('+sinvchan+')"'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=h'+sichan+'.m/crval3 value='+sifreq
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=h'+sichan+'.m/crpix3 value=1.'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=h'+sichan+'.m/cdelt3 value='+sdfreq
        subprocess.run(text, shell='True')
        # extract suitable true and imperfect continuum images from coarsely sampled continuum cubes
        # first make mini-cube template so that the center plane can be interpolated
        text = miriad_path+'imframe in="h'+sichan+'.m" frame="-'+scrp12m+','+scrp12m+',-'+scrp12m+','+scrp12m+',0,2." out="hf'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in=hf'+sichan+'.m/crpix3 value=2.'
        subprocess.run(text, shell='True')
        text = miriad_path+'regrid in="'+outdir+'cube_frcont_'+cver+'" tin="hf'+sichan+'.m" axes="3" tol=0.0001 out="cic'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'imsub in="cic'+sichan+'.m" region="images(2)" out="ci'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'regrid in="'+outdir+'cube_frncont_'+cver+'" tin="hf'+sichan+'.m" axes="3" tol=0.0001 out="ncic'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'imsub in="ncic'+sichan+'.m" region="images(2)" out="nci'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # combine hi and continuum into imperfect sky image
        text = miriad_path+'maths exp="<h'+sichan+'.m>+<ci'+sichan+'.m>-<nci'+sichan+'.m>" out="sc'+sichan+'.m"'
        subprocess.run(text, shell='True')

        # generate hi opacity image from hi emission and subtract opacity scaled continuum if red-shift of continuum is greater than current channel
        jy2nh = Cfact*(1.+redz)**4
        sjy2nh = '%.6e' % jy2nh
        nh19 = 1.e19
        jy19 = nh19/jy2nh
        sjy19 = '%.6e' % jy19
        # use the model outlined in 2012_ApJ_749_87 to relate columnn density and HI opacity
        nhmin = 1.25e20
        snhmin = '%.6e' % nhmin
        nhmax = 7.5e21
        snhmax = '%.6e' % nhmax
        dvnom = 15.
        sdvnom = '%.6f' % dvnom
        # arbitrarily boost apparent column densities (above 1e19) by power law index nhb to enhance hi absorption statistics (making up for low resolution)
        # but then apply an asymptotic filtering with a tanh function to make sure we stay below nhmax
        nhb = 1.9
        snhb = '%.6f' % nhb
        # restrict occurence of HI absorption to continuum sources with brightness temp greater than Tcmin, since otherewise can not distinguish from HI emission
        # in reality such mixing will actually occur, but our Tc is artifically low due to beam size
        Tcmin = 100.
        Scmin = Tcmin*db**2*1.e-3/(1.36*wlcm**2)
        sScmin = '%.6e' % Scmin
        text = miriad_path+'histo in="h'+sichan+'.m" > histm'+sichan+'.txt'
        subprocess.run(text, shell='True')
        # some channels are empty, so must check for this
        sempty = 'All pixels are   0.00'
        with open('histm'+sichan+'.txt') as f:
            content = f.readlines()
                

        scheck = (content[4][0:21])
        if (scheck == sempty):
            hmax = 0.
        else:
            hmax = float(content[5][19:31])

        if (hmax > jy19):
            text = miriad_path+'maths exp="10**(19.+(log10(<h'+sichan+'.m>*'+sjy2nh+')-19.)*'+snhb+')" mask="<h'+sichan+'.m>.gt.'+sjy19+'" out=nhr'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="'+snhmax+'*(exp(2.*<nhr'+sichan+'.m>/'+snhmax+')-1.)/(exp(2.*<nhr'+sichan+'.m>/'+snhmax+')+1.)" out=nh'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="log(('+snhmax+'-'+snhmin+')/('+snhmax+'-(<nh'+sichan+'.m>)))" mask="<nh'+sichan+'.m>.gt.'+snhmin+'" out=taui'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'imblr in="taui'+sichan+'.m" out="tau'+sichan+'.m"'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="tau'+sichan+'.m" mask="<'+indir+'sky_cont_'+cver+'_z>.gt.'+sredz+'" out=tauz'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'imblr in="tauz'+sichan+'.m" out="taur'+sichan+'.m"'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="((1.-exp(-'+sdvnom+'*<taur'+sichan+'.m>/'+sdvel+'))*<ci'+sichan+'.m>)" mask="<ci'+sichan+'.m>.gt.'+sScmin+'" out=hiai'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'imblr in="hiai'+sichan+'.m" out="hiap'+sichan+'.m"'
            subprocess.run(text, shell='True')
            # update header to approximately reflect that effective beam size is reduced by the non-linear scaling of intensity
            text = miriad_path+'puthd in="hiap'+sichan+'.m/bmaj" value='+sdbar
            subprocess.run(text, shell='True')
            text = miriad_path+'puthd in="hiap'+sichan+'.m/bmin" value='+sdbar
            subprocess.run(text, shell='True')
            # keep record of all HI absorption of continuum 
            text = miriad_path+'convol map=hiap'+sichan+'.m fwhm='+sdb+' options="final" out='+taudir+'hia'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="<sc'+sichan+'.m>-<'+taudir+'hia'+sichan+'.m>" out=s'+sichan+'.m'
            subprocess.run(text, shell='True')
        else:
            text = miriad_path+'maths exp="<sc'+sichan+'.m>" out=s'+sichan+'.m'
            subprocess.run(text, shell='True')
            text = miriad_path+'maths exp="<h'+sichan+'.m>*1.e-10" out='+taudir+'hia'+sichan+'.m'
            subprocess.run(text, shell='True')
        
        # extract peaks (> 3 RMS) both positive and negative from smooth (db arcsec) sky image and reset blanks to zero
        text = miriad_path+'maths exp="<s'+sichan+'.m>-'+sn+'" mask="<s'+sichan+'.m>.gt.'+sn+'" out="ppm'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'maths exp="<s'+sichan+'.m>+'+sn+'" mask="<s'+sichan+'.m>.lt.-'+sn+'" out="pnm'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'imblr in="ppm'+sichan+'.m" out="pp'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'imblr in="pnm'+sichan+'.m" out="pn'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # generate residuals
        text = miriad_path+'maths exp="<s'+sichan+'.m>-<pp'+sichan+'.m>-<pn'+sichan+'.m>" out="r'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # deconvolve beam from sky residuals and update header for "Jy/pixel" units
        text = miriad_path+'convol map="r'+sichan+'.m" fwhm='+sdb+' options="divide" sigma=1.e-1 out="rd'+sichan+'.m"'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in="rd'+sichan+'.m/bmaj" value=0.0'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in="rd'+sichan+'.m/bmin" value=0.0'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in="rd'+sichan+'.m/bpa" value=0.0'
        subprocess.run(text, shell='True')
        text = miriad_path+'puthd in="rd'+sichan+'.m/bunit" value="JY/PIXEL"'
        subprocess.run(text, shell='True')
        # convolve residuals with mosaic dirty beam
        text = miriad_path+'convol map="rd'+sichan+'.m" beam="mos'+sichan+'.b" out="rdc'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # add back the smooth peaks
        text = miriad_path+'imcomb in="pp'+sichan+'.m,pn'+sichan+'.m,rdc'+sichan+'.m" options="nonorm" out="'+prdir+'pr'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # add in the noise image
        text = miriad_path+'imcomb in="'+prdir+'pr'+sichan+'.m,n'+sichan+'.m" options="nonormalise,relax,fqaver" out="'+prndir+'prn'+sichan+'.m"'
        subprocess.run(text, shell='True')
        # clean up intermediate products
        text = '/bin/rm -rf '+'*'+sichan+'.* '
        subprocess.run(text, shell='True')
        text = '/bin/rm -rf '+'pb'+sichan+'.b c'+sichan+'.uv c'+sichan+'.b c'+sichan+'.m  histr'+sichan+'.txt n'+sichan+'.m pbr'+sichan+'.b s'+sichan+'.m h'+sichan+'.m'
        subprocess.run(text, shell='True')
        text = '/bin/rm -rf '+'mos'+sichan+'.b ppm'+sichan+'.m pp'+sichan+'.m r'+sichan+'.m rd'+sichan+'.m rdc'+sichan+'.m ci'+sichan+'.m nci'+sichan+'.m'
        subprocess.run(text, shell='True')
        text = '/bin/rm -rf '+'sc'+sichan+'.m tau'+sichan+'.m hiai'+sichan+'.m tauz'+sichan+'.m taui'+sichan+'.m pnm'+sichan+'.m pn'+sichan+'.m nh'+sichan+'.m'
        subprocess.run(text, shell='True')
        text = '/bin/rm -rf '+'hf'+sichan+'.m hiap'+sichan+'.m cic'+sichan+'.m taur'+sichan+'.m ncic'+sichan+'.m nhr'+sichan+'.m histm'+sichan+'.txt'
        subprocess.run(text, shell='True')
        text = 'Finished channel: '+sichan
            
        return text

if __name__ == '__main__':
    run_fhia(config)



