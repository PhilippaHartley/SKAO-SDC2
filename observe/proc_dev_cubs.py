import subprocess
import numpy as np
import math
import os


def run_cubs(config):
    cver = config.get("pipeline", "cver")
    outdir = config.get("pipeline", "outdir")
    prdir = cver+'_fpr/'
    prndir = cver+'_fprn/'
    taudir = cver+'_tau/'

    if os.path.exists(outdir+'sky_'+cver+'.fits'):
        print ('**** message from pipeline: '+outdir+'sky_'+cver+'.fits already exists')
        print ('**** message from pipeline: not running proc_cubs')
        return()

    os.chdir(outdir+prdir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr00*.m" out=sub00'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr01*.m" out=sub01'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr02*.m" out=sub02'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr03*.m" out=sub03'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr04*.m" out=sub04'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr05*.m" out=sub05'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr06*.m" out=sub06'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr07*.m" out=sub07'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr08*.m" out=sub08'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr09*.m" out=sub09'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr10*.m" out=sub10'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr11*.m" out=sub11'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr12*.m" out=sub12'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr13*.m" out=sub13'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr14*.m" out=sub14'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr15*.m" out=sub15'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr16*.m" out=sub16'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr17*.m" out=sub17'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr18*.m" out=sub18'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr19*.m" out=sub19'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr20*.m" out=sub20'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr21*.m" out=sub21'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr22*.m" out=sub22'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr23*.m" out=sub23'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr24*.m" out=sub24'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr25*.m" out=sub25'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr26*.m" out=sub26'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr27*.m" out=sub27'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr28*.m" out=sub28'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr29*.m" out=sub29'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr30*.m" out=sub30'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr31*.m" out=sub31'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr32*.m" out=sub32'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr33*.m" out=sub33'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr34*.m" out=sub34'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr35*.m" out=sub35'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr36*.m" out=sub36'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr37*.m" out=sub37'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr38*.m" out=sub38'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr39*.m" out=sub39'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr40*.m" out=sub40'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr41*.m" out=sub41'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr42*.m" out=sub42'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr43*.m" out=sub43'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr44*.m" out=sub44'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr45*.m" out=sub45'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr46*.m" out=sub46'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr47*.m" out=sub47'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr48*.m" out=sub48'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr49*.m" out=sub49'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr50*.m" out=sub50'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr51*.m" out=sub51'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr52*.m" out=sub52'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr53*.m" out=sub53'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr54*.m" out=sub54'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr55*.m" out=sub55'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr56*.m" out=sub56'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr57*.m" out=sub57'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr58*.m" out=sub58'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr59*.m" out=sub59'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr60*.m" out=sub60'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr61*.m" out=sub61'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr62*.m" out=sub62'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr63*.m" out=sub63'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr64*.m" out=sub64'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr65*.m" out=sub65'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="pr66*.m" out=sub66'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="sub*" out=fpr_cube'
    subprocess.run(text, shell='True')
    text = '/bin/rm -rf sub* '
    subprocess.run(text, shell='True')

    os.chdir('../'+prndir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn00*.m" out=sub00'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn01*.m" out=sub01'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn02*.m" out=sub02'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn03*.m" out=sub03'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn04*.m" out=sub04'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn05*.m" out=sub05'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn06*.m" out=sub06'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn07*.m" out=sub07'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn08*.m" out=sub08'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn09*.m" out=sub09'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn10*.m" out=sub10'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn11*.m" out=sub11'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn12*.m" out=sub12'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn13*.m" out=sub13'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn14*.m" out=sub14'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn15*.m" out=sub15'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn16*.m" out=sub16'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn17*.m" out=sub17'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn18*.m" out=sub18'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn19*.m" out=sub19'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn20*.m" out=sub20'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn21*.m" out=sub21'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn22*.m" out=sub22'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn23*.m" out=sub23'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn24*.m" out=sub24'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn25*.m" out=sub25'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn26*.m" out=sub26'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn27*.m" out=sub27'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn28*.m" out=sub28'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn29*.m" out=sub29'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn30*.m" out=sub30'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn31*.m" out=sub31'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn32*.m" out=sub32'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn33*.m" out=sub33'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn34*.m" out=sub34'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn35*.m" out=sub35'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn36*.m" out=sub36'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn37*.m" out=sub37'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn38*.m" out=sub38'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn39*.m" out=sub39'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn40*.m" out=sub40'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn41*.m" out=sub41'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn42*.m" out=sub42'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn43*.m" out=sub43'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn44*.m" out=sub44'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn45*.m" out=sub45'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn46*.m" out=sub46'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn47*.m" out=sub47'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn48*.m" out=sub48'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn49*.m" out=sub49'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn50*.m" out=sub50'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn51*.m" out=sub51'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn52*.m" out=sub52'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn53*.m" out=sub53'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn54*.m" out=sub54'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn55*.m" out=sub55'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn56*.m" out=sub56'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn57*.m" out=sub57'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn58*.m" out=sub58'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn59*.m" out=sub59'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn60*.m" out=sub60'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn61*.m" out=sub61'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn62*.m" out=sub62'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn63*.m" out=sub63'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn64*.m" out=sub64'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn65*.m" out=sub65'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="prn66*.m" out=sub66'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="sub*" out=fprn_cube'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/fits in=fprn_cube op=xyout out=../sky_'+cver+'.fits'
    subprocess.run(text, shell='True')
    text = '/bin/rm -rf sub* '
    subprocess.run(text, shell='True')
    os.system('pwd')
    os.chdir('../'+taudir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia00*.m" out=sub00'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia01*.m" out=sub01'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia02*.m" out=sub02'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia03*.m" out=sub03'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia04*.m" out=sub04'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia05*.m" out=sub05'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia06*.m" out=sub06'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia07*.m" out=sub07'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia08*.m" out=sub08'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia09*.m" out=sub09'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia10*.m" out=sub10'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia11*.m" out=sub11'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia12*.m" out=sub12'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia13*.m" out=sub13'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia14*.m" out=sub14'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia15*.m" out=sub15'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia16*.m" out=sub16'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia17*.m" out=sub17'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia18*.m" out=sub18'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia19*.m" out=sub19'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia20*.m" out=sub20'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia21*.m" out=sub21'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia22*.m" out=sub22'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia23*.m" out=sub23'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia24*.m" out=sub24'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia25*.m" out=sub25'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia26*.m" out=sub26'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia27*.m" out=sub27'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia28*.m" out=sub28'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia29*.m" out=sub29'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia30*.m" out=sub30'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia31*.m" out=sub31'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia32*.m" out=sub32'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia33*.m" out=sub33'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia34*.m" out=sub34'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia35*.m" out=sub35'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia36*.m" out=sub36'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia37*.m" out=sub37'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia38*.m" out=sub38'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia39*.m" out=sub39'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia40*.m" out=sub40'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia41*.m" out=sub41'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia42*.m" out=sub42'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia43*.m" out=sub43'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia44*.m" out=sub44'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia45*.m" out=sub45'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia46*.m" out=sub46'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia47*.m" out=sub47'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia48*.m" out=sub48'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia49*.m" out=sub49'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia50*.m" out=sub50'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia51*.m" out=sub51'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia52*.m" out=sub52'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia53*.m" out=sub53'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia54*.m" out=sub54'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia55*.m" out=sub55'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia56*.m" out=sub56'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia57*.m" out=sub57'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia58*.m" out=sub58'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia59*.m" out=sub59'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia60*.m" out=sub60'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia61*.m" out=sub61'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia62*.m" out=sub62'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia63*.m" out=sub63'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia64*.m" out=sub64'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia65*.m" out=sub65'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="hia66*.m" out=sub66'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/imcat in="sub*" out=hia_cube'
    subprocess.run(text, shell='True')
    text = '/bin/rm -rf sub* '
    subprocess.run(text, shell='True')

    os.chdir('../')


    os.chdir(taudir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="hia_cube" mom=0 out=hia_cube_m0'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="hia_cube" mom=1 out=hia_cube_m1'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="hia_cube" mom=2 out=hia_cube_m2'
    subprocess.run(text, shell='True')

    os.chdir('../')

    os.chdir(prndir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fprn_cube" mom=0 out=fprn_cube_m0'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fprn_cube" mom=1 out=fprn_cube_m1'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fprn_cube" mom=2 out=fprn_cube_m2'
    subprocess.run(text, shell='True')

    os.chdir('../')

    os.chdir(prdir)

    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fpr_cube" mom=0 out=fpr_cube_m0'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fpr_cube" mom=1 out=fpr_cube_m1'
    subprocess.run(text, shell='True')
    text = '/usr/local/miriad_ATNF_1024_ants/linux64/bin/moment in="fpr_cube" mom=2 out=fpr_cube_m2'
    subprocess.run(text, shell='True')

    os.chdir('../../..')


if __name__ == '__main__':
    run_cubs()

