# can compare with http://egg.astro.cornell.edu/alfalfa/pubs/a40_110825.pdf

import numpy as np
import os, glob, sys
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import astropy
from collections import Counter
from astropy.cosmology import LambdaCDM
import configparser


#need to escape path name ( and space)





rad2arcsec = 206265 


def make_plots(config):

    H = config.getfloat('cosmology', 'H')
    M = config.getfloat('cosmology', 'M')
    L = config.getfloat('cosmology', 'L')
    c = config.getfloat('cosmology', 'c')
    G = config.getfloat('cosmology', 'G')

    cosmo = LambdaCDM(H0 = H, Om0 = M, Ode0 = L)
    print (cosmo)
    dnu = config.getfloat('observation', 'channel_width')
    


    in_cat_name = config.get('field', 'catalogue')


    in_cat_fits = fits.open(in_cat_name)
    in_cat = in_cat_fits[1].data
    in_cols = in_cat_fits[1].columns.info

    in_cat = in_cat[in_cat['MHI']!=-100]
    in_cat = in_cat[in_cat['redshift']>0.00001]

    print (len(in_cat))






    dvel = 1 #km/s
    dvel_test = 25 # cones from dZ = ((config.getfloat('cosmology', 'c') * rest_freq)/(ref_freq**2))*dnu when ref_freq = rest_freq


    z = in_cat['redshift']



  #  diagnostics_path = data_path+'/'+config.get('pipeline', 'diagnostics_dir')
  #  diagnostics_dir = diagnostics_path+'output/'

    '''
    print (    diagnostics_dir)
    if not os.path.isdir(diagnostics_path):
 
        os.system('mkdir '+diagnostics_path)


    if not os.path.isdir(diagnostics_dir):
 
        os.system('mkdir '+diagnostics_dir)
    '''




    MHI = 10**in_cat['MHI']



    
    


    SHI = np.array([in_cat['HI flux']]).astype(np.float) # units state mJy but actually Jy?

    print (len(SHI[SHI>0.0002*24]))

    D_L = cosmo.luminosity_distance(z).value # Mpc
    D_A = cosmo.angular_diameter_distance(z).value # Mpc


    flux_integral = MHI / (49.8 * D_A**2)
    summed_flux = flux_integral/dnu



    flux_integral_kmpers = MHI/(2.36e5*D_A**2) # equal to int S dV, units jy-km/s
    # convert from Jy-km/s to Jy by dividing by channel width in km/s
    summed_flux_kmpers = flux_integral_kmpers/(dvel) # this enables normalisation to correct Jy value per velocity slice 

    print (SHI.min())


    flux_integral_kmpers_correct = MHI/(2.36e5*D_L**2) # equal to int S dV, units jy-km/s
    # convert from Jy-km/s to Jy by dividing by channel width in km/s
    summed_flux_kmpers_correct = flux_integral_kmpers_correct/(dvel_test) #

    print (summed_flux)
    print (SHI/summed_flux)
    '''

  #  plt.scatter( SHI, summed_flux)
    plt.title('SHI = MHI / (dnu * 49.8 * D_L**2), and dnu = 115e3 Hz')
    plt.scatter(  SHI,summed_flux)
    plt.xlabel('trecs SHI (Jy)')
    plt.ylabel('SHI from MHI (Duffy+12)')
    plt.savefig('SHI_from_MHI_D_L_and_delta_freq.png')
    plt.clf()
   # plt.scatter(summed_flux_kmper-SHI, summed_flux_kmpers)
   # plt.show()

    plt.title('SHI = MHI/(dvel * 2.36e5 * D_A**2), and dvel = 1 km/s')
    plt.scatter(  SHI,summed_flux_kmpers)
    plt.xlabel('trecs SHI (Jy)')
    plt.ylabel('SHI from MHI (Duffy+12)')
    plt.savefig('SHI_from_MHI_D_A_and_delta_vel.png')
    plt.clf()




    plt.title('both D_L, dnu vs dvel methods')
    plt.scatter(  summed_flux,summed_flux_kmpers_correct)
    plt.xlabel('trecs SHI (Jy)')
    plt.ylabel('SHI from MHI (Duffy+12)')
    plt.savefig('test.png')
    plt.clf()

    '''

    exit()
















    data_path = config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'/'
    out_cat_name = data_path+config.get('field', 'fits_prefix')+'_truthcat.fits'



    out_cat_fits = fits.open(out_cat_name)
    out_cat = out_cat_fits[1].data
    out_cols = out_cat_fits[1].columns.info


    MHI = 10**out_cat['MHI']

    z = out_cat['z']


    SHI = np.array([out_cat['line_flux_integral']]).astype(np.float)



    























    Anna_DHI = np.array([out_cat['HI_size']]).astype(np.float)*2
    DHI = np.array([out_cat['new_HI_size_arcsec']]).astype(np.float)
    DHI_kpc = np.array([out_cat['new_HI_size']]).astype(np.float)

    print (Anna_DHI, DHI)

    i = np.array([out_cat['i']]).astype(np.float)
    plt.hist(i[0])
    plt.show()
   





    D_A = cosmo.angular_diameter_distance(z).value # Mpc




    Anna_DHI_kpc = (Anna_DHI/rad2arcsec)*D_A*1000


    flux = MHI/(2.36e5*D_A**2) # equal to int S dV, units jy-km/s
    # convert from Jy-km/s to Jy by dividing by channel width in km/s
    summed_flux = flux/(dz/1000.) # this enables normalisation to correct Jy value per velocity slice 



    '''
    plt.figure(figsize = (8,6))


    #log10_atlas_D_HI = ((0.51*np.log10(atlas_M_HI)) -3.32)

    plt.scatter(np.log10(MHI), np.log10(Anna_DHI), alpha = 0.5, label = 'AB')
    plt.scatter(np.log10(MHI), np.log10(DHI), alpha = 0.5, label = 'PH')


 #   plt.scatter(z1,SHI1, alpha = 0.5)

    plt.xlabel(r'$\log10M_{\rm HI}$')
  #  plt.xlabel(r'S$_{HI}$')
    plt.ylabel(r'$\log10D_{\rm HI}$ (arcsec)')
    plt.show()



    plt.scatter(np.log10(MHI), np.log10(Anna_DHI_kpc), alpha = 0.5, label = 'AB')
    plt.scatter(np.log10(MHI), np.log10(DHI_kpc), alpha = 0.5, label = 'PH')


 #   plt.scatter(z1,SHI1, alpha = 0.5)

    plt.xlabel(r'$\log10M_{\rm HI}$')
  #  plt.xlabel(r'S$_{HI}$')
    plt.ylabel(r'$\log10D_{\rm HI}$ (kpc)')
    plt.legend()
    


    print (z)
    diff = Anna_DHI-DHI
    print (diff.shape)

#    plt.hist(diff[0], 100)



    plt.scatter(z, diff[0])
    plt.show()

   # plt.savefig('logDHI_vs_logMHIkpcbothcat2.png')

   # plt.scatter(MHI1, SHI1)
   # plt.scatter(MHI2, SHI2)
 


    #plot MHI-HIflux




    #plot MHI - HI size

    '''


    plt.hist(z, bins = 200)
    plt.xlabel('z')
    plt.ylabel('N')
    plt.show()
    




    #plot binned occurences of atlas sources
    count_atlas_sources = Counter(out_cat['Atlas_source'])
    nbars = np.arange(len(count_atlas_sources))
    heights = list(count_atlas_sources.values())
    names = list(count_atlas_sources.keys())
    plt.figure(figsize=(12,8))
    plt.bar(nbars, heights, align = 'center')
    plt.xticks(nbars, names, rotation=90)
    plt.xlabel('Atlas source')
    plt.ylabel('N counts')
    plt.savefig('Sample_counts.png')



if __name__=='__main__':


   # cosi = np.random.rand(1000000)
   # i = np.arcsin(cosi)

   # plt.hist(i, bins = 1000)
   # plt.show()
   # exit()
    config = configparser.ConfigParser()
    config.read(sys.argv[1])
    make_plots(config)



