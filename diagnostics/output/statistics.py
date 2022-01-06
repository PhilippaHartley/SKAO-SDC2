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
import matplotlib.colors as mcolors



#need to escape path name ( and space)





rad2arcsec = 206265 


def make_plots(config):


    clist = [(0, "red"),  (1, "blue")]
    rvb = mcolors.LinearSegmentedColormap.from_list("", clist)



    H = config.getfloat('cosmology', 'H')
    M = config.getfloat('cosmology', 'M')
    L = config.getfloat('cosmology', 'L')
    c = config.getfloat('cosmology', 'c')
    G = config.getfloat('cosmology', 'G')

    cosmo = LambdaCDM(H0 = H, Om0 = M, Ode0 = L)
   
    dnu = config.getfloat('observation', 'channel_width')
    diagnostics_dir = config.get('pipeline', 'diagnostics_dir')
    diagnostics_out = diagnostics_dir+'/output/'
    '''
    print (    diagnostics_dir)
    if not os.path.isdir(diagnostics_path):
        os.system('mkdir '+diagnostics_path)
    if not os.path.isdir(diagnostics_dir):
        os.system('mkdir '+diagnostics_dir)
    '''


   # data_path = config.get('pipeline', 'base_dir')+config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'/'
   # out_cat_name = data_path+config.get('field', 'fits_prefix')+'_truthcat.fits'
   # base_dir = config.get('pipeline', 'base_dir')
   # datacube_dir = base_dir +config.get('pipeline', 'datacube_dir')
   # prepared_dir = config.get('pipeline', 'prepared_dir')
   # prepared_metadata = config.get('pipeline', 'prepared_metadata')






    base_dir = config.get('pipeline', 'base_dir')
    data_path = config.get('pipeline', 'data_path')
    datacube_dir = base_dir +data_path+config.get('pipeline', 'datacube_dir')
    prepared_dir = config.get('pipeline', 'prepared_dir')
    prepared_metadata = config.get('pipeline', 'prepared_metadata')
    project_name = config.get('pipeline', 'project_name')
    out_cat_name = project_name+'/'+config.get('field', 'fits_prefix')+'_restricted_truthcat.fits'




    prepared_cubes = np.loadtxt(datacube_dir+prepared_metadata, dtype = 'str')
    atlas_M_HI = 10**prepared_cubes[:,8].astype(np.float)
    atlas_incl = prepared_cubes[:,4].astype(np.float)
    atlas_names = prepared_cubes[:,9]
    cubes = np.loadtxt(datacube_dir+prepared_metadata, dtype = 'str')


    out_cat_fits = fits.open(out_cat_name)
    out_cat = out_cat_fits[1].data
    out_cols = out_cat_fits[1].columns.info
  


    MHI = 10**out_cat['MHI']
    MHI_incl_scale = out_cat['MHI_incl_scale']
    MHI_incl_v_scale = out_cat['MHI_incl_v_scale']

    plt.hist(MHI_incl_scale, histtype=u'step', bins = np.linspace(0,5,50),label = 'spatial, minor axis')
    plt.hist(MHI_incl_v_scale, histtype=u'step', label = 'frequency', bins = np.linspace(0,5,50))

    plt.legend()
    plt.savefig(diagnostics_out+'binned_scalings.png')
    plt.clf()



    z = out_cat['z']
    SHI = np.array([out_cat['line_flux_integral']]).astype(np.float)
    Anna_DHI = np.array([out_cat['HI_size']]).astype(np.float)*2
    DHI = np.array([out_cat['new_HI_size_arcsec']]).astype(np.float)
    DHI_kpc = np.array([out_cat['new_HI_size']]).astype(np.float)

    incl = np.array([out_cat['i']]).astype(np.float)*np.pi/180
    plt.hist(incl[0], bins = 50)
    plt.xlabel('inclination (degrees)')
    plt.ylabel('N')
    plt.savefig(diagnostics_out+'binned_inclination.png')
    plt.clf()
   
    plt.hist(z, bins = 50)
    plt.xlabel('z')
    plt.ylabel('N')
    plt.savefig(diagnostics_out+'binned_z.png')
    plt.clf()

    atlas_sources = out_cat['Atlas_source']
    #plot binned occurences of atlas sources
    count_atlas_sources = Counter(out_cat['Atlas_source'])
    nbars = np.arange(len(count_atlas_sources))
    heights = list(count_atlas_sources.values())
    names = list(count_atlas_sources.keys())

    for j in range(len(atlas_names)):
        if atlas_names[j] not in names:
          
            names.append(atlas_names[j])
            heights.append(0)
    nbars = np.arange(len(heights))
            
    sorted_atlas_MHI = []
    sorted_atlas_incl = []
    for i in range(len(names)):     
        for j in range(len(atlas_names)):        
            if names[i] == atlas_names[j]:               
                sorted_atlas_MHI.append(atlas_M_HI[j])
                sorted_atlas_incl.append(atlas_incl[j])

            

    sour = 'H08'
   
    sorted_atlas_MHI = np.asarray(sorted_atlas_MHI)
    sorted_atlas_incl = np.asarray(sorted_atlas_incl)*np.pi/180
    names = np.asarray(names)

    plt.scatter((np.log10(sorted_atlas_MHI)), np.cos(sorted_atlas_incl),  color = 'blue' )

    plt.scatter((np.log10(MHI[atlas_sources!='p'])), np.cos(incl[0][atlas_sources!='p']), s = 0.8, color = 'red')
    plt.xlabel('log10 MHI')
    plt.ylabel('cos(incl)')

    plt.savefig(diagnostics_out+'MHIvsincl_counts.png')
    plt.clf()

    plt.figure(figsize=(12,8))
    plt.bar(nbars, heights, align = 'center', color= rvb(sorted_atlas_MHI/atlas_M_HI.max()) )
    plt.xticks(nbars, names, rotation=90)
    plt.xlabel('Atlas source')
    plt.ylabel('N counts')
    plt.savefig(diagnostics_out+'Sample_counts_MHI.png')
    plt.clf()


    plt.figure(figsize=(12,8))
    plt.bar(nbars, heights, align = 'center', color= rvb(sorted_atlas_incl/atlas_incl.max()) )
    plt.xticks(nbars, names, rotation=90)
    plt.xlabel('Atlas source')
    plt.ylabel('N counts')
    plt.savefig(diagnostics_out+'Sample_counts_incl.png')


if __name__=='__main__':




    config = configparser.ConfigParser()
    config.read(sys.argv[1])
    make_plots(config)



