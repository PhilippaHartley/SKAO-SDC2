import numpy as np
import os, glob, sys
from astropy.io import fits
import galsim
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import astropy
from astropy.cosmology import LambdaCDM
import configparser
import logging
import skymodel_tools
import time
import scipy
from skimage import transform 
from scipy.ndimage import zoom



def plot_properties(config):

   # datacube_dir = config.get('pipeline', 'datacube_dir')
   # data_path = config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'/'
   # cat_name = data_path+config.get('field', 'fits_prefix')+'_truthcat.fits'
   # diagnostics_dir = data_path+config.get('pipeline', 'data_path')+config.get('pipeline', 'diagnostics_dir')+'input/'



    base_dir = config.get('pipeline', 'base_dir')
    data_path = config.get('pipeline', 'data_path')
    datacube_dir = base_dir +data_path+config.get('pipeline', 'datacube_dir')
    prepared_dir = config.get('pipeline', 'prepared_dir')
    prepared_metadata = config.get('pipeline', 'prepared_metadata')

    cubes = np.loadtxt(datacube_dir+prepared_metadata, dtype = 'str')
    print(cubes)



    #if not os.path.isdir(diagnostics_dir):
     
    #    os.system('mkdir -p '+diagnostics_dir)

    names =  (cubes[:,8])
    survey = []

    for name in names:
        if "G" in name:
            survey = np.append(survey, "G")
        if "H" in name:
            survey = np.append(survey, "H")

    incl = (cubes[:,4]).astype(np.float)



    print (names)
    distance = cubes[:,2] # Mpc
    sys_vel = cubes[:,3].astype(np.float)*1000 # m/s, should be similar to rec vel but will include pec vel
  #  incl = (cubes[:,4].astype(np.float)/360.)*2*np.pi # from deg to radians
    rot_vel = cubes[:,5].astype(np.float) # this is corrected for inclination, according to hyperleda from which it was taken
    pa = cubes[:,6].astype(np.float) # make a new column for this : pas are in hyperleda thing
    log10_atlas_M_HI = cubes[:,7].astype(np.float)
    D_25 = (10**cubes[:,9].astype(np.float))/10
    log10_D_HI = cubes[:,10].astype(np.float)
    gal_class = cubes[:,1].astype(np.float)


    plt.scatter(log10_atlas_M_HI[survey=="H"], incl[survey=="H"], color= 'blue', label = 'HALOGAS sources')
    plt.scatter(log10_atlas_M_HI[survey=="G"], incl[survey=="G"], color= 'red', label = 'THINGS sources')

    plt.xlabel(r'log$_{10}$ M$_{\rm HI}$ (M$_{\odot})$')
    plt.ylabel('inclination (degrees)')
    plt.legend()
    plt.savefig('MHIvsinc.png')
# use mag abd sfr?
    plt.clf()

    M = np.linspace(7, 11, 2)
    D = (0.51*(M)) -3.32
    plt.scatter(log10_atlas_M_HI[gal_class<3], log10_D_HI[gal_class<3], color= 'red', label = 'spheroidal')
    plt.scatter(log10_atlas_M_HI[gal_class>7], log10_D_HI[gal_class>7], color= 'blue', label = 'late-type')
    plt.scatter(log10_atlas_M_HI[(gal_class<=7)&(gal_class>=3)], log10_D_HI[(gal_class<=7)&(gal_class>=3)], color= 'green', label = 'mid') 
    plt.plot(M,D,color = 'gray', label = 'Broeils and Rhee 1997')
    plt.xlabel(r'log$_{10}$ M$_{\rm HI}$ (M$_{\odot})$')
    plt.ylabel(r'log$_{10}$ D$_{\rm HI}$ (kpc)')
    plt.legend()
 
    plt.savefig('DHIvsMHI.png')


def prepare(config, diagnostics):    # rotates and also centres the cubes on dynamical centre

    t0 = time.time()


    doplot = 1#config.getfloat('pipeline', 'doplot')
    cmap_col = 'jet'
    base_dir = config.get('pipeline', 'base_dir')
    datacube_dir = config.get('pipeline', 'datacube_dir')
    prepared_dir = config.get('pipeline', 'prepared_dir')
    original_dir = config.get('pipeline', 'original_dir')
    coldens_dir = config.get('pipeline', 'coldens_dir')

    rms = config.getfloat('observation', 'rms_noise')
    channel_width_Hz = config.getfloat('observation', 'channel_width')
    rms_line_flux_Hz = rms*channel_width_Hz # Jy Hz
    rms_line_flux_km_s = (49.8/2.36e5)*rms_line_flux_Hz
    print ( rms_line_flux_Hz,rms_line_flux_km_s)
    


    if not os.path.exists(datacube_dir+prepared_dir):
        os.makedirs(datacube_dir+prepared_dir) # note that makedirs creates any intermediate dirs 
    if not os.path.exists(datacube_dir+'diagnostics'):
        os.makedirs(datacube_dir+'diagnostics')

    f = open(base_dir+datacube_dir+'metadata_appendednew.txt','a') 
    f.write('# ID\n') 
    f.write('# Class\n') 
    f.write('# Adopted distance\n')
    f.write('# Systemic velocity (km/s)\n')
    f.write('# Inclination (degrees)\n')
    f.write('# Rotation velocity (km/s)\n')
    f.write('# Position angle \n')
    f.write('# log10 HI mass from NEARGALCAT for HALOGAS sources, or from THINGS ascii for THINGS sources\n') 
    f.write('# Blanked cube name\n')
    f.write('# logD25, 0.1 arcmin, from hyperleda, isophotal level 25 mag/arcsec^2 in the B-band\n')
    f.write('# logDHI_maj, kpc, measured from atlas cube\n')
    f.write('# logDHI_min, kpc, measured from atlas cube\n')
    f.close() 

    cubes = np.loadtxt(base_dir+datacube_dir+'metadata_appended_LASs.txt', dtype = 'str')
    print (cubes)

    for count, i in enumerate(cubes):   # doesnt work for blanked due to unqual axes etc
        # collect together all preparation parts
        # acquire a sample from the atlas + things collection
        # get atlas source properties from metadata
        distance = float(i[2]) # Mpc
        sys_vel = float(i[3])*1000 # m/s, should be similar to rec vel but will include pec vel
        incldeg = (float(i[4])/360.) 
        incl = incldeg*2*np.pi # from deg to radians
        rot_vel = float(i[5]) # this is corrected for inclination, according to hyperleda from which it was taken
        pa = float(i[6]) # 
        atlas_M_HI = 10**float(i[7])
        D_25 = (10**float(i[9]))/10
        # get datacube and header properties
        cube_name = base_dir+datacube_dir+original_dir+'{}cr.fits'.format((i[8]))
        print ('cube name: ', cube_name)
        cube_fits = fits.open(cube_name)
        cube = np.copy(cube_fits[0].data)
        print ('cube shape: ',cube.shape)
        np.putmask(cube, cube<0, 0)
        if 'G' in i[8]:
            cube = cube[0]
        dx = 100 # pc   np.abs(cube_fits[0].header['CDELT1']*3600) # arcsec
        dy = 100 # pc   np.abs(cube_fits[0].header['CDELT2']*3600) # arcsec
        dz = cube_fits[0].header['CDELT3'] # m/s

        # test git
        print ('dx is: ', dx)
        velocities = cube_fits[0].header['CRVAL3']+(np.arange(cube.shape[0])*dz)
        atlas_bmaj =  cube_fits[0].header['BMAJ'] *3600
        atlas_bmin =  cube_fits[0].header['BMIN'] *3600
        dyn_centre_x = cube_fits[0].header['CRPIX1']
        dyn_centre_y = cube_fits[0].header['CRPIX2']
        im_centre_x = cube.shape[2]/2
        im_centre_y = cube.shape[1]/2

        pad_cube = np.copy(cube)

        while True:


            # shift to be centred on dynamical centre:
            #shift x axis:
            x_shift = int(im_centre_x-dyn_centre_x)
            shifted_cube = np.roll(pad_cube, x_shift, axis=2)

            #shift y axis:
            y_shift = int(im_centre_y-dyn_centre_y)
            shifted_cube = np.roll(shifted_cube, y_shift, axis=1)

            cube_fits[0].header['CRPIX1']+=x_shift# check this comes out to give centred crpix
            cube_fits[0].header['CRPIX2']+=y_shift
        
              
            # now do the rotation, need correct centre
            rot_cube = scipy.ndimage.rotate(shifted_cube, pa,  axes = (1,2))




            # do some preliminary resizing
            # highest frequency is 1.35GHz, which is equivalent to z = 0.052. Use z = 0.04 to start with
            # D_A at 0.04 is 164.2 Mpc
            # size = d_A * theta(rads) = 164.2 * 2/206265 = 0.001592 Mpc = 1592 pc pr pixel
            # could shrink by a factor of 10 to get 1000 pc per pixel, or 15 to get 1500 per pixel

            resize_v_scale = 1
            resize_scale = 100/1000 # use 1000 conservatively to allow for resizing due to mass and incl
            rot_cube = transform.resize(rot_cube, (np.ceil(rot_cube.shape[0]*resize_v_scale),\
                    np.ceil(rot_cube.shape[1]*resize_scale), np.ceil(rot_cube.shape[2]*resize_scale)),\
                    order = 3 ,preserve_range = True, anti_aliasing=True)
            


            # crop empty pixels
            # aim to crop pixels below line_flux_rms_km_s
            # distances were from Heald+2011
            flat_cube = np.sum(rot_cube, axis = 0) #Jy
            flat_cube_line_flux = flat_cube*(np.absolute(dz)/1000) # Jy km/s
            atlas_line_flux = atlas_M_HI/(2.36e5*distance**2)
            norm_factor = atlas_line_flux/np.sum(flat_cube_line_flux)
            flat_cube_line_flux*=norm_factor

            np.putmask(flat_cube_line_flux, flat_cube_line_flux<rms_line_flux_km_s, 0)

            sum_rot = np.copy(flat_cube_line_flux)
  
        



    
         #   # normalise the flattened cube to correct flux density, using M_HI
         #   atlas_line_flux = atlas_M_HI/(2.36e5*distance**2) # integration of F dz, Jy-km/s
         #   atlas_summed_flux = atlas_line_flux/(np.absolute(dz)/1000.) # convert dz from m/s to km/s: is what the summed flux should be using the atlas dv value
         #   print ('ska_summed_flux: ', atlas_summed_flux)
#            flat_cube_summed_flux = np.sum(flat_cube) # is what the summed flux of the atlas cube currently is
#            print ('flat_cube_summed_flux: ', flat_cube_summed_flux) 
#            norm_factor = atlas_summed_flux/flat_cube_summed_flux 
#            flat_cube*=norm_factor  





    

        ##    sum_rot = np.sum(rot_cube, axis = 0)
       
      ##      np.putmask(sum_rot,sum_rot<1e-10, 0)

            cube_size_array_maj =np.sum(sum_rot, axis = 1)
            #print (cube_size_array_maj)
            lhs = np.argmax(cube_size_array_maj>0 ) # finds first occurence of a True value
            rhs = np.argmax(cube_size_array_maj[::-1]>0) # same but backwards
         #   rhs = len(cube_size_array_maj) - rhs


            cube_size_array_min = np.sum(sum_rot, axis = 0)

    

            top = np.argmax(cube_size_array_min>0 ) # finds first occurence of a True value
            bottom = np.argmax(cube_size_array_min[::-1]>0) # same but backwards
           # bottom = len(cube_size_array_min) - bottom

            x_crop = np.min(([lhs, rhs]))
            y_crop = np.min(([top, bottom]))

            lhs = x_crop
            rhs = len(cube_size_array_maj)-x_crop

            top = y_crop
            bottom = len(cube_size_array_min)-y_crop
            print ('to crop, clockwise from top: ', top, rhs, bottom, lhs)

            rot_cube = rot_cube[:, lhs:rhs, top:bottom] # crop evenly to keep centred

    
            if doplot:


                sum_cube = np.sum(cube, axis = 0)

                plt.subplot(231)
                plt.imshow(sum_cube, cmap = cmap_col )
                plt.title('Original')
                plt.subplot(232)
                plt.imshow(np.sum(shifted_cube, axis = 0), cmap = cmap_col)  
                plt.title('Centred dynamically')     
                plt.subplot(233)
                plt.imshow(np.sum(rot_cube, axis = 0), cmap = cmap_col)
                plt.title('Rotated')


                plt.subplot(234)
                plt.imshow(sum_cube,norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min() ,
                                                  vmin=sum_cube.min(), vmax=sum_cube.max()), cmap = cmap_col )

                plt.subplot(235)
                plt.imshow(np.sum(shifted_cube, axis = 0),norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min() ,
                                                  vmin=sum_cube.min(), vmax=sum_cube.max()), cmap = cmap_col)  
           
                plt.subplot(236)
                plt.imshow(np.sum(rot_cube, axis = 0),norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min() ,
                                                  vmin=sum_cube.min(), vmax=sum_cube.max()), cmap = cmap_col)
   


                fig = plt.gcf()
                fig.savefig('diagnostics/input/'+i[8]+'_transform.pdf')
                plt.show()

          

                while True:

                    ok = input("Does it look OK (y/n)?: ")

                    if ok == 'y':



                        save_cube = 1
                        if save_cube:
                            yessave =  input('Writing new fits cube, are you sure (y/n)?')
                            if yessave == 'y':

                                prepared_cube_name = base_dir+datacube_dir+prepared_dir+'{}cr_rotated_shrunk.fits'.format((i[8]))
                                cube_fits[0].data = rot_cube
                                cube_fits.writeto(prepared_cube_name, overwrite = True)
                                print ('saved')
                            if yessave == 'n':
                                print ('not saved')
                        break

                    elif ok == 'n':
                        padding = np.int(input("Please enter amount of pixels for padding:"))                  
                        pad_cube = np.pad(cube, ((0,0),(padding,padding), (padding,padding)), 'constant', constant_values=(0)) # need to sort this padding out 
                        break

                    else:
                        continue
                if ok == 'n':
                    continue        
                else:
                    break





        atlas_x_maj = 0
        atlas_x_min = 0
        do_LAS_using_D_25 = 0
        do_LAS_using_M_HI_to_flux = 0 # move this to before cube is shrunk

        if do_LAS_using_M_HI_to_flux:# these values have been checked by eye and are in agreement with sizes of sources in blanked cubes
        # have deleted coldens cube and D_25 calculation methods

        # 13/03/20 this looks good now
        # need to account for inclination

      #  sin_i=ran_mwc(iseed)
      #     inclinations(i)=asin(sin_i)*180./pi
      #     ! all HI galaxies have spiral morphology. start from an intrinsic axis ratio of 0.2
      #     q=0.2+randgauss_boxmuller(iseed)*0.05
      #     q2=q**2 ! square alpha
      #     q=sqrt(q2+(cos(inclinations(i)))**2.*(1.-q2)) ! observed axis ratio, linked to inclination
      #     qrat(i)=q
      #     bmaj(i)=sqrt(sizes(i)**2./q) ! apparent bmaj
      #     bmin(i)=q*bmaj(i)     ! apparent bmin (edited) 



    #.  use cos**2(i) = (b/a)**2 - alpha**2 / 1 - (alpha)**2
  
    #   alpha = 0.2 # spiral galaxies; ellipticals = 0.5
    #       b_over_a = np.sqrt(((np.cos(incl)**2)*(1-alpha**2))+(alpha**2))

            
            # distances were from Heald+2011
            flat_cube = np.sum(rot_cube, axis = 0)
            print ('dist: ', distance, dz)
            # normalise the flattened cube to correct flux density, using M_HI
            atlas_line_flux = atlas_M_HI/(2.36e5*distance**2) # integration of F dz, Jy-km/s
            atlas_summed_flux = atlas_line_flux/(np.absolute(dz)/1000.) # convert dz from m/s to km/s: is what the summed flux should be using the atlas dv value
            print ('ska_summed_flux: ', atlas_summed_flux)
            flat_cube_summed_flux = np.sum(flat_cube) # is what the summed flux of the atlas cube currently is
            print ('flat_cube_summed_flux: ', flat_cube_summed_flux) 
            norm_factor = atlas_summed_flux/flat_cube_summed_flux 
            flat_cube*=norm_factor    # converts the current atlas flattened cube to what it should be, for the given dv of the unflattened cube



            # now measure D_HI using flux cut as a proxy for column density
            # convert flux density to solar masses. cut at 1 solar mass per parsec**2
            # each pixels in 100 pc along a side --> 10 000 pc**2
            # --> cut at 1*10 000 solar masses  = 1e4 per pixel

            flat_cube_M_HI = flat_cube* (2.36e5*distance**2)*(np.absolute(dz)/1000.) # converts each pixel in the flattened cube into a mass


            print ('check log10 summed mass: ', np.log10(np.sum(flat_cube_M_HI)))


            flat_cube_M_HI_masked = np.copy(flat_cube_M_HI)

            np.putmask(flat_cube_M_HI_masked, flat_cube_M_HI_masked<1e4, 0)# might need refinement due to wispy bits
            
            HI_size_array_maj =np.sum(flat_cube_M_HI_masked, axis = 1)
            HI_size_array_maj= HI_size_array_maj[HI_size_array_maj> 0]
       
            LAS_maj = len(HI_size_array_maj)*dx # pc
            print ('LAS_min: ', LAS_maj,  'pc')
   
            HI_size_array_min = np.sum(flat_cube_M_HI_masked, axis = 0)
            print (len(HI_size_array_min))
            HI_size_array_min= HI_size_array_min[HI_size_array_min> 0]    
            LAS_min = len(HI_size_array_min)*dx # pc
            print ('LAS_min: ', LAS_min,  'pc')
            logD_HI_maj_kpc = np.log10(len(HI_size_array_maj)*0.1) # convert from pixels to kpc
            logD_HI_min_kpc = np.log10(len(HI_size_array_min)*0.1)


            if doplot:




                plt.subplot(231)
                plt.imshow(flat_cube, cmap = cmap_col)
                plt.title('Normflux i %s'%incldeg)
                plt.subplot(232)
                plt.imshow(flat_cube_M_HI, cmap = cmap_col)

                plt.title('Mass/pixel')
                plt.subplot(233)
                plt.imshow(flat_cube_M_HI_masked,cmap = cmap_col)
         
                plt.title('Mask <1e4 M_sol')
                plt.subplot(234)
                plt.imshow(flat_cube, norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min()),cmap = cmap_col)
            
                plt.subplot(235)
                plt.imshow(flat_cube_M_HI,norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min() ,
                                              vmin=sum_cube.min(), vmax=sum_cube.max()), cmap = cmap_col)

 
                plt.subplot(236)
                plt.imshow(flat_cube_M_HI_masked,norm=colors.SymLogNorm(linthresh=sum_cube[sum_cube>0].min(), linscale=sum_cube[sum_cube>0].min() ,
                                              vmin=sum_cube.min(), vmax=sum_cube.max()), cmap = cmap_col)
         
                plt.text(5, -10, 'LAS = %s pix'%(LAS_maj/dx) )

                fig = plt.gcf()
                fig.savefig('../out/diagnostics/input/'+i[8]+'_LAS.pdf')
                plt.show()

            f = open(datacube_dir+'metadata_appended_LASs.txt','a') 
            f.write('%s %s %s %s %s %s %s %s %s %s %s %.2f %.2f\n'%(i[0], i[1], i[2], i[3], \
                                                i[4], i[5], i[6], i[7], i[8], i[9], logD_HI_maj_kpc, logD_HI_min_kpc  ))
            f.close()


        dofit = 0
        # fits a sigmoid to the velocity curce: probably not needed
        if dofit:# only do fit on original data - blanked will prob not fit as well since lower v res


            # now fit v_rot along major axis
            # parameterise velocity curve as sigmoid fit to major axis velocities; 
            # could use the method of fitting Gaussians along major axis, or:
            # TiRiFiC could be used 

            # make a cut along mojor axis to obtain fit constraints
            fit_args = np.argwhere(rot_cube[:,:,int((rot_cube.shape[2])/2)] > 200*np.median(rot_cube[:,:,int((rot_cube.shape[2])/2)])) # changed to get correct shape arg

            # get asymptotes from metadata
            sys_vel_channel = np.argmin(np.abs((sys_vel)-velocities))
            rot_vel_channel_from_sys_vel = rot_vel/(dz/1000.)
            print ('rot_vel_channel_from_sys_vel: ', rot_vel, dz,rot_vel_channel_from_sys_vel )
            print ('systemic velocity channel in cube: ', sys_vel_channel)
            xdata = fit_args[:,1]-(rot_cube.shape[1]/2) # using indexes, not values
            ydata = fit_args[:,0]-(rot_cube.shape[0]/2)
            sys_vel_channel-=rot_cube.shape[0]/2
            lower_asym = sys_vel_channel+rot_vel_channel_from_sys_vel 
            higher_asym = sys_vel_channel-rot_vel_channel_from_sys_vel # ! check right way round

            plt.plot(xdata*dx, ydata*(dz/1000), 'o', label='data')
            plt.show()

            popt, pcov = curve_fit(fsigmoid, xdata, ydata, method='dogbox',p0=[lower_asym, 0.1, higher_asym])#, bounds=([0., 600.],[0.01, 1200.])
            print (popt)# save popt for eacxh source: to metadata?
            x = np.linspace(-1, rot_cube.shape[1], 50)-(rot_cube.shape[1]/2)
            y = fsigmoid(x, *popt)

         
            if doplot:
                plt.clf()

                plt.plot(xdata*dx, ydata*(dz/1000), 'o', label='data')
                plt.plot(x*dx,y*(dz/1000.), label='fit')
                #plt.ylim(0, 1.05)
                plt.legend(loc='best')
                plt.xlabel('Radius from centre (arcsec)')
                plt.ylabel(r'$V_{\rm rot} {\rm (km/s)}$')
          #      plt.savefig(datacube_dir+'diagnostics/{}{:04d}-cube_rot_vel_fit.png'.format(i[0], int(i[1])))
                plt.show()

            ## save some plots from this for inspection and diagnostics
            f = open(datacube_dir+'metadata_appended.txt','a') 
            f.write('%s %s %s %s %s %s %s %s %s %s %s %.3f %.3f %.3f %.3f\n'%(i[0], i[1], i[2], i[3], \
                                                    i[4], i[5], i[6], i[7], i[8], i[9],  np.log10(atlas_M_HI), popt[0], popt[1], popt[2] )   )
            f.close()


        if diagnostics:
            plot_properties(config)    
                        

if __name__=='__main__':

    config = configparser.ConfigParser()
    config.read(sys.argv[1])



    diagnostics = 0

   # plot_properties(config)   

    prepare(config, diagnostics)

