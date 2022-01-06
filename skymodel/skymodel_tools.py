import numpy as np
import pdb
from astropy.io import fits
from astropy import wcs
import galsim
import time

def get_image_size(fov,  pixel_scale):   # ensure that image size is a power of 2
  image_size = np.ceil((fov*60)/(pixel_scale))
  power2 = 0
  if power2:
    exponent = np.floor(np.log2(image_size))
    # choose exponent to make new image smaller or equal to FoV specified in ini file (to avoid accidentally making giant files)
    image_size = np.int(2**exponent)
  fov = (image_size*pixel_scale)/60.
  return  fov, image_size


def get_spectral_size(bw, dnu): # linear in freq
  crpix3 =  np.ceil((bw/dnu)/2.) # round up
  n_chan = int(crpix3 *2)
  return crpix3, n_chan

def get_spectral_sampling(config, rest_freq, ref_freq, dnu):
  dZ = ((config.getfloat('cosmology', 'c') * rest_freq)/(ref_freq**2))*dnu
  return dZ # optical velocity interval

def add_padding(image_size, fov, base_freq, top_freq, pixel_scale, dnu):
  padding_arcsec = 0#360
  pad_cut_factor = 0.0 # this determines the fraction of the padding that is filled with sources (that will
                        # get cropped to the final FoV. The larger this fraction, the less empty space there is
                        # and the greater the chance of wrapping around
  pad_pixels = padding_arcsec/pixel_scale
  pad_image_size = int(np.ceil(image_size+pad_pixels)) # pad now to avoid subcube placement wraparound problems

  padding_arcmin = padding_arcsec/60
  pad_fov = fov+(padding_arcmin*pad_cut_factor) # pad fov slightly larger than final cube size, to include sources just outside 
  padding_freq = 0#0.1e9
  cut_freq = padding_freq*pad_cut_factor

  pad_base_freq = base_freq-padding_freq # same principle for freq, but done slightly differently
  pad_top_freq = top_freq+padding_freq
  cut_base_freq = base_freq-cut_freq
  cut_top_freq = top_freq+cut_freq



  return pad_image_size, pad_fov, pad_base_freq, pad_top_freq, cut_base_freq, cut_top_freq



def setup_wcs(config, ndim, cosmology=False,nu_axis=False):
  
  pixel_scale = config.getfloat('skymodel', 'pixel_scale')*galsim.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')

  fov, image_size = get_image_size((fov), pixel_scale/galsim.arcsec)

  dnu =  config.getfloat('observation', 'channel_width')


  rest_freq = config.getfloat('observation', 'rest_freq')
  base_freq = config.getfloat('observation', 'lowest_frequency')
  top_freq = config.getfloat('observation', 'highest_frequency')
  ref_freq = base_freq+ ((top_freq - base_freq) / 2)

  pad_image_size, pad_fov, pad_base_freq, pad_top_freq, cut_base_freq, cut_top_freq = add_padding(image_size, fov, base_freq, top_freq, pixel_scale/galsim.arcsec, dnu)

  
  #ra_field = config.get('field', 'field_ra') # now using deg in ini file
  #ra_field_gs = galsim.Angle.from_hms(ra_field)
 # dec_field = config.get('field', 'field_dec')
 # dec_field_gs = galsim.Angle.from_dms(dec_field)

  ra_field_gs =  config.getfloat('field', 'field_ra') 
  dec_field_gs =  config.getfloat('field', 'field_dec') 
  


  pad_bw = pad_top_freq - pad_base_freq



  dZ = get_spectral_sampling(config, rest_freq, ref_freq, dnu)

 
  pad_crpix3, pad_n_chan = get_spectral_size(pad_bw, dnu)

  ref_vel = config.getfloat('cosmology', 'c') * (rest_freq/ref_freq - 1)

  msname = config.get('pipeline', 'project_name')+'.ms'
  msname = config.get('pipeline', 'data_path')+msname
  imagename = msname+'.image'
  
  #channel_width = bw / (n_chan*n_ifs)
  #if_width = bw / n_ifs
  
  w = wcs.WCS(naxis=ndim)
  if ndim==4:
    if nu_axis:
      '''
      w.wcs.naxis = [float(image_size),
                     float(image_size),
                     1,
                     n_chan*n_ifs]
      '''
      w.wcs.crpix = [float(pad_image_size)/2,
                     float(pad_image_size)/2,
                     1,
                     1]
      w.wcs.cdelt = [-pixel_scale / galsim.degrees,
                     pixel_scale / galsim.degrees,
                     bw,
                     1]
      w.wcs.crval = [ra_field_gs / galsim.degrees,
                     dec_field_gs / galsim.degrees,
                     base_freq + bw/2,
                     1]
      w.wcs.ctype = ['RA---SIN',
                     'DEC--SIN',
                     'FREQ',
                     'STOKES']
      w.wcs.cunit = ['deg',
                     'deg',
                     'Hz',
                     '']
    
    else:
      '''
      w.wcs.naxis = [float(image_size),
                     float(image_size),
                     1,
                     1]
      '''
      w.wcs.crpix = [float(pad_image_size)/2,
                     float(pad_image_size)/2,
                     1,
                     1]
      w.wcs.cdelt = [-pixel_scale / galsim.degrees,
                     pixel_scale / galsim.degrees,
                     bw,
                     1]
      w.wcs.crval = [ra_field_gs ,
                     dec_field_gs ,
                     base_freq + bw/2,
                     1]
      w.wcs.ctype = ['RA---SIN',
                     'DEC--SIN',
                     'FREQ',
                     'STOKES']
      w.wcs.cunit = ['deg',
                     'deg',
                     'Hz',
                     '']

  elif ndim==3:
    

    # begin with a freq-based spectral axis
    w.wcs.crpix = [np.ceil((pad_image_size)/2),
                   np.ceil((pad_image_size)/2), pad_crpix3]
    w.wcs.cdelt = [-pixel_scale / galsim.degrees,
                   pixel_scale / galsim.degrees,
                   dZ]
    w.wcs.crval = [ra_field_gs ,
                   dec_field_gs , ref_vel]
    w.wcs.ctype = ['RA---SIN',
                   'DEC--SIN', 
                   'VOPT-F2W'] # see https://www.astro.rug.nl/software/kapteyn/spectralbackground.html (FELO_HEL behaves as VOPT_F2W)
    w.wcs.cunit = ['deg',
                   'deg',
                   'm/s']
    w.wcs.restfrq = rest_freq   




    '''

    w.wcs.crpix = [float(image_size)/2,
                   float(image_size)/2, crpix3]
    w.wcs.cdelt = [-pixel_scale / galsim.degrees,
                   pixel_scale / galsim.degrees,
                   channel_width_vel]
    w.wcs.crval = [ra_field_gs / galsim.degrees,
                   dec_field_gs / galsim.degrees, ref_vel]
    w.wcs.ctype = ['RA---SIN',
                   'DEC--SIN', 
                   'FELO-HEL'] # see https://www.astro.rug.nl/software/kapteyn/spectralbackground.html (FELO_HEL behaves as VOPT_F2W)
    w.wcs.cunit = ['deg',
                   'deg',
                   'm/s']
    w.wcs.restfrq = rest_freq       
    '''       

      
  elif ndim==2:
    w.wcs.crpix = [np.ceil((pad_image_size)/2),
                   np.ceil((pad_image_size)/2)]
    w.wcs.cdelt = [-pixel_scale / galsim.degrees,
                   pixel_scale / galsim.degrees]
    w.wcs.crval = [ra_field_gs ,
                   dec_field_gs ]
    w.wcs.ctype = ['RA---SIN',
                   'DEC--SIN']
    w.wcs.cunit = ['deg',
                   'deg']
    
  return w

def write4dImage(outname, image_data, pixel_scale,
                 obs_ra=150.e0, obs_dec=90.e0,
                 obs_freq=1.4e9, obs_bwidth = 1.e9,
                 clobber=True):
  '''Write a FITS file with minimum necessary information to use as input for a
  CASA observation simulation.

  Parameters
  ----------
  outname : string
    Output file name and location.
  image_data : ndarray
    Numpy array containing image data. Currently expects a 2d array.
  pixel_scale : float
    Scale of pixels in image ***IN DEGREES***
  obs_ra : float, optional, default: 150.e0
    Right Ascension of field **IN DEGREES***
  obs_dec : float, optional, default: 90.e0
    Declination of field ***IN DEGREES***
  obs_freq : float, optional, default: 1.4e9
    Frequency of emission in image ***in Hz***
  obs_bwidth : float, optional, default: 1.e9
    Bandwidth of image ***in Hz***. Should cover entire expected
    bandwidth of simulated observation.

  '''
  image_dim = image_data.ndim
  #TODO: some error checking on ndim (i.e. if already 4)

  new_image = np.zeros([1, 1, image_data.shape[0], image_data.shape[1]])
  new_image[0,0,:,:] = image_data
  new_image = np.asarray(new_image, dtype=np.float32)

  newhdr = fits.Header()

  newhdr['SIMPLE'] = True
  newhdr['BITPIX'] = -32
  newhdr['NAXIS'] = 4
  newhdr['NAXIS1'] = image_data.shape[0]
  newhdr['NAXIS2'] = image_data.shape[1]
  newhdr['NAXIS3'] = 1
  newhdr['NAXIS4'] = 1
  
  newhdr['EXTEND'] = True
  newhdr['BSCALE'] = 1.e0
  newhdr['BZERO'] = 0.e0
  newhdr['BMAJ'] = pixel_scale
  newhdr['BMIN'] = pixel_scale
  newhdr['BTYPE'] = 'Intensity'
  newhdr['OBJECT'] = '        '

  newhdr['BUNIT'] = 'JY/PIXEL'
  newhdr['EQUINOX'] = 2.e3
  newhdr['RADESYS'] = 'FK5'
  newhdr['LONPOLE'] = 1.8e2
  newhdr['LATPOLE'] = 6.8e1
  newhdr['PC01_01'] =   1.000000000000E+00
  newhdr['PC02_01'] =   0.000000000000E+00
  newhdr['PC03_01'] =   0.000000000000E+00
  newhdr['PC04_01'] =   0.000000000000E+00
  newhdr['PC01_02'] =   0.000000000000E+00
  newhdr['PC02_02'] =   1.000000000000E+00
  newhdr['PC03_02'] =   0.000000000000E+00
  newhdr['PC04_02'] =   0.000000000000E+00
  newhdr['PC01_03'] =   0.000000000000E+00
  newhdr['PC02_03'] =   0.000000000000E+00
  newhdr['PC03_03'] =   1.000000000000E+00
  newhdr['PC04_03'] =   0.000000000000E+00
  newhdr['PC01_04'] =   0.000000000000E+00
  newhdr['PC02_04'] =   0.000000000000E+00
  newhdr['PC03_04'] =   0.000000000000E+00
  newhdr['PC04_04'] =   1.000000000000E+00
  newhdr['CTYPE1'] = 'RA---SIN'
  newhdr['CRVAL1'] = obs_ra
  newhdr['CDELT1'] = pixel_scale
  newhdr['CRPIX1'] = float(new_image.shape[2]/2)
  newhdr['CUNIT1'] = 'deg'
  
  newhdr['CTYPE2'] = 'DEC--SIN'
  newhdr['CRVAL2'] = obs_dec
  newhdr['CDELT2'] = pixel_scale
  newhdr['CRPIX2'] = float(new_image.shape[3]/2)
  newhdr['CUNIT2'] = 'deg'
  
  newhdr['CTYPE3'] = 'STOKES'
  newhdr['CRVAL3'] = 1.e0
  newhdr['CDELT3'] = 1.e0
  newhdr['CRPIX3'] = 1.e0
  newhdr['CUNIT3'] = '        '
  
  newhdr['CTYPE4'] = 'FREQ'
  newhdr['CRVAL4'] = obs_freq
  newhdr['CDELT4'] = obs_bwidth
  newhdr['CRPIX4'] = 1.e0
  newhdr['CUNIT4'] = 'Hz      '

  newhdr['PV2_1'] = 0.e0
  newhdr['PV2_2'] = 0.e0

  newhdr['RESTFRQ'] =   obs_freq
  newhdr['SPECSYS'] = 'LSRK    '
  newhdr['ALTRVAL'] =  -0.e0
  newhdr['ALTRPIX'] =   1.000000000000E+00
  newhdr['VELREF'] =                  257
  newhdr['TELESCOP'] = 'SKYMODEL '
  newhdr['OBSERVER'] = 'A. R. Bitrary'
  newhdr['DATE-OBS'] = '2000-01-01T00:00:00.000100'
  newhdr['TIMESYS'] = 'UTC     '
  newhdr['OBSRA'] = obs_ra
  newhdr['OBSDEC'] = obs_dec
  newhdr['OBSGEO-X'] =   2.225061873184E+06
  newhdr['OBSGEO-Y'] =  -5.440061952280E+06
  newhdr['OBSGEO-Z'] =  -2.481682085791E+06
  newhdr['DATE'] = '2014-02-20T13:17:20.141840'
  newhdr['ORIGIN'] = 'write4dImage'

  #pdb.set_trace()

  hdu = fits.PrimaryHDU(new_image, header=newhdr)
  hdu.header.insert(8, ('EXTEND', 'BSCALE', 'BZERO'))
  hdu.header['EXTEND'] = True
  hdu.header['BSCALE'] = 1.e0
  hdu.header['BZERO'] = 0.e0
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(outname, clobber=clobber)

def ellip_pdf(x, B=0.19, C=0.58):
  '''Ellipticity probability distribution function from GREAT08
  (arXiv:0908.0945 equation A3).

  Parameters
  ----------
  x : float
    Ellipticity modulus at which to calculate PDF
  B : float
    Scale parameter
  C : float
    Index parameter
  '''
  retVar = x*(np.cos(np.pi*x/2.)**2.)*np.exp(-(2.*x/B)**C)
  return retVar
