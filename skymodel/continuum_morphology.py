import configparser
import logging
import sys
import time
import numpy as np
import scipy
import scipy.ndimage
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D
from astropy.modeling.models import Sersic2D
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from skimage import transform

def make_img(
    config,
    ska_flux,
    ska_alpha,
    ska_freqmin,
    freqs,
    ska_size,
    ska_min,
    ska_PA,
    ska_dx,
    psf_maj,
    radioclass,
    corefrac,
    ska_ranid,
):

    doplot = config.getboolean("pipeline", "doplot")

    datacube_dir = config.get("pipeline", "datacube_dir")
    prepared_dir = config.get("pipeline", "prepared_dir")
    prepared_metadata = config.get("pipeline", "prepared_metadata")

    nfreqs = len(freqs)

    # for the mild convolution with the beam
    # this is the simple psf beam
    FWHM = psf_maj / ska_dx
    sigma2FWHM = 2 * np.sqrt(2 * np.log(2))
    Gaus_sigma = FWHM / sigma2FWHM  # 'g'
    kernel = Gaussian2DKernel(x_stddev=Gaus_sigma)
    size_in_pixels = ska_size / ska_dx
    Gaussize_in_pixels = size_in_pixels  # FWHM - this is true for flat-spectrum AGN
    Gauss_minor_size_in_pixels = (
        Gaussize_in_pixels  # initialised for a circularly symmetric source
    )
    smoothing_done = 0  # initialise flag to control when mild convolution is done

    if radioclass < 4:
        Gaussize_in_pixels = (
            Gaussize_in_pixels * 1.4241
        )  # scale length to FHWM - true for SFGs
        Gauss_minor_size_in_pixels = (
            Gaussize_in_pixels * ska_min / ska_size
        )  # minor axis

    if radioclass == 6:
        Gaussize_in_pixels = (
            size_in_pixels / 3.0
        )  # LAS approximated with a 3sigma - it varies for different images and true distribution is not Gaussian
        Gauss_minor_size_in_pixels = (
            Gaussize_in_pixels  # initialised for a circularly symmetric source
        )

    ska_ellip = 1.0 - ska_min / ska_size
    beamsigma_in_pixels = psf_maj / sigma2FWHM / ska_dx

    # initialization for unresolved source
    sigma, mu = (
        beamsigma_in_pixels,
        0.0,
    )  # exponential scale lenght converted to sigma in pixels

    npixs = int(sigma * 8.0) + 1
    x, y = np.meshgrid(
        np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
        np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
    )

    d = np.sqrt(x * x + y * y)
    g = np.exp(-((d - mu) ** 2 / (2.0 * sigma ** 2)))
    g = g / np.sum(g)  # normalization
    cube2 = g
    is_unresolved = 1
    cube_name = "SFG unresolved"

    if (Gauss_minor_size_in_pixels > 3.0) and (radioclass < 4):
        # resolved SFG. exponential profile
        sigma, mu = Gaussize_in_pixels, 0.0

        npixs = int(sigma * 8.0) + 1
        x, y = np.meshgrid(
            np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
            np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
        )
        mod = Sersic2D(amplitude=1, r_eff=size_in_pixels, n=1, ellip=ska_ellip)
        # the PA is default, corresponding to positive x axis. the rotation is done later at the postage stamp level, CCW from positive x axis

        g = mod(x, y)
        g = g / np.sum(g)  # normalization
        cube2 = np.array(g)

        cube_name = "SFG resolved"
        is_unresolved = 0

    if (Gaussize_in_pixels > 3.0) and (radioclass >= 4) and (radioclass < 6):
        # flat-spectrum AGN source. A Gaussian lobe with a Gaussian core

        # the source is a Gaussian with FWHM=source size, plus a gaussian core of 2 pixels size
        # gauss2=randomn(seed,nrows)*0.1+0.75
        sigma, mu = Gaussize_in_pixels / sigma2FWHM, 0.0

        npixs = int(sigma * 8.0) + 1

        x, y = np.meshgrid(
            np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
            np.linspace(-npixs / 2.0, npixs / 2.0, npixs),
        )
        d = np.sqrt(x * x + y * y)
        g = np.exp(-((d - mu) ** 2 / (2.0 * sigma ** 2)))
        g = g / np.sum(g)  # normalization
        sigma, mu = beamsigma_in_pixels, 0.0  # beam size in pixels
        gcore = np.exp(-((d - mu) ** 2 / (2.0 * sigma ** 2)))
        gcore = gcore / np.sum(gcore)  # normalization
        g = g * (1.0 - corefrac) + gcore * corefrac
        g = g / np.sum(g)  # normalization
        cube2 = g
        cube_name = "Flat-AGN resolved"
        is_unresolved = 0

    if (radioclass == 6) and (Gaussize_in_pixels > 3.0):
        # SS resolvd AGN: use a postage stamp from a real image
     
        is_unresolved = 0
        prepared_cubes = np.loadtxt(datacube_dir + prepared_dir + prepared_metadata, dtype="str")
      
        # choose atlas cube
        atlas_ranid = prepared_cubes[:, 1].astype(np.float)  
        distance_to_sample = abs(atlas_ranid - ska_ranid)
        distance_min_arg = np.argmin(distance_to_sample)
     
        logging.info("Chosen atlas source: %s", prepared_cubes[distance_min_arg, 0])

        # get datacube and header properties
        cube_name = prepared_cubes[distance_min_arg, 0]  
          
        # PH: temporary path fix
        cube_name = datacube_dir + prepared_dir + cube_name.split(
            "/home/a.bonaldi/data-cold-for-backup/data_challenges/inputs/AGN_library/processed/"
        )[1]
       
        cube_fits = fits.open(cube_name)
        cube = cube_fits[0].data

        logging.info("Atlas sample shape: %s", cube.shape)
       
        if doplot:
            plt.scatter(np.arange(len(cube_summed_spectrum)), cube_summed_spectrum)
            plt.show()

        dx = np.abs(cube_fits[0].header["CDELT1"] * 3600)  # arcsec  pix reso
        atlas_bmaj = cube_fits[0].header["SIZE"]
        atlas_psf = np.abs(cube_fits[0].header["PSF"])  # arcsec
        atlas_bmaj_px = atlas_bmaj / dx
        zoom_factor = ska_size / ska_dx / atlas_bmaj_px

        logging.info("atlas agn size in arcsec %f", atlas_bmaj)
        logging.info("target agn size in arcsec %f", ska_size)
        logging.info("atlas agn size in pixels %f", atlas_bmaj_px)
        logging.info("target agn size in pixels %f", ska_size / ska_dx)
        logging.info("zoom factor %f", zoom_factor)

        new_a_size = np.ceil((cube.shape[0] * zoom_factor))
        new_b_size = np.ceil((cube.shape[1] * zoom_factor))

        logging.info("postage stamp size %f", new_a_size)
        is_unresolved = 0

        if new_a_size < 3.0:
            is_unresolved = 1

        halfs = new_a_size / 2.0
        if np.ceil(halfs) == np.floor(halfs):
            # a bit of padding on both sides to leave the image centered as before but to get an odd number of pixels in the end

            npad = int(
                1.0 / zoom_factor / 2.0
            )  # pixels in original size to correspond to 0.5 pixel in the new image

            cube = np.pad(
                cube, ((npad, npad), (npad, npad)), "constant", constant_values=(0)
            )

            new_a_size = np.ceil((cube.shape[0] * zoom_factor))
            new_b_size = np.ceil((cube.shape[1] * zoom_factor))

        cube2 = transform.resize(
            cube, (new_a_size, new_b_size), order=1, preserve_range=True
        )
        cube2 = cube2 / np.sum(cube2)

        # we want the final postage stamp to have an odd number of pixels and be centered on the central pixel. this is for consistency between continuum and HI

        # cube is now resized. Add smoothing that takes into account original resolution of atlas AGN
        logging.info("original atlas PSF %f", atlas_psf)
        logging.info("after scaling %f", atlas_psf * zoom_factor)
        logging.info("atlas psf in pixels %f", atlas_psf * zoom_factor / ska_dx)

        psf_smooth = np.sqrt(
            np.square(psf_maj / ska_dx) - np.square(atlas_psf * zoom_factor / ska_dx)
        )  # additional smoothing PSF in units of the map pixel

      
        if psf_smooth > 0.0:

            atlas_smoo = psf_smooth / sigma2FWHM
          
            atlas_kernel = Gaussian2DKernel(x_stddev=atlas_smoo)
            size_a = cube2.shape[0]
            size_b = cube2.shape[1]
            padding_a = int(np.ceil(size_a / 2.0))
            padding_b = int(np.ceil(size_b / 2.0))

            cube2_padded = np.pad(
                cube2,
                ((padding_a, padding_a), (padding_b, padding_b)),
                "constant",
                constant_values=(0),
            )

            astropy_conv = convolve(
                cube2_padded, atlas_kernel, boundary="extend", normalize_kernel=True
            )

            cube2 = astropy_conv[
                padding_a : padding_a + size_a, padding_b : padding_b + size_b
            ]
            cube2 = cube2 / np.sum(cube2)

        smoothing_done = 1  # flag so that it is not smoothed later (whether I did or not because in that case atlas was low reso)
     
    if smoothing_done == 0:
        # convolution done on all sources - resolved and unresolved - check
    
        size_a = cube2.shape[0]
        size_b = cube2.shape[1]
        padding_a = int(np.ceil(size_a / 2.0))
        padding_b = int(np.ceil(size_b / 2.0))
        # padding_a=cube2.shape[0]
        # padding_b=cube2.shape[1]

        cube2_padded = np.pad(
            cube2,
            ((padding_a, padding_a), (padding_b, padding_b)),
            "constant",
            constant_values=(0),
        )

        # source mild convolution
        astropy_conv = convolve(
            cube2_padded, kernel, boundary="extend", normalize_kernel=True
        )

        cube2 = astropy_conv[
            padding_a : padding_a + size_a, padding_b : padding_b + size_b
        ]
        cube2 = cube2 / np.sum(cube2)

    nfreqs = freqs.size
    new_a_size = cube2.shape[0]
    new_b_size = cube2.shape[1]

    cube2_allfreqs = np.zeros((nfreqs, int(new_a_size), int(new_b_size))).astype(
        np.float32
    )
    # print(ska_alpha,ska_flux)
    for ff in range(nfreqs):
        norm = (
            ska_flux * (freqs[ff] / ska_freqmin) ** ska_alpha
        )  # flux mormalization for frequency freq
        cube2_allfreqs[ff] = cube2 * norm
        # print(norm,freqs[ff],np.sum(cube2_allfreqs[ff]),ska_alpha)
   

    if doplot:
        plt.subplot(121)
        plt.imshow(np.sum(cube, axis=0))
        plt.subplot(122)
        plt.imshow(np.sum(cube2, axis=0))
        plt.show()

    cube3 = cube2_allfreqs  # this will be the smoothed version

    logging.info("Rotating to postion angle: %f", ska_PA)
    # ~~~~~~~~~~~~ rotate to PA ~~~~~~~~~~~~~~~~~~~~~~

    cube4 = scipy.ndimage.rotate(
        cube3, ska_PA - 90, axes=(1, 2), reshape=False
    )  # subract 90 to give PA anti-clockwise from North

    logging.info("Final shape of subcube %s", cube4.shape)
    np.putmask(cube4, cube4 < 0, 0)

 

    
  


    return (
        cube4,
        cube_name.split('/')[-1],
        is_unresolved,
        ska_flux,
    )  # sources are now centred on dynamical centre of galaxy (according to crpix values in original_blanked_cubes) so dont need to pass crpix values (assuming centred)


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(sys.argv[1])



      


 
