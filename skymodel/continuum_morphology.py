# using morphological parameters from catalogue to generate an image: HI sources


# to do

#

# check how many sources up to z = 0.5, resolved vs unresolved - what area of sky?
# normalise ---after--- rotation
# work out if unresolved first
# get unresolved profiles
# plot locations after cut


import configparser
import glob
import logging
import os
import sys
import time

import astropy
import galsim
import matplotlib
import numpy as np

# from skimage.filters import gaussian
import scipy

# from scipy.ndimage import zoom
import scipy.ndimage

# from scipy.signal import convolve as scipy_convolve
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D
from astropy.modeling.models import Sersic2D
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from skimage import transform

import skymodel.skymodel_tools as tools

t0 = time.time()


deg2rad = (2 * np.pi) / 360
rad2arcsec = 206265
M_HI2M_dyn = 10  # (from Robert's fit to ALFALFA data)
V_turb = 90  # (from Robert's fit to ALFALFA data)
one_solar_mass = 1.989e30
kpc2m = 3.086e19
atlas_physical_pixel_size = 100 * 1e-6  # Mpc
interp_order = 1


def fsigmoid(x, a, b, k):
    return a + (
        (k - a) / (1.0 + np.exp(-b * x))
    )  # could generalise with 1/nu as an index to denominator


def mxmul(a, b):
    output = np.zeros(4)
    output[0] = a[0] * b[0] + a[1] * b[2]
    output[1] = a[0] * b[1] + a[1] * b[3]
    output[2] = a[2] * b[0] + a[3] * b[2]
    output[3] = a[2] * b[1] + a[3] * b[3]
    return output


def mxinv(a):
    det = a[0] * a[3] - a[1] * a[2]
    output = np.array([a[3], -a[1], -a[2], a[0]]) / det
    return output


def mkgauss(naxes, pos, flux, fwhm, axrat=1.0, angle=0.0, ignore=4.0, dodist=False):

    # note that total flux = peak flux in a pixel * 1.1331*FWHM**2
    # angle is major axis East of North
    a = np.zeros(naxes[0] * naxes[1]).reshape(naxes[1], naxes[0])
    fwhm /= 1.66667
    if axrat == 1.0 and angle == 0.0:
        for i in range(naxes[1]):
            ydist = float(i) - pos[1]
            for j in range(naxes[0]):
                xdist = float(j) - pos[0]
                if xdist * xdist + ydist * ydist > ignore * ignore * fwhm * fwhm:
                    continue
                if not dodist:
                    a[i, j] = (
                        flux
                        * np.exp(-(xdist * xdist + ydist * ydist) / (fwhm * fwhm))
                        / (fwhm * fwhm * np.pi)
                    )
                else:
                    a[i, j] = np.hypot(xdist, ydist)
        return a
    sinth = np.sin(angle * np.pi / 180.0)
    costh = np.cos(angle * np.pi / 180.0)
    r = np.array([-sinth, costh, -costh, -sinth])
    rt = np.array([-sinth, -costh, costh, -sinth])
    sig = np.array([fwhm, 0.0, 0.0, fwhm * axrat])
    scr1 = mxmul(sig, r)
    scr2 = mxmul(rt, scr1)
    scr1 = mxinv(scr2)
    for i in range(naxes[1]):
        ydist = float(i) - pos[1]
        if abs(ydist) > ignore * fwhm:
            continue
        for j in range(naxes[0]):
            xdist = float(j) - pos[0]
            if abs(xdist) > ignore * fwhm:
                continue
            ex = scr1[0] * xdist + scr1[1] * ydist
            ey = scr1[2] * xdist + scr1[3] * ydist
            if not dodist:
                a[i, j] = (
                    (flux / axrat)
                    * np.exp(-(ex * ex + ey * ey))
                    / (fwhm * fwhm * np.pi)
                )
            else:
                a[i, j] = np.hypot(ex, ey) / 1.6666667

    return a


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

    #    ska_incl = ska_incl*deg2rad
    base_dir = config.get("pipeline", "base_dir")
    data_path = config.get("pipeline", "data_path")
    # datacube_dir = base_dir +data_path+config.get('pipeline', 'datacube_dir')
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
        )  # scale lenght to FHWM - true for SFGs
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

        # cut the image so it's not unnecessarily big
        # ix=0
        # newtot=1.
        # while (newtot ==1.):
        #    newtot=np.sum(cube2[ix:npixs-1-ix,ix:npixs-1-ix])
        #    print(ix,newtot)
        #    ix=ix+1

        # exit()
        # ix=ix-2

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
        print("This is a postage stamp")
        is_unresolved = 0
        prepared_cubes = np.loadtxt(datacube_dir + prepared_metadata, dtype="str")

        # find closest-matching atlas cube
        # properties that need matching:
        # keep incl ratios linear since v depends on sin i but b depends on ~cos a
        # make M_HI ratios square rooted, since D_HI ~ M_HI**0.5
        # have experimented with the relative ratios

        atlas_names = prepared_cubes[:, 0]
        atlas_ranid = prepared_cubes[:, 1].astype(np.float)

        print("Selecting atlas source")
        distance_to_sample = abs(atlas_ranid - ska_ranid)

        # print(ska_ranid)
        # print(np.min(atlas_ranid),np.max(atlas_ranid))
        # print(np.min(distance_to_sample),np.max(distance_to_sample))
        distance_min_arg = np.argmin(distance_to_sample)
        # print(distance_min_arg.shape)
        # print(distance_to_sample[distance_min_arg,0])

        # exit()
        logging.info("Chosen atlas source: %s", prepared_cubes[distance_min_arg, 0])

        # get datacube and header properties
        # cube_name = datacube_dir+prepared_dir+'{}cr_rotated_shrunk.fits'.format((i[9]))

        cube_name = prepared_cubes[distance_min_arg, 0]  # old format
        print(cube_name)

        # PH: temporary path fix
        #cube_name = cube_name.split(
        #    "/home/a.bonaldi/data-cold-for-backup/data_challenges/inputs/AGN_library/"
        #)[1]

        print("cube_name: ", cube_name)

        cube_fits = fits.open(cube_name)
        cube = cube_fits[0].data
        cube_summed_flux = np.sum(cube)

        logging.info("Atlas sample shape: %s", cube.shape)

        if doplot:
            plt.scatter(np.arange(len(cube_summed_spectrum)), cube_summed_spectrum)
            plt.show()

        dx = np.abs(cube_fits[0].header["CDELT1"] * 3600)  # arcsec  pix reso
        # dy = np.abs(cube_fits[0].header['CDELT2']*3600) # arcsec
        atlas_bmaj = cube_fits[0].header["SIZE"]

        #    atlas_bmin =  cube_fits[0].header['BMIN'] *3600
        #    atlas_bpa = cube_fits[0].header['BPA']  # degrees? this will be modified when when rotating in prep
        #    sys_vel_pix = cube_fits[0].header['CRPIX3']

        # c1 = np.abs(cube_fits[0].header['CRPIX1']*3600) # arcsec  pix reso
        # c2 = np.abs(cube_fits[0].header['CRPIX2']*3600) # arcsec
        # fn = np.abs(cube_fits[0].header['CRVAL3']) # frequency
        atlas_psf = np.abs(cube_fits[0].header["PSF"])  # arcsec
        atlas_bmaj_px = atlas_bmaj / dx
        zoom_factor = ska_size / ska_dx / atlas_bmaj_px

        logging.info("atlas agn size in arcsec %f", atlas_bmaj)
        logging.info("target agn size in arcsec %f", ska_size)
        logging.info("atlas agn size in pixels %f", atlas_bmaj_px)
        logging.info("target agn size in pixels %f", ska_size / ska_dx)
        logging.info("zoom factor %f", zoom_factor)

        # convolve postage stamp with beam to get to final resolution

        """
        psf_smooth=np.sqrt(np.square(psf_maj/dx)-np.square(atlas_psf/dx)) #additional smoothing PSF in units of the atlas pixel

    
        if (psf_smooth > 0. ):
            print('smoothing done on postage reso')
            print(psf_smooth)
            atlas_smoo = psf_smooth/sigma2FWHM
            atlas_kernel = Gaussian2DKernel(x_stddev=atlas_smoo)


#            padding_a=cube.shape[0]
#            padding_b=cube.shape[1]
            size_a=cube.shape[0]
            size_b=cube.shape[1]

            padding_a=int(np.ceil(size_a/2.))
            padding_b=int(np.ceil(size_b/2.))
            
            cube_padded = np.pad(cube, ((padding_a,padding_a), (padding_b,padding_b)), 'constant', constant_values=(0))

            astropy_conv = convolve(cube_padded, atlas_kernel,boundary='extend', normalize_kernel=True)
            
            #cube = astropy_conv[padding_a:2*padding_a,padding_b:2*padding_b]
            cube = astropy_conv[padding_a:padding_a+size_a,padding_b:padding_b+size_b]
            cube=cube/np.sum(cube)
            
            smoothing_done=1 #flag so that it is not smoothed later
        """

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

            print("padding on both sizes of", npad)
            cube = np.pad(
                cube, ((npad, npad), (npad, npad)), "constant", constant_values=(0)
            )

            new_a_size = np.ceil((cube.shape[0] * zoom_factor))
            new_b_size = np.ceil((cube.shape[1] * zoom_factor))

        cube2 = transform.resize(
            cube, (new_a_size, new_b_size), order=interp_order, preserve_range=True
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

        print("smoothing to be done in pixels", psf_smooth)
        if psf_smooth > 0.0:

            atlas_smoo = psf_smooth / sigma2FWHM
            print("sigma smoothing done on reduced beam", atlas_smoo)
            print("instead of", Gaus_sigma)
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
            cube2 = cube / np.sum(cube2)

        smoothing_done = 1  # flag so that it is not smoothed later (whether I did or not because in that case atlas was low reso)
        print(cube2.shape)

    if smoothing_done == 0:
        # convolution done on all sources - resolved and unresolved - check
        print("smoothing done on image reso")
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
    print("postage done")

    if doplot:
        plt.subplot(121)
        plt.imshow(np.sum(cube, axis=0))
        plt.subplot(122)
        plt.imshow(np.sum(cube2, axis=0))
        plt.show()

    cube5 = cube2_allfreqs  # this will be the smoothed version

    logging.info("Rotating to postion angle: %f", ska_PA)
    # ~~~~~~~~~~~~ rotate to PA ~~~~~~~~~~~~~~~~~~~~~~

    cube6 = scipy.ndimage.rotate(
        cube5, ska_PA - 90, axes=(1, 2), reshape=False
    )  # subract 90 to give PA anti-clockwise from North

    logging.info("Final shape of subcube %s", cube6.shape)
    np.putmask(cube6, cube6 < 0, 0)

    if doplot:
        plt.figure(10)
        plt.subplot(141)
        plt.imshow(np.sum(cube2, axis=0))
        plt.title("incl&mass")
        plt.subplot(142)
        plt.imshow(np.sum(cube4, axis=0))
        plt.title("redshift ")
        plt.subplot(143)
        plt.imshow(np.sum(cube5, axis=0))
        plt.title("mild conv")
        plt.subplot(144)
        plt.imshow(np.sum(cube6, axis=0))
        plt.title("rot %f deg" % ska_PA)
        plt.show()

    # print(np.sum(cube6[0]))
    # print(np.sum(cube6[1]))
    ### exit()
    print("shape:::::::::: ", cube6.shape, cube6.nbytes)
    #  cube6 = np.ones((5,100,100))
    return (
        cube6,
        cube_name,
        is_unresolved,
        ska_flux,
    )  # sources are now centred on dynamical centre of galaxy (according to crpix values in original_blanked_cubes) so dont need to pass crpix values (assuming centred)


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(sys.argv[1])

    analysis = 0
    prep = 1
    make = 0

    if analysis == 1:
        plot_properties(config)
    if prep == 1:
        prep_cube(config)
    if make == 1:

        """
        for i in cat:
            etc
        """

        # ska test params

        # ska params
        ska_bmaj = 5  # arcsec
        ska_bmin = ska_bmaj  # arcsec

        pixel_scale = ska_bmin / 2.0
        velo_resolution = 5000

        pixel_scale = ska_bmin / 3.0  # get these from input file
        velo_resolution = 5000  # m/s (ideally want 20 000 channels x 5 km/s)

        ska_z = 0.2
        ska_M_HI = 10 ** 9.3  # solarmass

        ska_Mh = ska_M_HI * 10  # log 10 solar
        # ska_env = ?  ska_env possibly provided

        ska_incl = (40 / 360.0) * 2 * np.pi

        ska_PA = 70  # degrees anti-clockwise from N-S (y axis)

        ska_radioclass = 1

        sub_cube, ska_D_HI, ska_D_HI_arcsec, atlas_source = make_cube(
            config,
            ska_M_HI,
            ska_Mh,
            ska_z,
            ska_incl,
            ska_radioclass,
            ska_PA,
            pixel_scale,
            velo_resolution,
        )

    # change in HI_mass

    # have/need:
    # mass vs rot vel
    # mass vs major axis
    # flux vs mass - tightly correlated (assumed optically thin)

# spatial scaling

# use log D_HI = 0.51 log M_HI -3.32 (Broeils and Rhee 1997)
