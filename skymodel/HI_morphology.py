# using morphological parameters from catalogue to generate an image: HI sources


# to do

#

# check how many sources up to z = 0.5, resolved vs unresolved - what area of sky?
# normalise ---after--- rotation
# work out if unreolved first
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
import scipy
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter, zoom
from scipy.optimize import curve_fit
from skimage import transform
from skimage.filters import gaussian

import skymodel.skymodel_tools as tools

t0 = time.time()


deg2rad = (2 * np.pi) / 360
rad2arcsec = 206265
M_HI2M_dyn = 10  # (from Robert's fit to ALFALFA data)
V_turb = 90  # (from Robert's fit to ALFALFA data)
one_solar_mass = 1.989e30
kpc2m = 3.086e19
# atlas_physical_pixel_size = 100*1e-6 # Mpc this has changed in function
interp_order = 1
V_turb = 90  # km/s


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


def make_cube(
    config,
    i,
    ska_M_HI,
    ska_Mh,
    ska_z,
    ska_incl,
    ska_PA,
    ska_dx,
    dnu,
    psf_maj,
    radioclass,
    optclass,
):

    # set up cosmology
    mainlog = logging.getLogger("main%d" % i)

    logging.root.setLevel(logging.DEBUG)
    mainlog.info("make_cube%s" % i)
    H = config.getfloat("cosmology", "H")
    M = config.getfloat("cosmology", "M")
    L = config.getfloat("cosmology", "L")
    c = config.getfloat("cosmology", "c")
    G = config.getfloat("cosmology", "G")

    cosmo = LambdaCDM(H0=H, Om0=M, Ode0=L)

    # HI line rest freq

    rest_freq = config.getfloat("observation", "rest_freq")

    # ~~~~~~~~~~~~~~~~` match the TRECS catalogue source with an atlas sample ~~~~~~~~~~~~~~~~~~~~~

    # this can be done by using either the measured LASs or the MHI-DHI correlation
    # currently uses measured values

    doplot = config.getboolean("pipeline", "doplot")

    # test git
    ska_incl = ska_incl * deg2rad
    base_dir = config.get("pipeline", "base_dir")
    data_path = config.get("pipeline", "data_path")
    datacube_dir = base_dir + data_path + config.get("pipeline", "datacube_dir")
    prepared_dir = config.get("pipeline", "prepared_dir")
    prepared_metadata = config.get("pipeline", "prepared_metadata")
    prepared_cubes = np.loadtxt(datacube_dir + prepared_metadata, dtype="str")

    # find closest-matching atlas cube
    # properties that need matching:
    # keep incl ratios linear since v depends on sin i but b depends on ~cos a
    #  make M_HI ratios square rooted, since D_HI ~ M_HI**0.5
    # have experimented with the relative ratios

    atlas_names = prepared_cubes[:, 8]
    atlas_M_HI = 10 ** prepared_cubes[:, 7].astype(np.float)
    atlas_incl = prepared_cubes[:, 4].astype(np.float) * deg2rad
    atlas_class = prepared_cubes[:, 1].astype(np.float)

    ska_match_properties = np.array([ska_M_HI, ska_incl])
    atlas_match_properties = np.vstack((atlas_M_HI, atlas_incl)).T.astype(np.float)

    logging.info("Selecting atlas source")
    logging.info("SKA properties to match: ")
    logging.info("MHI: %f", ska_M_HI)
    logging.info("Inclination: %f", ska_incl)
    logging.info("Radio class: %s", radioclass)

    # print ('radioclass: ', radioclass)
    # print ('optclass: ', optclass)

    if radioclass == 2.0:
        pass

    """
    HI sources all have OptClass = 2; use RadioClass values to select morphology of HI atlas source. 
    Some RadioClass are -100 (where there is no radio counterpart).. think about the hosts for these (radio quiet)
    """

    normalised_atlas_M_HI = (np.log10(atlas_M_HI) - np.min(np.log10(atlas_M_HI))) / (
        np.max(np.log10(atlas_M_HI)) - np.min(np.log10(atlas_M_HI))
    )
    normalised_atlas_incl = (np.cos(atlas_incl) - np.min(np.cos(atlas_incl))) / (
        np.max(np.cos(atlas_incl)) - np.min(np.cos(atlas_incl))
    )

    normalised_ska_M_HI = (np.log10(ska_M_HI) - np.min(np.log10(atlas_M_HI))) / (
        np.max(np.log10(atlas_M_HI)) - np.min(np.log10(atlas_M_HI))
    )
    normalised_ska_incl = (np.cos(ska_incl) - np.min(np.cos(atlas_incl))) / (
        np.max(np.cos(atlas_incl)) - np.min(np.cos(atlas_incl))
    )

    # should no HI sources go into AGN-- agn sources still have hi_pred

    distance_to_sample = np.hypot(
        normalised_ska_incl - normalised_atlas_incl,
        normalised_ska_M_HI - normalised_atlas_M_HI,
    )
    distance_min_args = np.argsort(np.abs(distance_to_sample))
    np.random.seed(i)

    while True:
        distance_min_arg = distance_min_args[np.random.randint(3)]

        if distance_min_arg == 22:
            #  print ('in while')

            if ska_incl / np.pi * 180 < 10:
                #  print ('in while if')
                break
            else:
                continue
        else:
            break

    """
    distance_min_args_cut = distance_min_args[distance_to_sample[distance_min_args]<0.0]


    np.random.seed(i) 
    if len(distance_min_args_cut):         
        distance_min_arg = distance_min_args[np.random.randint(len(distance_min_args_cut))]
    else:
        distance_min_arg = distance_min_args[0]

    """

    logging.info("Chosen atlas source: %s", prepared_cubes[distance_min_arg, 9])
    logging.info("Chosen atlas source MHI: %f", atlas_M_HI[distance_min_arg])
    logging.info("Chosen atlas source: %f", atlas_incl[distance_min_arg])

    # ~~~~~~~~~~~~~~~~~ load data and properties of matching source ~~~~~~~~~~~~~~~~~~

    i = prepared_cubes[distance_min_arg, :]

    # return np.array([0,0]), 0, 0, 0, 0, i[8], 0, 0, 0  # sources are now centred on dynamical centre of galaxy (according to crpix values in original_blanked_cubes) so dont need to pass crpix values (assuming centred)

    # get atlas source properties from metadata
    distance = float(
        i[2]
    )  # Mpc # this is irrelevant for data since flux scale incorrect and pixels resized, but may be relevabt for derived proprties
    sys_vel = (
        float(i[3]) * 1000
    )  # m/s, should be similar to rec vel but will include pec vel
    incl = (float(i[4]) / 360.0) * 2 * np.pi
    rot_vel = float(
        i[5]
    )  # this is corrected for inclination, according to hyperleda from which it was taken
    pa = (
        float(i[6]) * deg2rad
    )  # make a new column for this : pas are in hyperleda thing
    atlas_M_HI = 10 ** float(i[7])
    atlas_LAS_maj = 10 ** float(i[10])  # kpc
    atlas_LAS_min = 10 ** float(i[11])  # kpc
    atlas_w20 = float(i[12])

    # get datacube and header properties
    cube_name = datacube_dir + prepared_dir + "{}cr_rotated_shrunk.fits".format((i[8]))
    print(cube_name)
    cube_fits = fits.open(cube_name)
    cube = cube_fits[0].data
    cube_summed_flux = np.sum(cube)
    # get a summed spectrum so that unresolved sources can be created from scratch - to be replaced with ALFALFA spectra
    cube_summed_spectrum = np.sum(np.sum(cube, axis=1), axis=1)

    logging.info("Atlas sample shape: %s", cube.shape)

    if doplot:
        plt.scatter(np.arange(len(cube_summed_spectrum)), cube_summed_spectrum)
        plt.show()

    dx = np.abs(cube_fits[0].header["CDELT1"] * 3600)  # arcsec
    #  need to change to:  atlas_physical_pixel_size = 1000*1e-6 # Mpc # this should come from fits header - put in
    dy = np.abs(cube_fits[0].header["CDELT2"] * 3600)  # arcsec
    dz = np.abs(cube_fits[0].header["CDELT3"])  # m/s

    velocities = cube_fits[0].header["CRVAL3"] + (np.arange(cube.shape[0]) * dz)
    atlas_bmaj = cube_fits[0].header["BMAJ"] * 3600
    atlas_bmin = cube_fits[0].header["BMIN"] * 3600
    atlas_bpa = cube_fits[0].header[
        "BPA"
    ]  # degrees? this will be modified when when rotating in prep
    sys_vel_pix = cube_fits[0].header["CRPIX3"]

    """
    # calculate HI mass from cube using  mass to flux: M_HI = 2.36 x 10**5 D**2 * integration of F dv, D is distance in Mpc (distance in metatdata)
    atlas_M_HI = 2.36*10**5 *(distance**2)*np.sum(cube)*(dz/1000.) # solar masses
    print ('M_HI: ', atlas_M_HI )
    print ('dz: ', dz)
    solid_beam_angle = (np.pi*(atlas_bmaj/206265)*(atlas_bmin/206265))/(4*np.log(2))
    pix_area_over_beam_area = ((dx/206265)*(dy/206265))/solid_beam_angle # gives number of beams per pixel; need to do it on blanked sources to get accurate??
    # also gets rid of convolution artefacts
    print ('M_HI: ', atlas_M_HI )
    atlas_M_HI*=pix_area_over_beam_area
    print ('M_HI: ', atlas_M_HI )
    #this gets same number as https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_masscalculator.php
    # does not match recorded masses, but will be overestimate due to non-blanking 

    """

    t1 = time.time() - t0

    # ~~~~~~~~~~~~~~~~~~~~~ first calculate scalings and determine if resolved ~~~~~~~~~~~~~

    # Angular diameter distance in Mpc at a given redshift.

    ska_D_A = cosmo.angular_diameter_distance(ska_z).value  # Mpc

    ska_D_L = cosmo.luminosity_distance(ska_z).value  # Mpc

    # optical velocity sampling dZ at given redshift (freq sampling dnu is constant)
    # resampling to this interval should allow the subcube to be place at the given freq with the correct
    # velocity sampling for that freq. ???. The linear sampling of the subcube is ok, since the range of velocity
    # per subcube is relatively small

    central_freq = rest_freq / (ska_z + 1)
    ska_dz = tools.get_spectral_sampling(
        config, rest_freq, central_freq, dnu
    )  #  optical velocity interval
    print("ska_dz", ska_dz)
    logging.info("dZ at cat source redshift: %f", ska_dz)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ get scalings according to mass change ~~~~~~~~~~~~~~~~~~~~

    # first calculate diameter from  log D_HI = 0.51 log M_HI -3.32  D= kpc, MH1 = Msolar, Broeils and Rhee 1997

    # notes that the uncertainties are log10 MHI = (1.96+- 0.04)log DHI +(6.52+-0.06)
    # need to measure for our own objects to test it satisfies this relation
    # D_HI is defined at a surface density of 1 Mpc−2 and apparently independent of environment and other factors

    # have: rot_vel**2 = GM_dyn/r and M_dyn/M_HI = 10, V_turb = 90
    # assume same M_dyn/M_HI fraction for ska_M_HI:
    # rot_vel**2 = G * 10 * M_HI / r
    ###this might all be different in ellipticals (see gas poor galaxies on caltech page) - need to look through galaxy types in atlas

    # check atlas values to test: inputted M_HI and got out vrot, numbers agree:
    # atlas_M_HI = 10**9.23 # if value from metadata is used then get correct answer for vrot
    log10_atlas_D_HI = (0.51 * np.log10(atlas_M_HI)) - 3.32
    log10_ska_D_HI = (0.51 * np.log10(ska_M_HI)) - 3.32
    # this gets the values of the atlas and ska sizes from their masses, according to the relation
    # if the actual measured atlas size does not fall on the relation, then it could be said to fall on a different relation
    # we could say that the ska size also falls on this other relation
    # the ratio of sizes obtained from the ratio the relation as a function of mass is independent of the factor 10**-3.32
    # This factor contains the density of the sources; less dense types (SA, SAb) will have a higher (lower absolute) value of this factor
    # this would produce an offset, parallel, relation in log log space.
    #  A change in the exponent (0.51) could come about if the source was more spherical or less constant in density
    # Going to assume that the M**0.51 holds for all sources, and factor -3.32 could vary.

    ska_D_HI_kpc = (
        10 ** log10_ska_D_HI
    )  # kpc; retain for comaprison with D_HI measured from atlas cube
    ska_D_HI_arcsec = (
        ska_D_HI_kpc / (ska_D_A * 1000)
    ) * rad2arcsec  # convert D_A to kpc

    logging.info("SKA D_HI from Broeils and Rhee, kpc: %f", ska_D_HI_kpc)
    logging.info("Atlas D_HI from Broeils and Rhee, kpc: %f", 10 ** log10_atlas_D_HI)
    logging.info(
        "Atlas D_HI directly measured, kpc: %f", atlas_LAS_maj
    )  # ususally in agreement, sometimes up to a factor of ~1.7 out
    # it's not within 1 sigma but is within 2 sigma in the worst cases (eg the lowest z cat soucrse)
    # it is simply due to the morphology of the atlas source - some souces have big holes inbetween core and outer ring
    # see notes above - types SA and SAb are seen in Broeils and Rhee to have lower density.

    # need to convert D_HI into a value of r, or just scale:
    # v2 = v1* sqrt(M2/M1 * D1/D2)
    # v should go as ~M**0.25 based on combining equations
    # checks out according to sample test ie. M ratio 0.26, v ratio 0.26**0.25 =  0.72
    # i.e ~quartering mass makes size shrink by ~half and rot vel shrinks to 0.72%
    # mass times 8; size times ~4; rot vel time ~2

    M_HI_factor = ska_M_HI / atlas_M_HI
    MHI_scale = 10 ** (log10_ska_D_HI) / 10 ** (
        log10_atlas_D_HI
    )  # is equal to M_HI_factor**0.51
    MHI_v_scale = np.sqrt(M_HI_factor / MHI_scale)

    logging.info("Scalings due to mass change:")
    logging.info("MHI factor: %f", M_HI_factor)
    logging.info("DHI factor (scale along both spatial axes): %f", MHI_scale)
    logging.info("Velocity scale factor: %f", MHI_v_scale)

    # ~~~~~~~~~~~~~~~~~~~ get scalings according to incl change ~~~~~~~~~~~~~~~~~~

    ################ b ratios

    # use cos**2(i) = (b/a)**2 - alpha**2 / 1 - (alpha)**2

    alpha = 0.2  # spiral galaxies; ellipticals = 0.5
    old_b_over_a = np.sqrt(((np.cos(incl) ** 2) * (1 - alpha ** 2)) + (alpha ** 2))
    new_b_over_a = np.sqrt(((np.cos(ska_incl) ** 2) * (1 - alpha ** 2)) + (alpha ** 2))

    # a is a constant with varying incl; stretch factor is then (b2/a2)/(b1/a1) (a s cancel)

    incl_scale = new_b_over_a / old_b_over_a

    # rot_vel = rad_vel(r)/sin(i)
    # contributions to w20 are added in quadature:
    # w20 = sqrt(V_T**2 + 2V_rad**2)
    # new_w20 = old_w20*(new_w20/old_w20)

    # (new_w20/old_w20) = sqrt(v_turb**2+(2rot_vel*sin(new_i))**2)/sqrt(v_turb**2+(2rot_vel*sin(old_i))**2)

    # 2V_rot*sin(old_i) = sqrt(old_w20**2-V_turb**2)
    # atlas_V_rot = sqrt(old_w20**2-V_turb**2)/(2*sin(old_i))

    V_turb = 40  # this is lower than the mean alfalfa value, to avoid dampening the scaling between peaks in velocity too much
    if atlas_w20 < V_turb:
        V_turb = atlas_w20
    incl_v_scale = (
        np.sqrt(
            V_turb ** 2
            + (
                np.sqrt(atlas_w20 ** 2 - V_turb ** 2)
                * (np.sin(ska_incl) / np.sin(incl))
            )
            ** 2
        )
        / atlas_w20
    )

    incl_v_scale_no_vt = np.sin(ska_incl) / np.sin(incl)
    print("atlas_w20", atlas_w20)
    print(incl_v_scale, incl_v_scale_no_vt)

    logging.info("Scalings due to inclination change:")
    logging.info("Scale factor along b axis: %f", incl_scale)
    logging.info("Velocity scale factor %f: ", incl_v_scale)

    # ~~~~~~~~~~~~  get the scalings for new z ~~~~~~~~~~~

    logging.info("SKA D_A, Mpc: %f", ska_D_A)  # Mpc
    logging.info("SKA D_L, Mpc: %f", ska_D_L)  # Mpc

    atlas_physical_pixel_size = (
        1000 * 1e-6
    )  # Mpc # this should come from fits header - put in

    # D_A = physical_size/angular_size
    #  angular size is in radians
    # physical size and D_A in Mpc

    angular_pixel_scale_at_given_redshift = (
        atlas_physical_pixel_size / ska_D_A
    ) * rad2arcsec

    logging.info(
        "Angular scale (of a 100Mpc pixel) at given redshift, arcsec: %f",
        angular_pixel_scale_at_given_redshift,
    )

    # these are scalings due to pixel size changes (so are inverted wrt object size changes)
    redshift_scale = angular_pixel_scale_at_given_redshift / ska_dx

    redshift_v_scale = (dz * (1 + ska_z)) / ska_dz  # 1+z converts from restdV to optdV
    # v resizing is due to change in v resolution - which does vary with z

    logging.info("Scalings due to moving to redshift:")
    logging.info("Scale factor along both spatial axes: %f", redshift_scale)
    logging.info(
        "Velocity scale factor (due to change of v resolution only): %f",
        redshift_v_scale,
    )

    # ~~~~~~~~~~~~~~~ combine the scalings and do the resizing ~~~~~~~~~~~~~~~~~~~

    MHI_incl_scale = MHI_scale * incl_scale
    MHI_incl_v_scale = MHI_v_scale * incl_v_scale

    if MHI_incl_v_scale > 10:
        print("scalings: ", MHI_incl_scale, MHI_incl_v_scale)
        print("values MHI:", atlas_M_HI, ska_M_HI)
        print("values incl", incl / np.pi * 180, ska_incl / np.pi * 180)
        print("distance:", distance_to_sample[distance_min_arg])
        doplot = 1

    # axes order: v,a,b
    total_v_factor = MHI_v_scale * incl_v_scale * redshift_v_scale
    total_b_factor = MHI_scale * incl_scale * redshift_scale
    total_a_factor = MHI_scale * redshift_scale
    new_v_size = cube.shape[0] * total_v_factor
    new_a_size = cube.shape[1] * total_a_factor
    new_b_size = cube.shape[2] * total_b_factor
    logging.info(
        "Combined scale factors due to mass and inclination: %f, %f, %f",
        MHI_v_scale * incl_v_scale,
        MHI_scale * incl_scale,
        MHI_scale,
    )
    logging.info(
        "Combined scale factors due to mass and inclination and z: %f, %f, %f",
        total_v_factor,
        total_b_factor,
        total_a_factor,
    )

    while np.ceil(new_a_size) % 2 == 0:

        # a bit of padding on both sides to leave the image centered as before but to get an odd number of pixels in the end
        npad = int(
            (1.0 / total_a_factor / 2.0)
        )  # pixels in original size to correspond to 0.5 pixel in the new image

        cube = np.pad(
            cube, ((0, 0), (npad, npad), (0, 0)), "constant", constant_values=(0)
        )
        print((cube.shape[1] * total_a_factor))
        new_a_size = np.ceil((cube.shape[1] * total_a_factor))

    while np.ceil(new_b_size) % 2 == 0:

        # a bit of padding on both sides to leave the image centered as before but to get an odd number of pixels in the end
        npad = int(
            1.0 / total_b_factor / 2.0
        )  # pixels in original size to correspond to 0.5 pixel in the new image
        print("padding on both sizes of", npad)

        cube = np.pad(
            cube, ((0, 0), (0, 0), (npad, npad)), "constant", constant_values=(0)
        )

        new_b_size = np.ceil((cube.shape[2] * total_b_factor))

    # deal with nearly face-on sources
    ## V_turb = 90
    ## try_min_v_size = V_turb/ska_dz
    ## min_v_size = np.min([try_min_v_size,cube.shape[0]])#make sure we are not stretching if min size from rough calculation is actually larger than input cube size in v
    ## print ('new_v_size, try_min_v_size, min_v_size',new_v_size, try_min_v_size, min_v_size)
    ## new_v_size = np.max([new_v_size, min_v_size])# to prevent single-chaneel sources
    print("new_v_size", new_v_size)
    cube2 = transform.resize(
        cube,
        (np.ceil(new_v_size), np.ceil(new_a_size), np.ceil(new_b_size)),
        order=interp_order,
        preserve_range=True,
        anti_aliasing=True,
    )
    cube_summed_spectrum = transform.resize(
        cube_summed_spectrum,
        (np.ceil(new_v_size), 1),
        order=interp_order,
        preserve_range=True,
        anti_aliasing=True,
    )

    logging.info(
        "Subcube shape after scaling for mass and inclination and z: %s", cube2.shape
    )

    # scale physical properties
    ska_LAS_maj_kpc = atlas_LAS_maj * MHI_scale
    ska_LAS_min_kpc = atlas_LAS_min * MHI_scale * incl_scale

    ska_w20 = atlas_w20 * MHI_v_scale * incl_v_scale
    logging.info("SKA DHI after all physical scalings, kpc: %f", ska_LAS_maj_kpc)
    logging.info("SKA w20 after all physical scalings, km/s: %f", ska_w20)

    # ~~~~~~~~~~~~  work out sizes in arcsec, at new z ~~~~~~~~~~~

    # convert size scalings from kpc to pixels first (for determining if resolved)
    ska_LAS_maj_pix = (
        ska_LAS_maj_kpc / (atlas_physical_pixel_size * 1000)
    ) * redshift_scale  # 1000 converts to kpc
    ska_LAS_min_pix = (
        ska_LAS_min_kpc / (atlas_physical_pixel_size * 1000)
    ) * redshift_scale

    ska_D_HI_arcsec_direct = ska_LAS_maj_pix * ska_dx
    ska_D_HI_kpc_direct = ska_LAS_maj_kpc

    logging.info("SKA DHI after all scalings, arcsec: %f", ska_D_HI_arcsec_direct)
    logging.info("SKA DHI from Broeils and Rhee, arcsec: %f", ska_D_HI_arcsec)

    is_unresolved = 0

    if ska_LAS_maj_pix < 3.0:
        is_unresolved = 1

    # ~~~~~~~~~~~~~~~~ convolution ~~~~~~~~~~~~~~~~~

    padding = 5
    cube4 = np.pad(
        cube2,
        ((0, 0), (padding, padding), (padding, padding)),
        "constant",
        constant_values=(0),
    )  # need to sort this padding out

    FWHM = psf_maj / ska_dx
    sigma2FWHM = 2 * np.sqrt(2 * np.log(2))
    Gaus_sigma = FWHM / sigma2FWHM
    cube5 = gaussian(cube4, sigma=[1, Gaus_sigma, Gaus_sigma])

    # ~~~~~~~~~~~~ rotate to PA ~~~~~~~~~~~~~~~~~~~~~~

    # scipy rotates in anti-clockwise direction
    cube6 = scipy.ndimage.rotate(cube5, ska_PA, axes=(2, 1), reshape=False)
    # angle defins anti-clockwise from North

    logging.info("Rotating to postion angle: %f", ska_PA)

    if "G" in i[8]:
        ska_PA = 180 - ska_PA  # due to G sources starting w receeding major axis at 180

    elif "H" in i[8]:
        ska_PA = -ska_PA

    logging.info("Corrected postion angle: %f", ska_PA)
    if ska_PA < 0:
        ska_PA += 360
    if ska_PA > 360:
        ska_PA -= 360

    # ~~~~~~~~~~~~~ normalise to correct flux density ~~~~~~~~~~~~~~~~

    # normalise flux scale according to HI mass and z (d) (blanking subtraction of orginal cubes means that flux in those is not physically correct)
    # using M_HI = 2.36 x 10**5 D**2 * integration of F dz, D is distance in Mpc, M_HI in solar masses, dz in km/s
    # D is luminosity distance: see Duffy 2012 MNRAS 426, 3385, Roberts 1975
    flux_integral = ska_M_HI / (
        49.8 * ska_D_L ** 2
    )  # (Duffy+12, converted from below via c amd rest freq.. sampling is now fixed freq, should give same result? yes)
    # flux_integral = ska_M_HI/(2.36e5*ska_D_L**2) # equal to int S dV, units jy-km/s

    # convert from Jy-km/s to Jy by dividing by channel width in km/s
    ska_summed_flux = flux_integral / (
        dnu
    )  # this enables normalisation to correct Jy value per velocity slice
    logging.info("SKA summed flux from mass: %f", ska_summed_flux)
    cube6_ska_summed_flux = np.sum(cube6)
    logging.info("Subcube summed flux before normalisation: %f", cube6_ska_summed_flux)

    norm_factor = ska_summed_flux / cube6_ska_summed_flux
    cube6 *= norm_factor

    # normalise to Jy/beam instead of Jy/pix
    beam_area_in_pix = (np.pi * (FWHM) * (FWHM)) / (4 * np.log(2))
    logging.info("Beam area, pixels: %f", beam_area_in_pix)
    cube6_ska_summed_flux = np.sum(cube6)
    logging.info("Subcube flux after normalisation for mass: %f", cube6_ska_summed_flux)
    cube6 *= beam_area_in_pix
    cube6_ska_summed_flux = np.sum(cube6)
    logging.info(
        "Subcube flux after normalisation for pixels of beam: %f", cube6_ska_summed_flux
    )

    logging.info("Final shape of subcube %s", cube6.shape)
    np.putmask(cube6, cube6 < 0, 0)

    if doplot:

        plt.figure(10)

        plt.subplot(231)
        plt.imshow(np.sum(cube, axis=0))
        plt.title("original")
        plt.subplot(232)
        plt.imshow(np.sum(cube2, axis=0))
        plt.title("incl&mass")
        plt.subplot(233)
        plt.imshow(np.sum(cube4, axis=0))
        plt.title("redshift ")
        plt.subplot(234)
        plt.imshow(np.sum(cube5, axis=0))
        plt.title("mild conv")
        plt.subplot(235)
        plt.imshow(np.sum(cube6, axis=0))
        plt.title("rot %f deg" % ska_PA)

        plt.show()

    if doplot:
        plt.clf()

        plt.subplot(221)
        plt.imshow(np.sum(cube, axis=1))
        plt.subplot(222)
        plt.imshow(np.sum(cube, axis=2))
        plt.subplot(223)
        plt.imshow(np.sum(cube5, axis=1))
        plt.subplot(224)
        plt.imshow(np.sum(cube5, axis=2))
        plt.show()

    return (
        cube6,
        ska_D_HI_kpc_direct,
        ska_D_HI_arcsec_direct,
        ska_D_HI_kpc,
        ska_D_HI_arcsec,
        i[8],
        is_unresolved,
        flux_integral,
        ska_w20,
        MHI_incl_scale,
        MHI_incl_v_scale,
        ska_PA,
    )  # sources are now centred on dynamical centre of galaxy (according to crpix values in original_blanked_cubes) so dont need to pass crpix values (assuming centred)


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(sys.argv[1])
    make = 0
    if make == 1:

        """
        for i in cat:
            etc
        """

        # ska test params

        # ska params
        ska_bmaj = config.getfloat("skymodel", "simple_psf_maj")
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
