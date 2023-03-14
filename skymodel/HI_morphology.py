# using morphological parameters from catalogue to generate an image: HI sources

import configparser
import logging
import sys
import time
import numpy as np
import scipy
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from matplotlib import pyplot as plt
from skimage import transform
from skimage.filters import gaussian

import skymodel.skymodel_tools as tools


deg2rad = (2 * np.pi) / 360
rad2arcsec = 206265
V_turb = 90  # from ALFALFA data, km/s


def make_cube(
    config,
    i,
    ska_M_HI,
    ska_z,
    ska_incl,
    ska_PA,
    ska_dx,
    dnu,
    psf_maj,
    radioclass,
):

    # set up cosmology
    mainlog = logging.getLogger("main%d" % i)

    logging.root.setLevel(logging.DEBUG)
    mainlog.info("make_cube%s" % i)
    H = config.getfloat("cosmology", "H")
    M = config.getfloat("cosmology", "M")
    L = config.getfloat("cosmology", "L")


    cosmo = LambdaCDM(H0=H, Om0=M, Ode0=L)

    # HI line rest freq
    rest_freq = config.getfloat("observation", "rest_freq")

    # ~~~~~~~~~~~~~~~~ match the TRECS catalogue source with an atlas sample ~~~~~~~~~~~~~~~~~~~~~

    # this can be done by using either the measured LASs or the MHI-DHI correlation
    # currently uses measured values

    doplot = config.getboolean("pipeline", "doplot")

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
    # make M_HI ratios square rooted, since D_HI ~ M_HI**0.5
    # have experimented with the relative ratios


    atlas_M_HI = 10 ** prepared_cubes[:, 7].astype(np.float)
    atlas_incl = prepared_cubes[:, 4].astype(np.float) * deg2rad

    logging.info("Selecting atlas source")
    logging.info("SKA properties to match: ")
    logging.info("MHI: %f", ska_M_HI)
    logging.info("Inclination: %f", ska_incl)
    logging.info("Radio class: %s", radioclass)


    if radioclass == 2.0:
        pass


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

    distance_to_sample = np.hypot(
        normalised_ska_incl - normalised_atlas_incl,
        normalised_ska_M_HI - normalised_atlas_M_HI,
    )
    distance_min_args = np.argsort(np.abs(distance_to_sample))
    np.random.seed(i)

    while True:
        distance_min_arg = distance_min_args[np.random.randint(3)]
        if distance_min_arg == 22:          
            if ska_incl / np.pi * 180 < 10:
                break
            else:
                continue
        else:
            break

    logging.info("Chosen atlas source: %s", prepared_cubes[distance_min_arg, 9])
    logging.info("Chosen atlas source MHI: %f", atlas_M_HI[distance_min_arg])
    logging.info("Chosen atlas source: %f", atlas_incl[distance_min_arg])

    # ~~~~~~~~~~~~~~~~~ load data and properties of matching source ~~~~~~~~~~~~~~~~~~

    i = prepared_cubes[distance_min_arg, :]

    # get atlas source properties from metadata

    atlas_incl = (float(i[4]) / 360.0) * 2 * np.pi
    atlas_M_HI = 10 ** float(i[7])
    atlas_LAS_maj = 10 ** float(i[10])  # kpc
    atlas_LAS_min = 10 ** float(i[11])  # kpc
    atlas_w20 = float(i[12])

    # get datacube and header properties
    cube_name = datacube_dir + prepared_dir + "{}cr_rotated_shrunk.fits".format((i[8]))

    cube_fits = fits.open(cube_name)
    cube = cube_fits[0].data

    # get a summed spectrum so that unresolved sources can be created from scratch (not currently needed)
    cube_summed_spectrum = np.sum(np.sum(cube, axis=1), axis=1)

    logging.info("Atlas sample shape: %s", cube.shape)

    if doplot:
        plt.scatter(np.arange(len(cube_summed_spectrum)), cube_summed_spectrum)
        plt.show()

    dz = np.abs(cube_fits[0].header["CDELT3"])  # m/s

    # ~~~~~~~~~~~~~~~~~~~~~ first calculate scalings and determine if resolved ~~~~~~~~~~~~~

    # Angular diameter distance in Mpc at a given redshift.
    ska_D_A = cosmo.angular_diameter_distance(ska_z).value  # Mpc
    ska_D_L = cosmo.luminosity_distance(ska_z).value  # Mpc

    # optical velocity sampling dZ at given redshift (freq sampling dnu is constant)
    central_freq = rest_freq / (ska_z + 1)
    ska_dz = tools.get_spectral_sampling(
        config, rest_freq, central_freq, dnu
    )  #  optical velocity interval
    logging.info("dZ at cat source redshift: %f", ska_dz)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ get scalings according to mass change ~~~~~~~~~~~~~~~~~~~~

    log10_atlas_D_HI = (0.51 * np.log10(atlas_M_HI)) - 3.32
    log10_ska_D_HI = (0.51 * np.log10(ska_M_HI)) - 3.32

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
    )  

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

    # b ratios
    # use cos**2(i) = (b/a)**2 - alpha**2 / 1 - (alpha)**2
    alpha = 0.2  # spiral galaxies; ellipticals = 0.5
    old_b_over_a = np.sqrt(((np.cos(atlas_incl) ** 2) * (1 - alpha ** 2)) + (alpha ** 2))
    new_b_over_a = np.sqrt(((np.cos(ska_incl) ** 2) * (1 - alpha ** 2)) + (alpha ** 2))
    incl_scale = new_b_over_a / old_b_over_a
    V_turb = 40  # this is lower than the mean alfalfa value, to avoid dampening the scaling between peaks in velocity too much
    if atlas_w20 < V_turb:
        V_turb = atlas_w20
    incl_v_scale = (
        np.sqrt(
            V_turb ** 2
            + (
                np.sqrt(atlas_w20 ** 2 - V_turb ** 2)
                * (np.sin(ska_incl) / np.sin(atlas_incl))
            )
            ** 2
        )
        / atlas_w20
    )

    incl_v_scale_no_vt = np.sin(ska_incl) / np.sin(atlas_incl)

    logging.info("Scalings due to inclination change:")
    logging.info("Scale factor along b axis: %f", incl_scale)
    logging.info("Velocity scale factor %f: ", incl_v_scale)

    # ~~~~~~~~~~~~  get the scalings for new z ~~~~~~~~~~~

    logging.info("SKA D_A, Mpc: %f", ska_D_A)  # Mpc
    logging.info("SKA D_L, Mpc: %f", ska_D_L)  # Mpc

    atlas_physical_pixel_size = (
        1000 * 1e-6
    )  # Mpc # this should come from fits header - put in

    # angular size is in radians
    # physical size and D_A in Mpc

    angular_pixel_scale_at_given_redshift = (
        atlas_physical_pixel_size / ska_D_A
    ) * rad2arcsec

    logging.info(
        "Angular scale (of a 100Mpc pixel) at given redshift, arcsec: %f",
        angular_pixel_scale_at_given_redshift,
    )

    # these are scalings due to pixel size changes (so are inverted wrt object size changes)
    # v resizing is due to change in v resolution, which varies with z
    redshift_scale = angular_pixel_scale_at_given_redshift / ska_dx
    redshift_v_scale = (dz * (1 + ska_z)) / ska_dz  # 1+z converts from restdV to optdV
    # v resizing is due to change in v resolution, which varies with z

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
        print("values incl", atlas_incl / np.pi * 180, ska_incl / np.pi * 180)
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
        new_a_size = np.ceil((cube.shape[1] * total_a_factor))

    while np.ceil(new_b_size) % 2 == 0:
        # a bit of padding on both sides to leave the image centered as before but to get an odd number of pixels in the end
        npad = int(
            1.0 / total_b_factor / 2.0
        )  # pixels in original size to correspond to 0.5 pixel in the new image  
        cube = np.pad(
            cube, ((0, 0), (0, 0), (npad, npad)), "constant", constant_values=(0)
        )
        new_b_size = np.ceil((cube.shape[2] * total_b_factor))

    cube2 = transform.resize(
        cube,
        (np.ceil(new_v_size), np.ceil(new_a_size), np.ceil(new_b_size)),
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    )
    cube_summed_spectrum = transform.resize(
        cube_summed_spectrum,
        (np.ceil(new_v_size), 1),
        order=1,
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

    # convert size scalings from kpc to pixels first (for determining whether resolved)
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
    cube3 = np.pad(
        cube2,
        ((0, 0), (padding, padding), (padding, padding)),
        "constant",
        constant_values=(0),
    ) 

    FWHM = psf_maj / ska_dx
    sigma2FWHM = 2 * np.sqrt(2 * np.log(2))
    Gaus_sigma = FWHM / sigma2FWHM
    cube4 = gaussian(cube3, sigma=[1, Gaus_sigma, Gaus_sigma])

    # ~~~~~~~~~~~~ rotate to PA ~~~~~~~~~~~~~~~~~~~~~~

    # scipy rotates in anti-clockwise direction
    cube5 = scipy.ndimage.rotate(cube4, ska_PA, axes=(2, 1), reshape=False)
    # angle defined anti-clockwise from North

    logging.info("Rotating to postion angle: %f", ska_PA)

    if "G" in i[8]:
        ska_PA = 180 - ska_PA  # fix due to 'G' atlas sources starting w receeding major axis at 180
    elif "H" in i[8]:
        ska_PA = -ska_PA

    logging.info("Corrected postion angle: %f", ska_PA)
    if ska_PA < 0:
        ska_PA += 360
    if ska_PA > 360:
        ska_PA -= 360

    # ~~~~~~~~~~~~~ normalise to correct flux density ~~~~~~~~~~~~~~~~

    # normalise flux scale according to HI mass and z (d) (blanking subtraction of orginal cubes means that flux in those is not physically correct)
    # using M_HI = 2.36 x 10**5 D**2 * integration of F dz; D is luminosity distance in Mpc, M_HI in solar masses, dz in km/s
    # see Duffy 2012 MNRAS 426, 3385, Roberts 1975
    flux_integral = ska_M_HI / (
        49.8 * ska_D_L ** 2
    )  # (Duffy+12, converted from below via c amd rest freq.. sampling is now fixed freq, should give same result? yes)
    # flux_integral = ska_M_HI/(2.36e5*ska_D_L**2) # equal to int S dV, units jy-km/s

    # convert from Jy-km/s to Jy by dividing by channel width in km/s
    ska_summed_flux = flux_integral / (
        dnu
    )  # this enables normalisation to correct Jy value per velocity slice
    logging.info("SKA summed flux from mass: %f", ska_summed_flux)
    cube5_ska_summed_flux = np.sum(cube5)
    logging.info("Subcube summed flux before normalisation: %f", cube5_ska_summed_flux)

    norm_factor = ska_summed_flux / cube5_ska_summed_flux
    cube5 *= norm_factor

    # normalise to Jy/beam instead of Jy/pix
    beam_area_in_pix = (np.pi * (FWHM) * (FWHM)) / (4 * np.log(2))
    logging.info("Beam area, pixels: %f", beam_area_in_pix)
    cube5_ska_summed_flux = np.sum(cube5)
    logging.info("Subcube flux after normalisation for mass: %f", cube5_ska_summed_flux)
    cube5 *= beam_area_in_pix
    cube5_ska_summed_flux = np.sum(cube5)
    logging.info(
        "Subcube flux after normalisation for pixels of beam: %f", cube5_ska_summed_flux
    )

    logging.info("Final shape of subcube %s", cube5.shape)
    np.putmask(cube5, cube5 < 0, 0)

    if doplot:

        plt.figure(10)

        plt.subplot(231)
        plt.imshow(np.sum(cube, axis=0))
        plt.title("original")
        plt.subplot(232)
        plt.imshow(np.sum(cube2, axis=0))
        plt.title("incl&mass")
        plt.subplot(233)
        plt.imshow(np.sum(cube3, axis=0))
        plt.title("redshift ")
        plt.subplot(234)
        plt.imshow(np.sum(cube4, axis=0))
        plt.title("mild conv")
        plt.subplot(235)
        plt.imshow(np.sum(cube5, axis=0))
        plt.title("rot %f deg" % ska_PA)

        plt.show()

    if doplot:
        plt.clf()
        plt.subplot(221)
        plt.imshow(np.sum(cube, axis=1))
        plt.subplot(222)
        plt.imshow(np.sum(cube, axis=2))
        plt.subplot(223)
        plt.imshow(np.sum(cube4, axis=1))
        plt.subplot(224)
        plt.imshow(np.sum(cube4, axis=2))
        plt.show()

    return (
        cube5,
        ska_D_HI_kpc_direct,
        ska_D_HI_arcsec_direct,
        i[8],
        flux_integral,
        ska_w20,
        MHI_incl_scale,
        MHI_incl_v_scale,
        ska_PA,
    ) 


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(sys.argv[1])
