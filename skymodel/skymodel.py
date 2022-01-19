"""
Script to convert a T-RECS catalogue into a sky model FITS file.

Usage:
python skymodel.py example.ini

Author:
Philippa Hartley
p.hartley@skatelescope.org

Credit:
Ian Harrison's simuCLASS pipeline
ian.harrison-2@manchester.ac.uk

"""


import configparser
import inspect
import logging
import multiprocessing
import os
import pdb
import pickle
import sys
import time

import fitsio
import galsim

###from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uns
from astropy import wcs as ast_wcs
from astropy.coordinates import SkyCoord

# from spectral_cube import SpectralCube
from astropy.cosmology import LambdaCDM
from astropy.io import fits as astfits
from astropy.table import Table
from fitsio import FITS, FITSHDR
from numpy.core.defchararray import add as stradd
from numpy.core.defchararray import multiply as strmultiply
from scipy import random

import skymodel.skymodel_tools as tools
from skymodel.HI_morphology import make_cube
from skymodel.skymodel_tools import setup_wcs

# from primarybeam.primarybeam import *


tstart = time.time()


big_fft_params = galsim.GSParams(maximum_fft_size=84188)
arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value


def add_source(
    i,
    cat_gal,
    nobj,
    w_spectral,
    config,
    pixel_scale_str,
    dnu,
    psf_maj_arcsec,
    arr_dims,
    all_gals_fname,
    cat,
):

    logging.info(
        "..........Adding source {0} of {1} to skymodel..........".format(i + 1, nobj)
    )

    x, y, v = w_spectral.wcs_world2pix(
        cat_gal["RA"], cat_gal["Dec"], cat_gal["opt_vel"], 1
    )
    #    for three axes:   pixels = wcs.world_to_pixel(coord, 3000 * u.m / u.s)

    logging.info(
        "RA, Dec, opt_vel: %f %f %f  ",
        cat_gal["RA"],
        cat_gal["Dec"],
        cat_gal["opt_vel"],
    )
    logging.info("freq, z: %f %f", cat_gal["central_freq"], cat_gal["z"])
    logging.info("x, y, v: %f %f %f ", x, y, v)
    logging.info("HI_size from cat:  %f", cat_gal["HI_size"])

    # Create the sub-image for this galaxy

    #### get the sub cube from modified cube library sample:
    (
        sub_cube,
        D_HI,
        D_HI_arcsec,
        D_HI_relation,
        D_HI_arcsec_relation,
        atlas_source,
        unresolved,
        flux,
        w20,
        MHI_incl_scale,
        MHI_incl_v_scale,
        ska_PA,
    ) = make_cube(
        config,
        i,
        10 ** cat_gal["MHI"],
        10 ** cat_gal["Mh"],
        cat_gal["z"],
        cat_gal["i"],
        cat_gal["PA"],
        float(pixel_scale_str),
        dnu,
        psf_maj_arcsec,
        cat_gal["RadioClass"],
        cat_gal["OptClass"],
    )
    # unresolveds+=unresolved

    print("subcube shape and size:", sub_cube.shape, sub_cube.nbytes)
    sub_cube_shape = sub_cube.shape

    # test indexing:
    # sub_cube, D_HI,  D_HI_arcsec, D_HI_relation,  D_HI_arcsec_relation,atlas_source, unresolved, flux, w20, MHI_incl_scale,MHI_incl_v_scale \
    #                 =np.ones((1,1,1)), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1

    # sub_cube_shape = sub_cube.shape

    # have set cube to correct resolution in make_cube by reading pixel_scale from config

    l_bounds = np.array([v, y, x]) - np.array(
        [(sub_cube.shape[0] / 2), (sub_cube.shape[1] / 2), (sub_cube.shape[2] / 2)]
    )
    u_bounds = np.array([v, y, x]) + np.array(
        [
            sub_cube_shape[0] - (sub_cube.shape[0] / 2),
            sub_cube_shape[1] - (sub_cube.shape[1] / 2),
            sub_cube_shape[2] - (sub_cube.shape[2] / 2),
        ]
    )

    logging.info("Lower bounds, upper bounds: %s, %s", l_bounds, u_bounds)

    l_bounds = np.floor(l_bounds).astype(np.int)
    u_bounds = np.floor(u_bounds).astype(np.int)

    logging.info("Lower bounds, upper bounds, int: %s, %s", l_bounds, u_bounds)
    logging.info("Subcube shape: %s", sub_cube.shape)

    #### add it to the large cube:

    img3 = sub_cube
    blc0 = l_bounds[0]
    blc1 = l_bounds[1]
    blc2 = l_bounds[2]
    trc0 = u_bounds[0] - 1
    trc1 = u_bounds[1] - 1
    trc2 = u_bounds[2] - 1

    # the upper bounds are all -1 the true values, since the large cube is added using Fortran indexing

    # check whether subcube touches the edge of main cube
    # don't add if it does, in order to avoid wrapping problems
    trcs = np.array([trc0, trc1, trc2])

    blcs = np.array([blc0, blc1, blc2])
    # the top bounds are all -1 the true values, since the large cube is added using Fortran indexing

    top_excess = arr_dims - (
        trcs + 1
    )  # pixels from the image to top coordinates of field
    bottom_excess = blcs
    excess = np.hstack((bottom_excess, top_excess))

    overlap = False
    for coord in excess:
        if coord < 0:

            overlap = True
            logging.info("Subcube is overlapping the edge: not adding to full cube")
            logging.info("Adding no property values to table")
            # cat.remove_row(i)  this messes up table; deal with it by deleting at end
            break
    if overlap:

        print("img3 shape:", img3.shape)
        # plt.subplot(121)
        # plt.imshow(img3[0,:,:],origin = 'lower')
        # print ('overlapping')
        start_list = np.copy(bottom_excess)
        end_list = np.copy(top_excess)
        np.putmask(start_list, bottom_excess < 0, (-bottom_excess))
        np.putmask(start_list, bottom_excess >= 0, 0)
        start0, start1, start2 = start_list

        np.putmask(end_list, top_excess >= 0, img3.shape)
        end0, end1, end2 = end_list
        print(start0, start1, start2, end0, end1, end2)
        img3 = img3[start0:end0, start1:end1, start2:end2]

        # plt.subplot(122)
        # plt.imshow(img3[0,:,:],origin ='lower')
        # plt.show()
        # plt.clf()
        print("old coords:", blc0, blc1, blc2, trc0, trc1, trc2)

        np.putmask(blcs, bottom_excess < 0, 0)
        np.putmask(trcs, top_excess < 0, arr_dims - 1)
        print(arr_dims)
        blc0, blc1, blc2 = blcs
        trc0, trc1, trc2 = trcs
        print("new coords:", blc0, blc1, blc2, trc0, trc1, trc2)

    logging.info("BLC, TRC: %f %f %f, %f %f %f ", blc0, blc1, blc2, trc0, trc1, trc2)

    # open cube file for writing
    fitsf = FITS(all_gals_fname, "rw")
    # retrieve any signal that is already in source location
    region = fitsf[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
    if np.sum(region) > 0:
        logging.info("Flux already exists in region; adding new flux")
    img3 += region
    print("subcube shape and size:", img3.shape, img3.nbytes)
    fitsf[0].write(img3, blc0, blc1, blc2, trc0, trc1, trc2)

    t_source = time.time() - tstart

    # append catalogue with any new properties
    cat["HI_size_kpc"][i] = D_HI
    cat["HI_size"][i] = D_HI_arcsec
    cat["Atlas_source"][i] = atlas_source

    cat["line_flux_integral"][i] = flux
    cat["w20"][i] = w20
    cat["MHI_incl_v_scale"][i] = MHI_incl_v_scale
    cat["MHI_incl_scale"][i] = MHI_incl_scale
    print("pa fix:")
    print(cat["PA"][i])
    cat["PA"][i] = ska_PA
    print(cat["PA"][i])
    fitsf.close()
    logging.info("")

    return


def runSkyModel(config):
    """Simulate a sky model from a T-RECS catalogue.

    Parameters
    ----------
    config : configparser
        ConfigParser configuration containing necessary sections.

    """

    # logger = multiprocessing.get_logger()
    # logger.setLevel(logging.INFO)
    # logger.FileHandler('test.log')

    """
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    DEFAULT_LOGGING_FORMAT = '[%(levelname)s/%(processName)s] %(message)s'
    formatter = logging.Formatter(DEFAULT_LOGGING_FORMAT)
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    """

    # Set up logging
    logfilename = "logs/%s.log" % config.get("field", "fits_prefix")
    os.system("rm %s" % logfilename)
    log = logfilename
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
    )
    logging.info("Beginning simulation")

    n_cores = int(config.getfloat("pipeline", "n_cores"))
    logging.info("Running with %d cores", n_cores)

    doplot = config.getboolean("pipeline", "doplot")

    # set up cosmology

    H = config.getfloat("cosmology", "H")
    M = config.getfloat("cosmology", "M")
    L = config.getfloat("cosmology", "L")
    c = config.getfloat("cosmology", "c")
    G = config.getfloat("cosmology", "G")

    cosmo = LambdaCDM(H0=H, Om0=M, Ode0=L)

    data_path_large_files = (
        config.get("pipeline", "data_path_large_files")
        + config.get("pipeline", "project_name")
        + "/"
    )
    data_path = (
        config.get("pipeline", "base_dir")
        + config.get("pipeline", "data_path")
        + config.get("pipeline", "project_name")
        + "/"
    )
    # set image properties
    psf_maj_arcsec = config.getfloat("skymodel", "simple_psf_maj")
    psf_maj = psf_maj_arcsec * galsim.arcsec

    psf_min = config.getfloat("skymodel", "simple_psf_min") * galsim.arcsec
    psf_pa = config.getfloat("skymodel", "simple_psf_pa") * galsim.degrees
    pixel_scale = config.getfloat("skymodel", "pixel_scale")
    pixel_scale_str = str(pixel_scale).split()[0]
    fov = config.getfloat("skymodel", "field_of_view")
    logging.info("FoV from ini file, arcmin: %f", fov)

    # set spectral properties
    HI_line = config.getfloat("observation", "rest_freq")
    base_freq = config.getfloat("observation", "lowest_frequency")
    top_freq = config.getfloat("observation", "highest_frequency")
    bw = top_freq - base_freq
    fov, image_size = tools.get_image_size(fov, pixel_scale)
    logging.info("Unpadded image_size, power of two, pixels: %f", image_size)
    logging.info("FoV reduced to a power of two, arcmin: %f", fov)
    logging.info("Final base_freq, Hz: %f", base_freq)
    logging.info("Final top_freq, Hz: %f", top_freq)

    dnu = config.getfloat("observation", "channel_width")
    crpix3, n_chan = tools.get_spectral_size(bw, dnu)

    z_min = HI_line / top_freq - 1
    z_max = HI_line / base_freq - 1
    logging.info("Final cropped min and max redshifts: %f %f" % (z_min, z_max))

    # using optical definition throughout
    # need velocity values since the wcs is velocity-based (still linear in freq)
    velo_base = config.getfloat("cosmology", "c") * z_min
    velo_max = config.getfloat("cosmology", "c") * z_max
    logging.info(
        "Final cropped min and max velocities, m/s: %f %f" % (velo_base, velo_max)
    )
    arr_dims = np.array([n_chan, image_size, image_size]).astype(np.int)

    logging.info("Final array dimensions: %d %d %d" % (n_chan, image_size, image_size))
    logging.info(
        "Final array size, elements: %.3e" % (n_chan * image_size * image_size)
    )
    logging.info(
        "Final array size, bytes: %.3e" % (n_chan * image_size * image_size * 4)
    )

    # pad now, crop later
    # pad the three dimensions, and also cut the sources slightly wider than the final size but less than padding,
    # to retain those which overlap edges
    (
        pad_image_size,
        pad_fov,
        pad_base_freq,
        pad_top_freq,
        cut_base_freq,
        cut_top_freq,
    ) = tools.add_padding(image_size, fov, base_freq, top_freq, pixel_scale, dnu)

    logging.info("Padding values:")
    logging.info("pad_image_size, pixels: %f", pad_image_size)
    logging.info("pad_fov, arcmin: %f", pad_fov)
    logging.info("pad_base_freq, Hz: %f", pad_base_freq)
    logging.info("pad_top_freq, Hz: %f", pad_top_freq)
    logging.info("cut_base_freq, Hz: %f", cut_base_freq)
    logging.info("cut_top_freq, Hz: %f", cut_top_freq)

    cut_z_min = HI_line / cut_top_freq - 1
    cut_z_max = HI_line / cut_base_freq - 1
    logging.info("Source cut min and max redshifts: %f %f" % (cut_z_min, cut_z_max))
    # using optical definition throughout
    # need velocity values since the wcs is velocity-based (still linear in freq)
    cut_velo_base = config.getfloat("cosmology", "c") * cut_z_min
    cut_velo_max = config.getfloat("cosmology", "c") * cut_z_max
    logging.info(
        "Source cut min and max velocities, m/s: %f %f" % (cut_velo_base, cut_velo_max)
    )

    pad_bw = pad_top_freq - pad_base_freq
    pad_crpix3, pad_n_chan = tools.get_spectral_size(pad_bw, dnu)

    # set sky coordinates
    ra_field_gs = config.getfloat("field", "field_ra")
    # convert to range +/- 180 to enable cutoffs later
    if ra_field_gs > 180.0:
        ra_field_gs -= 360.0

    # dec_field = config.get('field', 'field_dec') ini file now uses deg format
    # dec_field_gs = galsim.Angle.from_dms(dec_field)  /galsim.degrees
    dec_field_gs = config.getfloat("field", "field_dec")

    # get wcs for fits header
    w_spectral = setup_wcs(config, ndim=3, cosmology=cosmo)
    w_twod = setup_wcs(config, ndim=2)
    #  w_fourd = setup_wcs(config, ndim=4)

    header_spectral = w_spectral.to_header()
    header_twod = w_twod.to_header()
    #  header_fourd = w_fourd.to_header()
    header_spectral["BUNIT"] = "Jy/beam"

    bmaj = psf_maj / galsim.degrees
    bmin = psf_min / galsim.degrees
    bpa = psf_pa / galsim.radians

    header_spectral["BMAJ"] = bmaj
    header_spectral["BMIN"] = bmin
    header_spectral["BPA"] = bpa

    header_twod["BMAJ"] = bmaj
    header_twod["BMIN"] = bmin
    header_twod["BPA"] = bpa

    # initialse an empty cube
    all_gals_fname = (
        data_path_large_files + config.get("field", "fits_prefix") + ".fits"
    )
    all_gals_summed_fname = (
        data_path + config.get("field", "fits_prefix") + "_image.fits"
    )

    os.system("rm {0}".format(all_gals_fname))
    logging.info("Creating empty image file, {0} ...".format(all_gals_fname))

    # first use fitsio to make huge cube by exploting the expand_if_needed functionality of the image write function

    logging.info(
        "Padded array dimensions: %d %d %d"
        % (pad_n_chan, pad_image_size, pad_image_size)
    )
    logging.info(
        "Padded array size, elements: %.2e"
        % (pad_n_chan * pad_image_size * pad_image_size)
    )
    logging.info(
        "Padded array size, bytes: %.2e"
        % (pad_n_chan * pad_image_size * pad_image_size * 4)
    )

    test_array_size = 0
    if test_array_size:
        data = np.zeros((pad_n_chan, pad_image_size, pad_image_size)).astype(np.float32)
        print("nbytes, 32 bits precision: ", data.nbytes)

    header_dict = {}
    logging.info("Making header:")
    # convert astropy header to dict for fitsio
    for key in header_spectral:
        header_dict[key] = header_spectral[key]
        logging.info("%s = %s" % (key, header_dict[key]))

    img = np.zeros((10, 10, 10)).astype(np.float32)

    fitsf = FITS(all_gals_fname, "rw")
    fitsf.write(img, header=header_dict)

    # start pixel coords need to be supplied in reverse order when expanding - fiddly but works
    # fitsio code changed to reverse axes ordering when writing (going from C to fortran order)

    img2 = np.zeros((1, 1, 1))

    blc0 = pad_image_size - 1  # -1 since the large cube is added using Fortran indexing
    blc1 = pad_image_size - 1
    blc2 = pad_n_chan - 1
    print(
        blc2,
        blc1,
        blc0,
        blc2 + img2.shape[2] - 1,
        blc1 + img2.shape[1] - 1,
        blc0 + img2.shape[0] - 1,
    )
    fitsf[0].write(
        img2,
        blc2,
        blc1,
        blc0,
        blc2 + img2.shape[2] - 1,
        blc1 + img2.shape[1] - 1,
        blc0 + img2.shape[0] - 1,
    )

    fitsf.close()

    cr_x, cr_y = w_twod.wcs_world2pix(
        ra_field_gs,
        dec_field_gs,
        1,
    )

    logging.info("Check world2pix crpix1, crpix2: %f %f", cr_x, cr_y)

    if config.getboolean("skymodel", "doHI"):
        # Load the catalogue
        cat_file_name = (
            config.get("pipeline", "base_dir")
            + config.get("pipeline", "data_path")
            + config.get("field", "catalogue")
        )

        #  #AB:file where to record the actual galaxy positions
        #  outf=cat_file_name+'_pos_offsets'
        #  f= open(outf,"w+")

        logging.info("Loading catalogue from {0} ...".format(cat_file_name))
        cat = Table()

        cat_read = Table.read(cat_file_name)  # remove ascii

        source_prefix = "TRECS-"
        source_name_ra = np.asarray(cat_read["longitude"], dtype=str)
        source_name_dec = np.asarray(cat_read["latitude"], dtype=str)

        source_prefix_arr = strmultiply(
            source_prefix, np.ones_like(source_name_ra, dtype=int)
        )
        source_l_arr = strmultiply("l", np.ones_like(source_name_ra, dtype=int))
        source_b_arr = strmultiply("b", np.ones_like(source_name_ra, dtype=int))

        source_name_pos = stradd(
            source_prefix_arr,
            stradd(
                source_l_arr,
                stradd(source_name_ra, (stradd(source_b_arr, source_name_dec))),
            ),
        )
        cat["id"] = np.zeros(len(cat_read))

        cat["RA"] = cat_read["longitude"]  # deg

        cat["ra_offset"] = cat["RA"] - ra_field_gs  # deg
        cat["ra_offset"].unit = "deg"

        cat["Dec"] = cat_read["latitude"]  # deg

        cat["dec_offset"] = cat_read["latitude"] - dec_field_gs  # deg
        cat["dec_offset"].unit = "deg"
        dec_abs_radians = cat["Dec"] * galsim.degrees / galsim.radians

        cat["HI_size2"] = cat_read["HI size"]  # arcsec, not necessarily correct atm
        cat["HI_size2"].unit = "arcsec"

        cat["HI_flux"] = cat_read["HI flux"]
        cat["HI_flux"].unit = "Jy ?"

        cat["MHI"] = cat_read["MHI"]
        cat["MHI"].unit = "M_solar"  # ?

        cat["Mh"] = cat_read["Mh"]
        cat["Mh"].unit = "M_solar"  # ?

        cat["z"] = cat_read["redshift"]
        cat["z"].unit = "none"

        cat["i"] = cat_read["inclination"]

        cat["PA"] = cat_read["PA"]
        cat["PA"].unit = "deg"

        cat["RadioClass"] = cat_read["RadioClass"]
        cat["OptClass"] = cat_read["OptClass"]

        cat["opt_vel"] = config.getfloat("cosmology", "c") * (cat["z"])
        cat["opt_vel"].unit = "m/s"

        cat["central_freq"] = HI_line / (cat["z"] + 1)
        cat["central_freq"].unit = "Hz"

        # new cat attributes
        cat["Atlas_source"] = np.zeros(len(cat_read)).astype(np.str)

        cat["HI_size_kpc"] = np.zeros(len(cat_read))
        cat["HI_size_kpc"].unit = "kpc"

        cat["HI_size"] = np.zeros(len(cat_read))
        cat["HI_size"].unit = "arcsec"

        cat["HI_size_relation_kpc"] = np.zeros(len(cat_read))
        cat["HI_size_relation_kpc"].unit = "kpc"

        cat["HI_size_relation"] = np.zeros(len(cat_read))
        cat["HI_size_relation"].unit = "arcsec"

        cat["line_flux_integral"] = np.zeros(len(cat_read))
        cat[
            "line_flux_integral"
        ].unit = "Jy Hz"  # using freq version of MHI-line flux conversion (Duffy+12)
        cat["w20"] = np.zeros(len(cat_read))
        cat["w20"].unit = "km/s"

        cat["MHI_incl_scale"] = np.zeros(len(cat_read))
        cat["MHI_incl_v_scale"] = np.zeros(len(cat_read))

        logging.info("Cat length before fov cut: %d", len(cat))

        logging.info(
            "ra offset min max: %f %f",
            np.min(cat["ra_offset"]),
            np.max(cat["ra_offset"]),
        )
        logging.info(
            "dec offset min max: %f %f",
            np.min(cat["dec_offset"]),
            np.max(cat["dec_offset"]),
        )

        # fov cut, put cos(dec) factor into ra offset
        cosdec = np.cos(dec_field_gs * 2 * np.pi / 360)

        # need to think about removing cosdec since higher-up sources are not filling plane
        ra_offset_max = (1 / cosdec) * (
            (pad_fov / 60) / 2
        )  # this is now a fraction of the written image size
        dec_offset_max = (pad_fov / 60) / 2  # convert fov to degrees

        fov_cut = (abs(cat["ra_offset"]) < ra_offset_max) * (
            abs(cat["dec_offset"]) < dec_offset_max
        )

        cat = cat[fov_cut]

        logging.info(
            "ra offset min max: %f %f",
            np.min(cat["ra_offset"]),
            np.max(cat["ra_offset"]),
        )
        logging.info(
            "dec offset min max: %f %f",
            np.min(cat["dec_offset"]),
            np.max(cat["dec_offset"]),
        )

        # using optical definition throughout
        cut_z_min = cut_velo_base / config.getfloat("cosmology", "c")
        cut_z_max = cut_velo_max / config.getfloat("cosmology", "c")

        logging.info("source cut min and max redshifts: %f %f" % (z_min, z_max))
        logging.info(
            "source cut min and max velocities, m/s: %f %f"
            % (cut_velo_base, cut_velo_max)
        )
        vel_cut = ((cat["opt_vel"]) > cut_velo_base) * ((cat["opt_vel"]) < cut_velo_max)
        cat = cat[vel_cut]

        MHI_cut = cat["MHI"] > 0

        cat = cat[MHI_cut]

        #  cat = cat[cat['MHI']>10.5]

        logging.info("Cat length after fov and z cut: %d", len(cat))

        # number of sources, on grid if requested
        if config.getboolean("skymodel", "grid"):
            nobj = int(np.sqrt(config.getint("skymodel", "ngals"))) ** 2.0
            cat["ra_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
            cat["dec_offset"] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
        else:
            nobj = len(cat)

            if config.getint("skymodel", "ngals") > -1:
                nobj = config.getint("skymodel", "ngals")

                cat = cat[:nobj]

        ix_arr = np.ones(nobj)
        iy_arr = np.ones(nobj)
        iv_arr = np.ones(nobj)

        #    unresolveds = 0

        # plot some catalogue statistics

        if doplot:
            H, xedges, yedges = np.histogram2d(
                (cat["MHI"]), np.cos(cat["i"] / 180 * np.pi) * 10
            )  # , bins=(xedges, yedges))
            H = H.T  # Let each row list bins with common y range
            fig = plt.figure(figsize=(7, 7))
            ax = fig.add_subplot(111, title="imshow: square bins")
            plt.imshow(
                H,
                interpolation="nearest",
                origin="low",
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            )
            plt.ylabel("10*np.cos(inclination (degrees))")
            plt.xlabel("log10MHI (M_solar)")
            plt.show()

        for i, cat_gal in enumerate(cat):
            add_source(
                i,
                cat_gal,
                nobj,
                w_spectral,
                config,
                pixel_scale_str,
                dnu,
                psf_maj_arcsec,
                arr_dims,
                all_gals_fname,
                cat,
            )
        # if i == 100000:
        #     break
        """

        print ('going into loop')
        pool = multiprocessing.Pool(n_cores) 
        for i,cat_gal in enumerate(cat): 
        
    
        # Draw the galaxies onto the galsim image

            pool.apply_async(add_source, args = (i,cat_gal, nobj, w_spectral, config, pixel_scale_str, dnu, psf_maj_arcsec, arr_dims, all_gals_fname, cat ))#, callback = log_result)
        pool.close()
        pool.join()
        """

        atlas_sources = cat["Atlas_source"]
        filled_rows = np.argwhere(atlas_sources != "0.0")[:, 0]
        cat = cat[filled_rows]
        cat["id"] = np.arange(len(cat))

        # write out catalogues

        # restricted catalogue: for ska use
        truthcat_name = (
            data_path + config.get("field", "fits_prefix") + "_restricted_truthcat.fits"
        )
        logging.info(
            "Writing restricted truth catalogue to: {0} ...".format(truthcat_name)
        )
        cat.write(truthcat_name, format="fits", overwrite=True)

        # public catalogue: for participant use
        public_cat = Table()
        public_cat["id"] = cat["id"]
        public_cat["RA"] = cat["RA"]
        public_cat["RA"].unit = "deg"
        public_cat["Dec"] = cat["Dec"]
        public_cat["Dec"].unit = "deg"
        #    public_cat['HI_size_kpc'] = cat['HI_size_kpc']
        #    public_cat['HI_size_kpc'].unit = 'kpc'
        public_cat["HI_size"] = cat["HI_size"]
        public_cat["HI_size"].unit = "arcsec"
        public_cat["line_flux_integral"] = cat["line_flux_integral"]
        public_cat[
            "line_flux_integral"
        ].unit = "Jy Hz"  # using freq version of MHI-line flux conversion (Duffy+12)
        public_cat["central_freq"] = cat["central_freq"]
        public_cat["central_freq"].unit = "Hz"
        public_cat["PA"] = cat["PA"]
        public_cat["PA"].unit = "deg"
        public_cat["i"] = cat["i"]
        public_cat["i"].unit = "deg"
        #    public_cat['MHI'] = cat['MHI']
        #    public_cat['MHI'].unit = 'M_solar'
        #    public_cat['opt_vel'] =  cat['opt_vel']
        #    public_cat['opt_vel'].unit = 'm/s'
        #    public_cat['z'] = cat['z']
        #    public_cat['z'].unit = 'none'
        public_cat["w20"] = cat["w20"]
        public_cat["w20"].unit = "km/s"

        public_truthcat_name = (
            data_path + config.get("field", "fits_prefix") + "_public_truthcat.fits"
        )
        logging.info(
            "Writing public truth catalogue to: {0} ...".format(public_truthcat_name)
        )
        public_cat.write(public_truthcat_name, format="fits", overwrite=True)

        #   logging.info('Unresolveds: %d', unresolveds)
        logging.info("Writing summed cube to: {0} ...".format(all_gals_summed_fname))
        # save a 2d image using astropy instead of fitsio
        flux = astfits.getdata(all_gals_fname)
        summed_line_flux = np.sum(flux, axis=0)
        hdu = astfits.PrimaryHDU(summed_line_flux, header=header_twod)

        hdulist = astfits.HDUList([hdu])

        hdulist.writeto(all_gals_summed_fname, overwrite=True)
        print(all_gals_summed_fname)
        # logging.info('Compressing cube')
        # need to escape ( in path os.system('zip %s %s'%(data_path+config.get('field', 'fits_prefix')+'.fits'.strip('fits')+'zip', all_gals_fname))

        # os.system('zip %s %s'%(all_gals_fname.strip('fits')+'zip', all_gals_fname))

        #  logging.info('Producing diagnostic plots')

        tend = time.time()
        logging.info("...done in {0} seconds.".format(tend - tstart))
        print(tend - tstart)

        print("preparing check by doing simple sum of abs values")
        print("sum of absolute 3d array values is %f" % np.abs(flux).sum())

        print("preparing checksum")
        os.system(
            "md5sum %s > %s" % (all_gals_fname, all_gals_fname.strip("fits") + "md5")
        )

        exit()

    all_gals_fname = data_path + config.get("field", "fitsname")
    print("Writing image data to {0} ...".format(all_gals_fname))

    f.close()  # AB: close the file with true positions

    mom0 = np.sum(data, axis=0)

    mom0 = astfits.PrimaryHDU(mom0)
    mom0.writeto(data_path + "mom0.fits", overwrite=True)
    # Extract the numpy array from the galsim image
    image_data = data

    if config.getboolean("primarybeam", "dopb"):

        nstokes = config.getint("primarybeam", "nstokes")
        nfreq = config.getint("primarybeam", "nfreq")
        bw = config.getfloat("observation", "total_bandwidth")
        base_freq = config.getfloat("observation", "lowest_frequency")
        freq_width = bw / nfreq

        image_cube = np.empty((nstokes, nfreq) + image_data.shape)

        stokes_list = ["I", "Q", "U"][:nstokes]
        freq_list = base_freq + np.arange(nfreq) * freq_width

        for i_stokes, stokes in enumerate(stokes_list):
            for i_freq, freq in enumerate(freq_list):
                image_cube[i_stokes, i_freq] = image_data * primary_beam(config, freq)

        hdu = astfits.PrimaryHDU(image_cube, header=header_fourd)

    else:
        if config.getboolean("skymodel", "dosimple_psf"):
            header_fourd["BMAJ"] = psf_maj / galsim.degrees
            header_fourd["BMIN"] = psf_min / galsim.degrees
            header_fourd["BPA"] = psf_pa / galsim.radians5

            noise_sigma = image_data.max() / config.getfloat(
                "skymodel", "simple_psf_snr"
            )
            noise_data = np.random.normal(
                loc=0.0, scale=noise_sigma, size=image_data.shape
            )
            image_data = image_data + noise_data

            hdu_noise = astfits.PrimaryHDU(
                np.expand_dims(np.expand_dims(noise_data, axis=0), axis=0),
                header=header_spectral,
            )
            hdulist_noise = astfits.HDUList([hdu_noise])
            all_gals_fname_noise = all_gals_fname.split(".")[0] + ".noise.fits"
            hdulist_noise.writeto(all_gals_fname_noise, clobber=True)

            psf_image = galsim.ImageF(
                image_size, image_size, scale=pixel_scale / galsim.arcsec
            )
            psf.draw(psf_image)
            psf_data = psf_image.array
            hdu_psf = astfits.PrimaryHDU(
                np.expand_dims(np.expand_dims(psf_data, axis=0), axis=0),
                header=header_spectral,
            )
            hdulist_psf = astfits.HDUList([hdu_psf])
            all_gals_fname_psf = all_gals_fname.split(".")[0] + ".psf.fits"
            hdulist_psf.writeto(all_gals_fname_psf, clobber=True)

        hdu = astfits.PrimaryHDU(
            np.expand_dims(np.expand_dims(image_data, axis=0), axis=0),
            header=header_spectral,
        )

    hdulist = astfits.HDUList([hdu])
    hdulist.writeto(all_gals_fname, overwrite=True)

    print("...done.")

    if config.getboolean("skymodel", "im3cat"):
        np.savetxt(
            data_path + config.get("pipeline", "project_name") + "_im3cat.txt",
            np.column_stack([np.arange(nobj), ix_arr, iy_arr]),
        )


if __name__ == "__main__":

    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])

    runSkyModel(config)
