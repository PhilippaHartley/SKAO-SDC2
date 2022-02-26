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
from skymodel.continuum_morphology import make_img
from skymodel.HI_morphology import make_cube
from skymodel.skymodel_tools import setup_wcs

# from primarybeam.primarybeam import *


tstart = time.time()

# mother_seed=5820743 #seed for random number generation - fullcube
mother_seed = 6879432  # seed for random number generation - smallcube1
# mother_seed=7984532 #seed for random number generation -smallcube2


arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value

"""
def initialise_file(filename, ):
    
def make_header(header_type,):

def retrieve_flux():
"""


def log_result_continuum(result):
    return


def add_source_continuum():

    return


def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    #  i = result[0]
    #  D_HI = result[1]
    #  cat["HI_size_kpc"][i] = D_HI

    global cat

    # cat["HI_size"][i] = result[2]
    print("************************")
    print(result)
    (
        i,
        D_HI,
        D_HI_arcsec,
        atlas_source,
        flux,
        w20,
        MHI_incl_v_scale,
        MHI_incl_scale,
        ska_PA,
    ) = result

    print("*****************", D_HI)
    cat["HI_size_kpc"][i] = D_HI
    cat["HI_size"][i] = D_HI_arcsec
    cat["Atlas_source"][i] = atlas_source
    cat["line_flux_integral"][i] = flux
    cat["w20"][i] = w20
    cat["MHI_incl_v_scale"][i] = MHI_incl_v_scale
    cat["MHI_incl_scale"][i] = MHI_incl_scale
    cat["PA"][i] = ska_PA


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
    root,
):
    print("making source")
    mainlog = logging.getLogger("main%d" % i)
    h = logging.FileHandler("log%d.log" % i)
    mainlog.addHandler(h)
    logging.root.setLevel(logging.DEBUG)
    mainlog.info("result%s" % i)

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
    # clip if it does
    trcs = np.array([trc0, trc1, trc2])
    blcs = np.array([blc0, blc1, blc2])

    top_excess = arr_dims - (
        trcs + 1
    )  # pixels from the image to top coordinates of field
    bottom_excess = blcs
    excess = np.hstack((bottom_excess, top_excess))

    overlap = False
    for coord in excess:
        if coord < 0:
            overlap = True
            logging.info("Subcube is overlapping the edge: clipping")

            break
    if overlap:
        start_list = np.copy(bottom_excess)
        end_list = np.copy(top_excess)
        np.putmask(start_list, bottom_excess < 0, (-bottom_excess))
        np.putmask(start_list, bottom_excess >= 0, 0)
        start0, start1, start2 = start_list
        np.putmask(end_list, top_excess >= 0, img3.shape)
        end0, end1, end2 = end_list
        img3 = img3[start0:end0, start1:end1, start2:end2]
        np.putmask(blcs, bottom_excess < 0, 0)
        np.putmask(trcs, top_excess < 0, arr_dims - 1)
        blc0, blc1, blc2 = blcs
        trc0, trc1, trc2 = trcs

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

    fitsf.close()
    logging.info("")
    print(i)
    mainlog.info("done_make_cube%s" % i)
    with open("log%d.log" % i, "r") as f:
        a = f.readlines()
        print(a)

    with open("all_log.txt", "a") as f:
        f.writelines(a)
    os.system("rm log%d.log" % i)

    return (
        i,
        D_HI,
        D_HI_arcsec,
        atlas_source,
        flux,
        w20,
        MHI_incl_v_scale,
        MHI_incl_scale,
        ska_PA,
    )


def runSkyModel(config):
    """Simulate a sky model from a T-RECS catalogue.

    Parameters
    ----------
    config : configparser
        ConfigParser configuration containing necessary sections.

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

    # set sky coordinates
    ra_field_gs = config.getfloat("field", "field_ra")
    # convert to range +/- 180 to enable cutoffs later
    if ra_field_gs > 180.0:
        ra_field_gs -= 360.0

    dec_field_gs = config.getfloat("field", "field_dec")

    if config.getboolean("skymodel", "doHI"):

        # set spectral properties
        HI_line = config.getfloat("observation", "rest_freq")
        base_freq = config.getfloat("observation", "lowest_frequency")
        top_freq = config.getfloat("observation", "highest_frequency")
        bw = top_freq - base_freq
        fov, image_size = tools.get_image_size(fov, pixel_scale)
        logging.info("Image_size, power of two, pixels: %f", image_size)
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

        logging.info(
            "Final array dimensions: %d %d %d" % (n_chan, image_size, image_size)
        )
        logging.info(
            "Final array size, elements: %.3e" % (n_chan * image_size * image_size)
        )
        logging.info(
            "Final array size, bytes: %.3e" % (n_chan * image_size * image_size * 4)
        )

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

        test_array_size = 0
        if test_array_size:
            data = np.zeros((n_chan, image_size, image_size)).astype(np.float32)
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

        blc0 = image_size - 1  # -1 since the large cube is added using Fortran indexing
        blc1 = image_size - 1
        blc2 = n_chan - 1
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
        global cat
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
            (fov / 60) / 2
        )  # this is now a fraction of the written image size
        dec_offset_max = (fov / 60) / 2  # convert fov to degrees

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

        vel_cut = ((cat["opt_vel"]) > velo_base) * ((cat["opt_vel"]) < velo_max)
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
        pool = 1
        if not pool:
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
        if pool:
            print("going into loop")
            multiprocessing.set_start_method("fork")
            pool = multiprocessing.Pool(n_cores)
            for i, cat_gal in enumerate(cat):
                mainlog = logging.getLogger("main%d" % i)
                h = logging.FileHandler("log%d.log" % i)
                mainlog.addHandler(h)
                logging.root.setLevel(logging.DEBUG)
                mainlog.info("test%s" % i)

                pool.apply_async(
                    add_source,
                    args=(
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
                        mainlog,
                    ),
                    callback=log_result,
                )

            pool.close()
            pool.join()
            print(cat)
            print(cat["Atlas_source"][i])

        atlas_sources = cat["Atlas_source"]

        filled_rows = np.argwhere(atlas_sources != "0.0")[
            :, 0
        ]  # this should no longer be needed as clipping
        len1 = len(cat)

        cat = cat[filled_rows]
        len2 = len(cat)
        print(len1, len2)

        if len1 != len2:
            print("some sources were not added")
            exit()

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
    if config.getboolean("skymodel", "doagn"):
        # set spectral properties
        HI_line = config.getfloat("observation", "rest_freq")
        base_freq = config.getfloat("observation", "lowest_frequency")  # *1.e6 #Hz
        base_freqname = config.get("observation", "lowest_frequency")
        top_freq = config.getfloat("observation", "highest_frequency")  # *1.e6 #Hz
        top_freqname = config.get("observation", "highest_frequency")

        bw = top_freq - base_freq
        fov, image_size = tools.get_image_size(fov, pixel_scale)
        print(fov, pixel_scale, image_size)

        logging.info("Image_size, power of two, pixels: %f", image_size)
        logging.info("FoV reduced to a power of two, arcmin: %f", fov)
        logging.info("Final base_freq, Hz: %f", base_freq)
        logging.info("Final top_freq, Hz: %f", top_freq)

        # Continuum version:
        # definition of bandwidth and channels for continuum is different.
        # Explicit definition of low and up frequency via input file
        # a cube is generated with continuum channels

        dnu = config.getfloat("observation", "channel_width")
        nfreqs = int((top_freq - base_freq) / dnu) + 1
        print("nfreqs", nfreqs)
        freqs = np.zeros(nfreqs).astype(np.float32)
        print(freqs.shape)

        freq = base_freq
        for ff in range(nfreqs):
            freqs[ff] = freq
            print(freq)
            freq = freq + dnu

        crpix3 = 1
        n_chan = nfreqs

        arr_dims = np.array([n_chan, image_size, image_size]).astype(np.int)

        logging.info(
            "Final array dimensions: %d %d %d" % (n_chan, image_size, image_size)
        )
        logging.info(
            "Final array size, elements: %.3e" % (n_chan * image_size * image_size)
        )
        logging.info(
            "Final array size, bytes: %.3e" % (n_chan * image_size * image_size * 4)
        )

        # ANNA: here I am calling the appropriate wcs for continuum.
        # Based on header_4D plus extra keywords
        # TODO: Carta does not like something about this header. Probably something not standard.
        # get wcs for fits header
        w_spectral = setup_wcs(config, ndim=3, cosmology=cosmo)
        w_twod = setup_wcs(config, ndim=2, cosmology=cosmo)
        w_fourd = setup_wcs(config, ndim=4)

        header_spectral = w_spectral.to_header()
        header_twod = w_twod.to_header()
        header_fourd = w_fourd.to_header()
        header_fourd["BUNIT"] = "JY/PIXEL"

        header_spectral["BUNIT"] = "JY/BEAM"

        bmaj = psf_maj / galsim.degrees
        bmin = psf_min / galsim.degrees
        bpa = psf_pa / galsim.radians

        header_spectral["BMAJ"] = bmaj
        header_spectral["BMIN"] = bmin
        header_spectral["BPA"] = bpa

        header_fourd["BMAJ"] = bmaj
        header_fourd["BMIN"] = bmin
        header_fourd["BPA"] = bpa

        # initialse an empty cube
        all_gals_fname = (
            data_path_large_files + config.get("field", "fits_prefix") + ".fits"
        )

        # the summed cube is not meaningful for continuum
        # all_gals_summed_fname = data_path+config.get('field', 'fits_prefix')+'_image.fits'

        os.system("rm {0}".format(all_gals_fname))
        os.system("rm {0}".format(all_gals_fname + "_z.fits"))
        os.system("rm {0}".format(all_gals_fname + "_maxflux.fits"))
        logging.info("Creating empty image file, {0} ...".format(all_gals_fname))

        # first use fitsio to make huge cube by exploting the expand_if_needed functionality of the image write function

        logging.info(
            "Padded array dimensions: %d %d %d" % (n_chan, image_size, image_size)
        )
        logging.info(
            "Padded array size, elements: %.2e" % (n_chan * image_size * image_size)
        )
        logging.info(
            "Padded array size, bytes: %.2e" % (n_chan * image_size * image_size * 4)
        )

        test_array_size = 0
        if test_array_size:
            data = np.zeros((n_chan, image_size, image_size)).astype(np.float32)
            print("nbytes, 32 bits precision: ", data.nbytes)

        header_dict = {}
        logging.info("Making header:")

        # convert astropy header to dict for fitsio
        for key in header_fourd:

            header_dict[key] = header_fourd[key]
            logging.info("%s = %s" % (key, header_dict[key]))

        ###ANNA: qui metterre if per scegliere se spectral o continuum
        #    for key in header_spectral:

        #        header_dict[key] = header_spectral[key]
        #        logging.info('%s = %s'%(key,header_dict[key]))

        # img=np.zeros((10,10,10)).astype(np.float32)

        # fitsf = FITS(all_gals_fname,'rw')
        # fitsf.write(img, header = header_dict)

        # exit()
        # start pixel coords need to be supplied in reverse order when expanding - fiddly but works
        # fitsio code changed to reverse axes ordering when writing (going from C to fortran order)

        blc0 = image_size - 1  # -1 since the large cube is added using Fortran indexing
        blc1 = image_size - 1
        blc2 = n_chan - 1

        print("array dimensions", blc0, blc1, blc2)

        n_x = blc0  # naxis1 needed to decide whether to process source
        n_y = blc1  # naxis2 needed to decide whether to process source
        print(n_x, n_y)
        # ANNA: this feature seemed to be not present in Philippa's code. If that's the case, an image is prepared for every source in the catalogue and discarded at the end if out of FoV. For efficiency, if can be a good idea to check that beforehand for the centre of the source and process only sources whose centre is within FoV.

        # here we create 4 maps:
        # 1) continuum cube with all sources for the chosen continuum channels
        # 2) redshift map containing the redshift of the brightest source oon the LoS - needed for HI absorption
        # 3) map of the brightest fluxes along the LoS - This is needed for obtaining 3 and could just be a trowaway dummy map. For the moment I am keeping for debugging purposes

        img2 = np.zeros((blc2 + 1, blc0 + 1, blc1 + 1)).astype(
            np.float32
        )  # empty array for the continuum cube
        img2_1D = np.zeros((1, blc0 + 1, blc1 + 1)).astype(
            np.float32
        )  # empy array for 1D maps

        fitsf = FITS(all_gals_fname, "rw")
        # files storing the brightest object and the redshift of the brightest object
        fitsf_f = FITS(all_gals_fname + "_maxflux.fits", "rw")
        fitsf_z = FITS(all_gals_fname + "_z.fits", "rw")

        # initialise files
        fitsf.write(img2, header=header_dict)
        fitsf_f.write(img2_1D, header=header_dict)
        fitsf_z.write(img2_1D, header=header_dict)

        fitsf.close()
        fitsf_f.close()
        fitsf_z.close()

        fitsf = FITS(all_gals_fname, "rw")

        cr_x, cr_y = w_twod.wcs_world2pix(
            ra_field_gs,
            dec_field_gs,
            1,
        )

        logging.info("Check world2pix crpix1, crpix2: %f %f", cr_x, cr_y)

        if config.getboolean("skymodel", "doagn"):
            # reading the catalogue. The keywork could become 'docontinuum' an SFGs and AGNs done together here

            # Load the catalogue
            # cat_file_name = config.get('pipeline', 'base_dir')+ config.get('pipeline', 'data_path')+config.get('field', 'catalogue')

            HI_cross = False  # assume that the catalogue is con cross-matched with HI
            cat_file_name = config.get("field", "catalogue")

            outf = cat_file_name + "_pos_offsets"
            f = open(outf, "w+")

            logging.info("Loading catalogue from {0} ...".format(cat_file_name))
            cat = Table()

            cat_read = Table.read(cat_file_name)  # remove ascii
            keywords = cat_read.colnames

            if "M_HI" in keywords:
                HI_cross = True

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
            cat["Source_id"] = source_name_pos

            cat["RA"] = cat_read["longitude"]  # deg

            cat["ra_offset"] = cat["RA"] - ra_field_gs  # deg
            cat["ra_offset"].unit = "deg"

            cat["DEC"] = cat_read["latitude"]  # deg

            cat["dec_offset"] = cat_read["latitude"] - dec_field_gs  # deg
            cat["dec_offset"].unit = "deg"

            dec_abs_radians = cat["DEC"] * galsim.degrees / galsim.radians

            z = cat_read["redshift"]
            if HI_cross == True:
                z_1 = cat_read["redshift_1"]
                z[z == -100] = z_1[z == -100]
            cat["z"] = z
            cat["z"].unit = "none"

            # each source is approximated as a power law within the cube. A spectral index is computed between the lowest and highest specified frequencies.
            # This approximation is OK for channels within the same band.
            # For frequencies belonging do different bands, perform multiple runs of the code.

            cat["Total_flux"] = cat_read["I" + base_freqname] * 1.0e-3  # Jy
            cat["Total_flux"].unit = "Jy"

            cat["flux2"] = cat_read["I" + top_freqname] * 1.0e-3  # Jy
            cat["flux2"].unit = "Jy"

            cat["spectral_index"] = np.log10(
                cat["Total_flux"] / cat["flux2"]
            ) / np.log10(base_freq / top_freq)

            # read the relevant quantities to implement flux cuts
            if config.getboolean("skymodel", "highfluxcut") == True:
                highflux = config.getfloat("skymodel", "highfluxcut_value")
                flux_sel_freq = config.get("skymodel", "fluxcut_frequency")
                print(highflux, flux_sel_freq)
                cat["flux_selection"] = cat_read["I" + flux_sel_freq] * 1.0e-3  # Jy

            if config.getboolean("skymodel", "lowfluxcut") == True:
                lowflux = config.getfloat("skymodel", "lowfluxcut_value")
                flux_sel_freq = config.get("skymodel", "fluxcut_frequency")
                cat["flux_selection"] = cat_read["I" + flux_sel_freq] * 1.0e-3  # Jy

            maj = cat_read["size"]  # arcsec
            cat["Maj"] = maj

            cat["Maj"].unit = "arcsec"

            q = cat_read["axis ratio"]
            if HI_cross == True:
                q1 = cat_read["axis ratio_1"]
                q[q == -100] = q1[q == -100]
            cat["Min"] = maj * q

            cat["Min"].unit = "arcsec"

            # ANNA: check if those are still needed
            scale_radius_to_hlr = 1.67834699
            cat["Maj_halflight"] = cat_read["size"] * scale_radius_to_hlr
            cat["Maj_halflight"].unit = "arcsec"

            cat["Min_halflight"] = cat_read["size"] * scale_radius_to_hlr
            cat["Min_halflight"].unit = "arcsec"

            cat["Peak_flux"] = cat["Total_flux"]  # / (2.*cat['Maj']*arcsectorad)
            cat["Peak_flux"].unit = "Jy"
            # ANNA: end check if those are still needed

            cat["Rs"] = cat_read["Rs"]

            rdcl = cat_read["RadioClass"]  # to be used in source selection
            cat["RadioClass"] = rdcl

            pa = cat_read[
                "PA"
            ]  # this is the HI PA. rotate of 90 degs for AGN counterparts

            if HI_cross == True:
                cat["MHI"] = cat_read[
                    "MHI"
                ]  # this needed to select only Hi counterparts - for test purposes
                # PA in continuum is the HI PA rotated by 90 degs
                pa_copy = pa + 90.0
                pa_copy[pa_copy > 359.0] = pa_copy[pa_copy > 359.0] - 360.0
                pa_copy[pa == -100] = -100.0

                pa[rdcl > 3] = pa_copy[rdcl > 3]  # AGN PA 90degs from HI PA.
                pa_1 = cat_read["PA_1"]

                np.random.seed(mother_seed + 1)
                pa_2 = np.random.uniform(
                    low=0, high=359.0, size=len(cat)
                )  # PA not defined for AGN, here generate random

                pa[pa == -100] = pa_1[pa == -100]
                pa[pa == -100] = pa_2[pa == -100]

                # free
                pa_1 = 0
                pa_2 = 0
                pa_copy = 0

            cat["PA"] = pa
            cat["PA"].unit = "deg"

            # selects only continuum, AGNs
            #  cat = cat[(cat['RadioClass']>3)*(cat['RadioClass']!=-100)*(cat['Maj']<3600.)] #exclude too big - memory problem and not realistic

            # cat = cat[(cat['RadioClass']<4)*(Cat['RadioClass']!=-100)] #sfg
            # cat = cat[(cat['RadioClass']<4)*(cat['RadioClass']!=-100)*(cat['Maj']>10.)] #resolved sfg

            # select continuum sources
            cat = cat[(cat["RadioClass"] != -100) * (cat["Maj"] < 200.0)]
            # exclude too big - memory problem and not realistic] #select only continuum
            # print(len(cat))
            # cat = cat[(cat['MHI']!=-100.)] # only to get continuum only #change change
            # cat = cat[(cat['RadioClass']<4)] # only to get SFG only #change change
            # print(len(cat))
            #  exit()

            if config.getboolean("skymodel", "highfluxcut") == True:
                print("applying high flux cut")
                len_old = len(cat)

                cat = cat[(cat["flux_selection"] < highflux)]

                print("number of sources excluded")
                print(len_old - len(cat))

            if config.getboolean("skymodel", "lowfluxcut") == True:
                print("applying low flux cut")
                len_old = len(cat)

                cat = cat[(cat["flux_selection"] > lowflux)]

                print("number of sources excluded")
                print(len_old - len(cat))

            # define additional source attributes not contained in the TRECS cat
            cat["Atlas_source"] = np.zeros(len(cat)).astype(np.str)
            np.random.seed(mother_seed + 100)

            corefrac = np.random.normal(
                loc=0.75, scale=0.1, size=len(cat)
            )  # initialise core fraction. Steep-spectrum AGN dont use it as it is determined by the postage stamp.
            cat["corefrac"] = corefrac
            np.random.seed(mother_seed + 1000)

            # this random number is used later to associate sources to postage stamps
            ranid = np.random.uniform(low=0, high=1, size=len(cat))

            ranid[cat["Rs"] <= 0.5] = ranid[cat["Rs"] <= 0.5] - 10
            ranid[cat["Rs"] > 0.5] = ranid[cat["Rs"] > 0.5] + 10
            cat["ranid"] = ranid

            if config.get("skymodel", "sizescale") == "constant":
                cat["Maj"] = np.ones_like(cat["Maj"]) * config.getfloat(
                    "skymodel", "sizescale_constant_value"
                )

                #      scale_radius_to_hlr = galsim.Exponential(scale_radius=1., flux=1.).calculateHLR()
                scale_radius_to_hlr = 1.67834699
                cat["Maj_halflight"] = cat["Maj"] * scale_radius_to_hlr
                cat["Maj_halflight"].unit = "arcsec"

                cat["Min"] = cat["Maj"] * cat["q"]
                cat["Min"].unit = "arcsec"

                cat["Min_halflight"] = cat["Maj_halflight"] * cat["q"]
                cat["Min_halflight"].unit = "arcsec"

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

            # flux range
            if config.get("skymodel", "fluxscale") == "constant":
                cat["Total_flux"] = np.ones_like(cat["Total_flux"]) * config.getfloat(
                    "skymodel", "fluxscale_constant_value"
                )
                cat["Peak_flux"] = cat["Total_flux"] / (2.0 * cat["Maj"] * arcsectorad)

            # scale flux
            cat["Total_flux"] = cat["Total_flux"] * config.getfloat(
                "skymodel", "flux_factor"
            )
            cat["Peak_flux"] = cat["Peak_flux"] * config.getfloat(
                "skymodel", "flux_factor"
            )

            # scale size
            cat["Maj"] = cat["Maj"] * config.getfloat("skymodel", "sizefactor")
            cat["Maj_halflight"] = cat["Maj_halflight"] * config.getfloat(
                "skymodel", "sizefactor"
            )
            cat["Min"] = cat["Min"] * config.getfloat("skymodel", "sizefactor")
            cat["Min_halflight"] = cat["Min_halflight"] * config.getfloat(
                "skymodel", "sizefactor"
            )

            ix_arr = np.ones(nobj)
            iy_arr = np.ones(nobj)
            iv_arr = np.ones(nobj)

            unresolveds = 0

            # fov cut, put cos(dec) factor into ra offset
            cosdec = np.cos(dec_field_gs * 2 * np.pi / 360)

            # need to think about removing cosdec since higher-up sources are not filling plane
            ra_offset_max = (1 / cosdec) * (
                (fov / 60) / 2
            )  # this is now a fraction of the written image size
            dec_offset_max = (fov / 60) / 2  # convert fov to degrees

            fov_cut = (abs(cat["ra_offset"]) < ra_offset_max) * (
                abs(cat["dec_offset"]) < dec_offset_max
            )

            cat = cat[fov_cut]

            # Draw the galaxies onto the galsim image
            for i, cat_gal in enumerate(cat):

                logging.info(
                    "..........Adding source {0} of {1} to skymodel..........".format(
                        i + 1, nobj
                    )
                )

                # ANNA: not sure if this does anything as called before the postage stamp is created. I have inserted the convolution step inside the sub_img.
                """
                if config.getboolean('skymodel', 'dosimple_psf'):
                    psf_maj = config.getfloat('skymodel', 'simple_psf_maj')*galsim.arcsec
                    psf_min = config.getfloat('skymodel', 'simple_psf_min')*galsim.arcsec
                    psf_pa = config.getfloat('skymodel', 'simple_psf_pa')*galsim.degrees
                    q = (psf_min/galsim.arcsec)/(psf_maj/galsim.arcsec)
                    psf = galsim.Gaussian(fwhm=psf_maj/galsim.arcsec)
                    psf_shear = galsim.Shear(q=q, beta=psf_pa)
                    psf = psf.shear(psf_shear)
                    
                    gal = galsim.Convolve(gal, psf)
            """

                x, y = w_twod.wcs_world2pix(
                    cat_gal["RA"],
                    cat_gal["DEC"],
                    1,
                )

                x = float(x)
                y = float(y)
                # print(cat_gal["RA"],
                #      cat_gal["DEC"],x,y)
                # exit()
                #    for three axes:   pixels = wcs.world_to_pixel(coord, 3000 * u.m / u.s)

                logging.info("RA, Dec: %f %f ", cat_gal["RA"], cat_gal["DEC"])
                # logging.info('freq, z: %f %f',  cat_gal['central_freq'], cat_gal['z'])
                logging.info("x, y,: %f %f ", x, y)

                # ! determine whether need this offset bit - is it so that it's not all shifted by 0.5 pix??
                # Account for the fractional part of the position:
                ix = int(np.floor(x + 0.5))
                iy = int(np.floor(y + 0.5))
                ix_arr[i] = ix
                iy_arr[i] = iy

                # iv = int(np.floor(v+0.5))
                # iv_arr[i] = iv

                offset = galsim.PositionD(x - ix, y - iy)  # AB: original
                # offset is used by galsim

                logging.info("Continuum size from cat:  %f", cat_gal["Maj"])

                # Create the sub-image for this galaxy

                #            if source centre is in FoV:

                if (ix > 0) and (iy > 0) and (ix <= n_x) and (iy <= n_y):
                    print(
                        "source",
                        cat_gal["PA"],
                        cat_gal["Total_flux"],
                        cat_gal["Maj"],
                        cat_gal["Min"],
                        cat_gal["spectral_index"],
                        ix,
                        iy,
                    )

                    #### get the postage for the source. It can be AGN from library, Gaussian lobe and Gaussian core, Sersic of simple Gaussian

                    sub_img, atlas_source, unresolved, flux = make_img(
                        config,
                        cat_gal["Total_flux"],
                        cat_gal["spectral_index"],
                        base_freq,
                        freqs,
                        cat_gal["Maj"],
                        cat_gal["Min"],
                        cat_gal["PA"],
                        float(pixel_scale_str),
                        psf_maj_arcsec,
                        cat_gal["RadioClass"],
                        cat_gal["corefrac"],
                        cat_gal["ranid"],
                    )
                    unresolveds += unresolved
                    sub_img_shape = sub_img.shape

                    sub_img_size = sub_img.shape[1]
                    logging.info("postage stamp size %f", sub_img_size)
                    print(i)
                    if 1:  # (sub_img_size <= 300.):    #this is for memory issues!!!
                        print(i)
                        # have set cube to correct resolution in make_cube by reading pixel_scale from config

                        # works out the bounds for the postage stamp in the FoV image
                        l_bounds = np.array([0, y, x]) - np.array(
                            [0, (sub_img.shape[1] / 2), (sub_img.shape[2] / 2)]
                        )
                        u_bounds = np.array([nfreqs, y, x]) + np.array(
                            [
                                0,
                                sub_img_shape[1] - (sub_img.shape[1] / 2),
                                sub_img_shape[2] - (sub_img.shape[2] / 2),
                            ]
                        )

                        logging.info(
                            "Lower bounds, upper bounds: %s, %s", l_bounds, u_bounds
                        )

                        l_bounds = np.floor(l_bounds).astype(np.int)
                        u_bounds = np.floor(u_bounds).astype(np.int)

                        logging.info(
                            "Lower bounds, upper bounds, int: %s, %s",
                            l_bounds,
                            u_bounds,
                        )
                        logging.info("Subcube shape: %s", sub_img.shape)

                        #### add it to the large cube:

                        img3 = sub_img
                        blc0 = l_bounds[0]
                        blc1 = l_bounds[1]
                        blc2 = l_bounds[2]
                        trc0 = u_bounds[0] - 1
                        trc1 = u_bounds[1] - 1
                        trc2 = u_bounds[2] - 1

                        print(blc0, blc1, blc2, trc0, trc1, trc2)

                        trcs = np.array([trc0, trc1, trc2])

                        blcs = np.array([blc0, blc1, blc2])
                        # the top bounds are all -1 the true values, since the large cube is added using Fortran indexing

                        top_excess = arr_dims - (
                            trcs + 1
                        )  # pixels from the image to top coordinates of field
                        bottom_excess = blcs
                        excess = np.hstack((bottom_excess, top_excess))
                        print(excess)
                        overlap = False  # initialise indicator to say if the galaxy clips the edges
                        #  address the below
                        ### the galaxy is written only if it fits in its entirety - not clipped. This is possibly not desiderable as AGN can be big
                        for coord in excess:
                            if coord < 0:
                                overlap = True
                                logging.info(
                                    "Subcube is overlapping the edge: cropping to fit"
                                )
                                print(
                                    "Subcube is overlapping the edge: cropping to fit"
                                )
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

                        logging.info(
                            "BLC, TRC: %f %f %f, %f %f %f ",
                            blc0,
                            blc1,
                            blc2,
                            trc0,
                            trc1,
                            trc2,
                        )

                        # write the info for this object to files
                        fitsf = FITS(all_gals_fname, "rw")
                        fitsf_f = FITS(all_gals_fname + "_maxflux.fits", "rw")
                        fitsf_z = FITS(all_gals_fname + "_z.fits", "rw")

                        # the redshift map contains the redshift of the brightest source on the LoS. This is judged by looking at the dummy map _maxflux and comparing if with the postage stamp. The z map is updated only where the postage stamp is brighter than what recorder in _maxflux

                        # read the recorded values for flux and redshift at the postage location
                        flux_old = fitsf_f[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]
                        z_old = fitsf_z[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]

                        # initialise the new arrays
                        flux_new = z_old * 0.0
                        print(flux_new.shape)
                        flux_new[0] = img3[0]  # at the lowest frequency
                        zvalue = cat_gal["z"]

                        img_z = z_old
                        img_f = flux_old

                        # if (np.sum(flux_new)>np.sum(flux_old)):
                        #    img_z=img_z*0.+zvalue
                        #    img_f=flux_new

                        # update only where postage brighter than record
                        img_z[flux_new > flux_old] = zvalue
                        img_f[flux_new > flux_old] = flux_new[flux_new > flux_old]

                        fitsf_f[0].write(img_f, 0, blc1, blc2, 0, trc1, trc2)
                        fitsf_z[0].write(img_z, 0, blc1, blc2, 0, trc1, trc2)

                        # adding the source to the total map
                        region = fitsf[0][
                            blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1
                        ]
                        #             print(region.shape)
                        print(img3.shape)
                        img3 += region
                        fitsf[0].write(img3, blc0, blc1, blc2, trc0, trc1, trc2)

                        t_source = time.time() - tstart

                        fitsf.close()

                        logging.info("")

        print("loop finished")
        # write out catalogue
        truthcat_name = (
            data_path + config.get("field", "fits_prefix") + "_truthcat.fits"
        )
        logging.info("Writing truth catalogue to: {0} ...".format(truthcat_name))
        cat.write(truthcat_name, format="fits", overwrite=True)

        logging.info("Unresolveds: %d", unresolveds)

        tend = time.time()
        logging.info("...done in {0} seconds.".format(tend - tstart))


if __name__ == "__main__":

    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])

    runSkyModel(config)
