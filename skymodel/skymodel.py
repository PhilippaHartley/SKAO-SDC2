"""
Script to convert a T-RECS catalogue into a sky model FITS file.

"""


import logging
import multiprocessing
import os
import sys
import time
import galsim
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uns
from astropy.cosmology import LambdaCDM
from astropy.io import fits as astfits
from astropy.table import Table
from fitsio import FITS, FITSHDR
from numpy.core.defchararray import add as stradd
from numpy.core.defchararray import multiply as strmultiply


import skymodel.skymodel_tools as tools
from skymodel.continuum_morphology import make_img
from skymodel.HI_morphology import make_cube
from skymodel.skymodel_tools import setup_wcs

tstart = time.time()

arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value

"""
def initialise_file(filename, ):
    
def make_header(header_type,):

def retrieve_flux():
"""


def log_result(result):
    global cat
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
    # add new property values to catalogue
    cat["id"][i] = i
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

):
    '''
    print("making source")
    mainlog = logging.getLogger("main%d" % i)
    h = logging.FileHandler("log%d.log" % i)
    mainlog.addHandler(h)
    logging.root.setLevel(logging.DEBUG)
    mainlog.info("result%s" % i)
    '''

    logging.info(
        "..........Adding source {0} of {1} to skymodel..........".format(i + 1, nobj)
    )

    print ('source', i)
    x, y, v = w_spectral.wcs_world2pix(
        cat_gal["RA"], cat_gal["Dec"], cat_gal["opt_vel"], 1
    )

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
    (
        sub_cube,
        D_HI,
        D_HI_arcsec,
        atlas_source,
        flux,
        w20,
        MHI_incl_scale,
        MHI_incl_v_scale,
        ska_PA,
    ) = make_cube(
        config,
        i,
        10 ** cat_gal["MHI"],
        cat_gal["z"],
        cat_gal["i"],
        cat_gal["PA"],
        float(pixel_scale_str),
        dnu,
        psf_maj_arcsec,
        cat_gal["RadioClass"],
    )

    sub_cube_shape = sub_cube.shape

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
    fitsf[0].write(img3, blc0, blc1, blc2, trc0, trc1, trc2)
    fitsf.close()
    logging.info("")

    '''
    mainlog.info("done_make_cube%s" % i)
   
    with open("log%d.log" % i, "r") as f:
        a = f.readlines()
        print(a)

    with open("all_log.txt", "a") as f:
        f.writelines(a)
    os.system("rm log%d.log" % i)
    '''

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

    global cat

    # set up cosmology
    H = config.getfloat("cosmology", "H")
    M = config.getfloat("cosmology", "M")
    L = config.getfloat("cosmology", "L")
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

    header_spectral = w_spectral.to_header()
    header_twod = w_twod.to_header()
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

    # fov cut
    cosdec = np.cos(dec_field_gs * 2 * np.pi / 360)
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
    # cut according to frequency range
    vel_cut = ((cat["opt_vel"]) > velo_base) * ((cat["opt_vel"]) < velo_max)
    cat = cat[vel_cut]

    # retain only sources with HI signal
    MHI_cut = cat["MHI"] > 0
    cat = cat[MHI_cut]

    logging.info("Cat length after fov and z cut: %d", len(cat))

    nobj = len(cat)

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

    if pool:
        print("going into loop")
        multiprocessing.set_start_method("fork")
        pool = multiprocessing.Pool(n_cores)
        for i, cat_gal in enumerate(cat):
            '''
            mainlog = logging.getLogger("main%d" % i)
            h = logging.FileHandler("log%d.log" % i)
            mainlog.addHandler(h)
            logging.root.setLevel(logging.DEBUG)
            mainlog.info("test%s" % i)
            '''

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
                    
                ),
                callback=log_result,
            )

        pool.close()
        pool.join()
        print(cat)
        print(cat["Atlas_source"][i])



    # check all sources have been created and added
    filled_rows = np.argwhere(cat["Atlas_source"] != "0.0")[:,0]  
    len1 = len(cat)
    cat = cat[filled_rows]
    len2 = len(cat)
    if len1 != len2:
        print("warning: some sources were not added")
        print('input sources', len1,' vs output sources', len2)
        exit()


    print (cat["id"])
    cat["id"] = np.arange(len(cat))
    print (cat["id"])

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
    public_cat["w20"] = cat["w20"]
    public_cat["w20"].unit = "km/s"
    public_truthcat_name = (
        data_path + config.get("field", "fits_prefix") + "_public_truthcat.fits"
    )
    logging.info(
        "Writing public truth catalogue to: {0} ...".format(public_truthcat_name)
    )
    public_cat.write(public_truthcat_name, format="fits", overwrite=True)

    logging.info("Writing summed cube to: {0} ...".format(all_gals_summed_fname))
    # save a 2d image using astropy instead of fitsio
    flux = astfits.getdata(all_gals_fname)
    summed_line_flux = np.sum(flux, axis=0)
    hdu = astfits.PrimaryHDU(summed_line_flux, header=header_twod)
    hdulist = astfits.HDUList([hdu])
    hdulist.writeto(all_gals_summed_fname, overwrite=True)

    tend = time.time()
    logging.info("...done in {0} seconds.".format(tend - tstart))
    print(tend - tstart)

    print("preparing check by doing simple sum of abs values")
    print("sum of absolute 3d array values is %f" % np.abs(flux).sum())
    print("preparing checksum")
    os.system(
        "md5sum %s > %s" % (all_gals_fname, all_gals_fname.strip("fits") + "md5")
    )
    

if __name__ == "__main__":

    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])

    runSkyModel(config)
