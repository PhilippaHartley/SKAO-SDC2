"""
Script to convert a T-RECS catalogue into a continuum sky model FITS file.

Usage:
python skymodel.py example.ini
"""


import logging
import multiprocessing
import os
import sys
import time
import galsim
import numpy as np
from astropy import units as uns
from astropy import wcs as ast_wcs
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
from astropy.io import fits as astfits
from astropy.table import Table
from fitsio import FITS, FITSHDR
from numpy.core.defchararray import add as stradd
from numpy.core.defchararray import multiply as strmultiply

import skymodel.skymodel_tools as tools
from skymodel.continuum_morphology import make_img
from skymodel.skymodel_tools import setup_wcs


tstart = time.time()

# mother_seed=5820743 #seed for random number generation - fullcube
mother_seed = 6879432  # seed for random number generation - smallcube1 (dev)
# mother_seed=7984532 #seed for random number generation -smallcube2 (eval)

arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value

"""
def initialise_file(filename, ):
    
def make_header(header_type,):

def retrieve_flux():
"""


def log_result(result):

    global cat


    (i, atlas_source, flux, unresolved) = result
    print (i)
    cat["id"][i] = i
    cat["Atlas_source"][i] = atlas_source
    cat["New_flux"][i] = flux
    cat["Unresolved"][i] = unresolved
    print (i)




def add_source_continuum(
    i,
    cat_gal,
    nobj,
    w_twod,
    config,
    pixel_scale_str,
    psf_maj_arcsec,
    arr_dims,
    all_gals_fname,
    base_freq,
    freqs,
):
    """
    mainlog = logging.getLogger("main%d" % i)
    h = logging.FileHandler("log%d.log" % i)
    mainlog.addHandler(h)
    logging.root.setLevel(logging.DEBUG)
    mainlog.info("result%s" % i)
    """
   

    logging.info(
        "..........Adding source {0} of {1} to skymodel..........".format(i + 1, nobj)
    )

    x, y = w_twod.wcs_world2pix(
        cat_gal["RA"],
        cat_gal["DEC"],
        1,
    )

    x = float(x)
    y = float(y)

    logging.info("RA, Dec: %f %f ", cat_gal["RA"], cat_gal["DEC"])
    # logging.info('freq, z: %f %f',  cat_gal['central_freq'], cat_gal['z'])
    logging.info("x, y,: %f %f ", x, y)

    logging.info("Continuum size from cat:  %f", cat_gal["Maj"])

    # get the postage for the source
    # it can be AGN from library, Gaussian lobe and Gaussian core, Sersic of simple Gaussian
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

    sub_img_shape = sub_img.shape
    sub_img_size = sub_img.shape[1]
    logging.info("postage stamp size %f", sub_img_size)

    # works out the bounds for the postage stamp in the FoV image
    l_bounds = np.array([0, y, x]) - np.array(
        [0, (sub_img.shape[1] / 2), (sub_img.shape[2] / 2)]
    )
    u_bounds = np.array([len(freqs), y, x]) + np.array(
        [
            0,
            sub_img_shape[1] - (sub_img.shape[1] / 2),
            sub_img_shape[2] - (sub_img.shape[2] / 2),
        ]
    )

    logging.info("Lower bounds, upper bounds: %s, %s", l_bounds, u_bounds)

    l_bounds = np.floor(l_bounds).astype(np.int)
    u_bounds = np.floor(u_bounds).astype(np.int)

    logging.info(
        "Lower bounds, upper bounds, int: %s, %s",
        l_bounds,
        u_bounds,
    )
    logging.info("Subcube shape: %s", sub_img.shape)

    # add it to the large cube
    img3 = sub_img
    blc0 = l_bounds[0]
    blc1 = l_bounds[1]
    blc2 = l_bounds[2]
    trc0 = u_bounds[0] - 1
    trc1 = u_bounds[1] - 1
    trc2 = u_bounds[2] - 1

    # the top bounds are all -1 the true values, since the large cube is added using Fortran indexing
    trcs = np.array([trc0, trc1, trc2])
    blcs = np.array([blc0, blc1, blc2])

    # pixels from the image to top coordinates of field
    top_excess = arr_dims - (
        trcs + 1
    )  
    bottom_excess = blcs
    excess = np.hstack((bottom_excess, top_excess))
    # initialise indicator to say whether the galaxy overlaps an edge
    overlap = False  
    # the galaxy is clipped if it overlaps an edge
    for coord in excess:
        if coord < 0:
            overlap = True
            logging.info("Subcube is overlapping the edge: cropping to fit")            
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

    # the redshift map contains the redshift of the brightest source on the LoS. 
    # this is judged by looking at the dummy map _maxflux and comparing if with the postage stamp. 
    # the z map is updated only where the postage stamp is brighter than what recorder in _maxflux

    # read the recorded values for flux and redshift at the postage location
    flux_old = fitsf_f[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]
    z_old = fitsf_z[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]

    # initialise the new arrays
    flux_new = z_old * 0.0
    flux_new[0] = img3[0]  # at the lowest frequency
    zvalue = cat_gal["z"]
    img_z = z_old
    img_f = flux_old

    # update only where postage brighter than record
    img_z[flux_new > flux_old] = zvalue
    img_f[flux_new > flux_old] = flux_new[flux_new > flux_old]

    fitsf_f[0].write(img_f, 0, blc1, blc2, 0, trc1, trc2)
    fitsf_z[0].write(img_z, 0, blc1, blc2, 0, trc1, trc2)

    # adding the source to the total map
    # if running in parallel, this step can go out of synch
    region = fitsf[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
    img3 += region
    fitsf[0].write(img3, blc0, blc1, blc2, trc0, trc1, trc2)

    t_source = time.time() - tstart

    fitsf.close()

    logging.info("")

    return (i, atlas_source, flux, unresolved)

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
    global cat
    
    # continuum specific
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

    # calling the appropriate wcs for continuum.
    # based on header_4D plus extra keywords
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

    # initialse empty cubes
    all_gals_fname = (
        data_path_large_files + config.get("field", "fits_prefix") + ".fits"
    )

    os.system("rm {0}".format(all_gals_fname))
    os.system("rm {0}".format(all_gals_fname + "_z.fits"))
    os.system("rm {0}".format(all_gals_fname + "_maxflux.fits"))
    logging.info("Creating empty image file, {0} ...".format(all_gals_fname))

    # exploting the expand_if_needed functionality of the image write function in fitsio

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

    # start pixel coords need to be supplied in reverse order when expanding - fiddly but works
    # fitsio code changed to reverse axes ordering when writing (going from C to fortran order)

    blc0 = image_size - 1  # -1 since the large cube is added using Fortran indexing
    blc1 = image_size - 1
    blc2 = n_chan - 1

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
    fitsf.close()
    cr_x, cr_y = w_twod.wcs_world2pix(
        ra_field_gs,
        dec_field_gs,
        1,
    )
    logging.info("Check world2pix crpix1, crpix2: %f %f", cr_x, cr_y)

    # Load the catalogue
    HI_cross = False  # assume that the catalogue is con cross-matched with HI
    cat_file_name = config.get("field", "catalogue")
    logging.info("Loading catalogue from {0} ...".format(cat_file_name))
    cat = Table()
    cat_read = Table.read(cat_file_name)  # remove ascii
    keywords = cat_read.colnames

    if "MHI" in keywords:
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
    cat["id"] = np.zeros(len(cat_read))
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

    # select continuum sources
    # exclude too big; memory problem and not realistic
    cat = cat[(cat["RadioClass"] != -100) * (cat["Maj"] < 200.0)]

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

    cat["Unresolved"] = np.zeros(len(cat)).astype(np.str)
    cat["New_flux"] = np.zeros(len(cat)).astype(np.str)
    np.random.seed(mother_seed + 100)

    # initialise core fraction. Steep-spectrum AGN dont use it as it is determined by the postage stamp.
    corefrac = np.random.normal(
        loc=0.75, scale=0.1, size=len(cat)
    )  
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
        scale_radius_to_hlr = 1.67834699
        cat["Maj_halflight"] = cat["Maj"] * scale_radius_to_hlr
        cat["Maj_halflight"].unit = "arcsec"
        cat["Min"] = cat["Maj"] * cat["q"]
        cat["Min"].unit = "arcsec"
        cat["Min_halflight"] = cat["Maj_halflight"] * cat["q"]
        cat["Min_halflight"].unit = "arcsec"
        # number of sources, on grid if requested

    nobj = len(cat)

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


    # fov cut, put cos(dec) factor into ra offset
    cosdec = np.cos(dec_field_gs * 2 * np.pi / 360)
    ra_offset_max = (1 / cosdec) * (
        (fov / 60) / 2
    )  
    dec_offset_max = (fov / 60) / 2  # convert fov to degrees
    fov_cut = (abs(cat["ra_offset"]) < ra_offset_max) * (
        abs(cat["dec_offset"]) < dec_offset_max
    )
    cat = cat[fov_cut]
    cat = cat[
        "id",
        "RA",
        "DEC",
        "Total_flux",
        "z",
        "spectral_index",
        "Maj",
        "Min",
        "PA",
        "RadioClass",
        "corefrac",
        "ranid",
        "Atlas_source",
        "New_flux",
        "Unresolved",
    ]


    print("going into loop")
    multiprocessing.set_start_method("fork")
    pool = multiprocessing.Pool(n_cores)
    for i, cat_gal in enumerate(cat):
        p = 0
        """
        mainlog = logging.getLogger("main%d" % i)
        h = logging.FileHandler("log%d.log" % i)
        mainlog.addHandler(h)
        logging.root.setLevel(logging.DEBUG)
        mainlog.info("test%s" % i)
        """
        

        pool.apply_async(
            add_source_continuum,
            args=(
                i,
                cat_gal,
                nobj,
                w_twod,
                config,
                pixel_scale_str,
                psf_maj_arcsec,
                arr_dims,
                all_gals_fname,
                base_freq,
                freqs,
            ), callback=log_result,
        )

    pool.close()
    pool.join()
    print (cat)
    print (cat["Atlas_source"])  
    print  (cat["ranid"])
    print (cat["id"])    
    print(time.time() - tstart)

    # write out continuum catalogue
    truthcat_name = (
        data_path + config.get("field", "fits_prefix") + "_truthcat.fits"
    )
    logging.info("Writing truth catalogue to: {0} ...".format(truthcat_name))
    cat.write(truthcat_name, format="fits", overwrite=True)

    flux = astfits.getdata(all_gals_fname) 
    summed_line_flux = np.sum(flux, axis=0)
    print (np.sum(summed_line_flux))
    
    tend = time.time()
    logging.info("...done in {0} seconds.".format(tend - tstart))


if __name__ == "__main__":

    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])

    runSkyModel(config)
