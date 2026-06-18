'''
    Filename: compute_psf_matched_maps.py
    Notes: Computes effective radius for all objects in a given PASSAGE field/s, and stores the result as a FITS table
           Run this code within pyt312 conda environment
    Author : Ayan
    Created: 16-01-26
    Example: run compute_psf_matched_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run compute_psf_matched_maps.py --system ssd --field Par28 --id 2865 --clobber --debug_psf
             run compute_psf_matched_maps.py --system ssd --field Par28 --do_all_obj --clobber
             run compute_psf_matched_maps.py --system ssd --do_all_fields --do_all_obj --clobber
'''

import os
import re
import sys
import numpy as np
import argparse
import tempfile
import gc
import subprocess
import pandas as pd
from scipy.ndimage import rotate, shift
from astropy.io import fits
from drizzlepac.astrodrizzle import adrizzle
from grizli.utils import get_wcs_pscale, transform_wcs
from astropy.convolution import convolve, Box2DKernel
from astropy.table import Table
from astropy.time import Time
from astropy import wcs as pywcs
from astropy import units as u
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from pathlib import Path

import warnings
warnings.filterwarnings("ignore")

import stpsf
stpsf.setup_logging(level='ERROR') # to suppress chatty-ness of webbpsf

from datetime import datetime, timedelta
start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')

    parser.add_argument('--input_dir', metavar='input_dir', type=str, action='store', default=None, help='Where do the input files reside?')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, action='store', default=None, help='Where do you want to store the outputs?')
    parser.add_argument('--system', metavar='system', type=str, action='store', default='hd', help='Which file system is the code being run on?')
    parser.add_argument('--code_dir', metavar='code_dir', type=str, action='store', default='/Users/acharyya/Work/astro/ayan_codes/passage/', help='Where is the source code?')
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='Over-write existing plots? Default is no.')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress some generic print statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')
    parser.add_argument('--drv', metavar='drv', type=str, action='store', default='v0.6', help='Which data reduction version? Default v0.6')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=15, help='fontsize of plot labels, etc.; default is 15')

    parser.add_argument('--field', metavar='field', type=str, action='store', default='Par3', help='Which passage field? Default is Par50')
    parser.add_argument('--id', metavar='id', type=str, action='store', default=None, help='Object ID. Default is None')

    parser.add_argument('--pixscale', metavar='pixscale', type=float, action='store', default=0.04, help='Pixel scale (in arcsec/pixel) of the thumbnails produced; default is 0.04')

    parser.add_argument('--arcsec_limit', metavar='arcsec_limit', type=float, action='store', default=-99, help='Half box size (in arcsec) of the thumbnails to plot; default is 1.5')
    parser.add_argument('--re_limit', metavar='re_limit', type=float, action='store', default=None, help='Effective radius to limit all analysis to; default is None (uses arcsec_limit)')
    parser.add_argument('--only_seg', dest='only_seg', action='store_true', default=False, help='Cut out the emission line plots corresponding to the grizli segmentation map? Default is no.')
    parser.add_argument('--test_cutout', dest='test_cutout', action='store_true', default=False, help='Plot the cutout 2D clear image as a testing phase? Default is no.')
    parser.add_argument('--no_text_on_plot', dest='no_text_on_plot', action='store_true', default=False, help='Skip putting text annotations on plot2D? Default is no.')

    parser.add_argument('--kernel_mode', metavar='kernel_mode', type=str, action='store', default='linear_interp', help='Mode to be used for astropy Box2DKernel? Default is linear_interp')
    parser.add_argument('--kernel_size', metavar='kernel_size', type=float, action='store', default=3, help='Kernel size to be used for astropy Box2DKernel? Default is 3')

    parser.add_argument('--debug_offset', dest='debug_offset', action='store_true', default=False, help='Do extra plots and prints for debugging offset calculation from center? Default is no.')
    parser.add_argument('--debug_psf', dest='debug_psf', action='store_true', default=False, help='Debug the PSF-matching step? Default is no.')
    parser.add_argument('--plot_psf', dest='plot_psf', action='store_true', default=False, help='Save plots for the PSF-matching step? Default is no.')
    parser.add_argument('--do_all_fields', dest='do_all_fields', action='store_true', default=False, help='Include ALL available fields? Default is no.')
    parser.add_argument('--do_all_obj', dest='do_all_obj', action='store_true', default=False, help='Reduce spectra and make beam files for ALL detected objects? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()

    args.field_arr = args.field.split(',')
    for index in range(len(args.field_arr)):
        if 'Par' in args.field_arr[index]: args.field_arr[index] = f'Par{int(args.field_arr[index].split("Par")[1]):03d}'
    args.field = args.field_arr[0]

    if args.id is not None: args.id = [int(item) for item in args.id.split(',')]
    if args.system == 'hd' and not os.path.exists('/Volumes/Elements/'): args.system = 'local'
    if args.system == 'ssd' and not os.path.exists('/Volumes/Ayan_SSD/'): args.system = 'local'

    if 'local' in args.system: root_dir =  '/Users/acharyya/Work/astro/passage'
    elif 'hd' in args.system: root_dir = '/Volumes/Elements/acharyya_backup/Work/astro/passage'
    elif 'ssd' in args.system: root_dir = '/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage'
    elif 'gdrive' in args.system: root_dir = '/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/passage'
    else: root_dir = ''
    args.root_dir = Path(root_dir)

    if 'glass' in args.field: survey_name = 'glass'
    else: survey_name = 'passage'
    if args.input_dir is None:
        args.input_dir = args.root_dir / f'{survey_name}_data/'
    if args.output_dir is None:
        args.output_dir = args.root_dir / f'{survey_name}_output/'

    if 'glass' in args.field:
        args.drv = args.glass_version
    elif 'Par' in args.field and 'v' not in args.drv:
        args.drv = 'v' + args.drv

    args.input_dir = Path(args.input_dir) / args.drv
    args.output_dir = Path(args.output_dir) / args.drv
    (args.output_dir / 'catalogs').mkdir(exist_ok=True, parents=True)
    (args.output_dir / 'plots').mkdir(exist_ok=True, parents=True)

    args.input_dir = Path(args.input_dir)
    args.output_dir = Path(args.output_dir)
    args.code_dir = Path(args.code_dir)

    if args.arcsec_limit == -99: # i.e. it has not been explicitly set by the user
        if args.debug_psf:
            args.arcsec_limit = 0.2
        else:
            args.arcsec_limit = 1.0 # default value

    return args

# ----------------------------------------------------------------
def setup_plot_style():
    '''
    Function to set default style for all plots made in this project
    '''
    # plt.rcParams['pdf.fonttype']	= 42
    # plt.rcParams['ps.fonttype'] 	= 42
    # plt.rcParams['savefig.dpi'] 	= 600
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, fig_dir, figname, args, dpi=100):
    '''
    Saves a given figure handle as a given output filename
    '''
    fig_dir.mkdir(exist_ok=True, parents=True)
    figname = fig_dir / figname
    fig.savefig(figname, transparent=False, dpi=dpi)
    print(f'\t\tSaved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def read_passage_sed_catalog(filename):
    '''
    Read the combined master catalog from PASSAGE SED fits, rename a few columns, and only keep the mass and SFR columns
    Return pandas dataframe
    '''
    if Path(filename).suffix == '.fits':
        print(f'Reading master PASSAGE SED catalog from {filename}..')
        full_df = Table.read(filename).to_pandas()
    else:
        print(f'Reading PASSAGE line finding catalog from {filename}..')
        full_df = pd.read_csv(filename, header=0, sep='\t')
        full_df = full_df[full_df['mass_50'] > 0].reset_index(drop=True) # to get only those sources that have stellar mass measured
        full_df = full_df.drop('id', axis=1)

    if 'cosmoswebid_1' in full_df: full_df.drop('cosmoswebid_1', axis=1, inplace=True)
    full_df.rename(columns={'Par':'field', 'passage_id':'id', 'id_photcat':'id', 'objid':'id', 'field_id':'id', 'z_best':'redshift', 
                            'zbest':'redshift', 'mass_50':'log_mass', 'stellar_mass_50':'log_mass', 'ssfr_50':'log_ssfr', 'sfr_50':'sfr', 
                            'ra_obj':'ra', 'dec_obj':'dec', 'cosmoswebid_2': 'cosmoswebid'}, inplace=True)
    
    # -----------computing new columns------------
    full_df['log_sfr'] = np.log10(full_df['sfr'])
    
    full_df['delta_tform_50_10'] = full_df['tform50_50'] - full_df['tform10_50']
    full_df['delta_tform_90_50'] = full_df['tform90_50'] - full_df['tform50_50']
    full_df['delta_tform_90_10'] = full_df['tform90_50'] - full_df['tform10_50']

    # --------extracting relevant columns-----------
    columns_to_extract = ['field', 'id', 'redshift', 'log_mass', 'log_sfr', 'log_ssfr', 'cosmoswebid', 'delta_tform_50_10', 'delta_tform_90_50', 'delta_tform_90_10']
    df = full_df[columns_to_extract]
    df['field'] = df['field'].astype(str)
    nobj1 = len(df)

    df = df.dropna(subset=['log_mass', 'log_sfr'])
    nobj2 = len(df)
    print(f'\nOut of initial {nobj1} objects, {nobj2} had mass available.')

    return df

# --------------------------------------------------------------------------------------------------------------------
def myimshow(data, ax, contour=None, re_pix=None, label='', cmap='viridis', fontsize=10, col='w', cmin=None, cmax=None):
    '''
    Utility function to plot a 2D data array on to ax, and plot center and segmentation map
    Returns ax
    '''
    p = ax.imshow(data, origin='lower', cmap=cmap, vmin=cmin, vmax=cmax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.0)
    plt.colorbar(p, cax=cax)
    cen_x, cen_y = int(np.shape(data)[0] / 2), int(np.shape(data)[1] / 2)
    #ax.scatter(cen_y, cen_x, marker='x', c='k')
    if contour is not None: ax.contour(contour, levels=0, colors=col, linewidths=0.5)
    ax.text(0.9, 0.9, label, c=col, fontsize=fontsize, ha='right', va='top', transform=ax.transAxes)
    if re_pix is not None: ax.add_patch(plt.Rectangle((-re_pix + cen_y, -re_pix + cen_x), 2 * re_pix, 2 * re_pix, lw=0.5, color='r', fill=False))

    return ax

# --------------------------------------------------------------------------------------------------------------------
def read_direct_image_from_extension(full_hdu, args, filter='F150W', for_offset=False):
    '''
    Reads in the direct image of an object in the given filter, by trying out a few combinations of how the filter name might be in the header
    Returns direct image and the final filter name it found in the header
    '''
    dummy_filter = 'F140W'
    try:
        dir_img = full_hdu['DSCI', f'{filter}-{filter}-CLEAR'].data
    except:
        try:
            dir_img = full_hdu['DSCI', f'{filter}-CLEAR'].data
        except:
            try: 
                dir_img = full_hdu['DSCI', filter].data
            except:
                try:
                    dir_img = full_hdu['DSCI', dummy_filter].data
                    filter = dummy_filter
                except:
                    try:
                        dir_img = full_hdu['DSCI', 'CLEAR'].data
                        filter = 'CLEAR'
                    except:
                        dir_img = None

    dir_img = trim_image(dir_img, args, skip_re_trim=True)

    return dir_img, filter

# --------------------------------------------------------------------------------------------------------------
def trim_image(image, args=None, arcsec_limit=None, pix_size_arcsec=None, re_limit=None, re_arcsec=None, skip_re_trim=False):
    '''
    Trim a given 2D image to a given arcsecond dimension
    Returns 2D map
    '''
    if args is not None:
        pix_size_arcsec = args.pix_size_arcsec
        if args.re_limit is None or skip_re_trim:
            arcsec_limit = args.arcsec_limit
        else:
            arcsec_limit = args.re_limit * args.re_arcsec
    elif re_limit is not None and re_arcsec is not None:
        arcsec_limit = re_limit * re_arcsec

    image_shape = np.shape(image)
    center_pix = int(image_shape[0] / 2.)
    farthest_pix = int(arcsec_limit / pix_size_arcsec) # both quantities in arcsec

    image = image[center_pix - farthest_pix : center_pix + farthest_pix, center_pix - farthest_pix : center_pix + farthest_pix]
    # print(f'Trimming image of original shape {image_shape} to {args.arcsec_limit} arcseconds, which is from pixels {center_pix - farthest_pix} to {center_pix + farthest_pix}, so new shape is {np.shape(image)}')

    return image

# --------------------------------------------------------------------------------------------------------------------
def get_offsets_from_center(full_hdu, args, filter='F200W', silent=False):
    '''
    Computes the offset from the original center of image to the brightest pixel in the direct image with the given filter
    Returns two integers (offset in x and y axes)
    '''
    dir_img, filter = read_direct_image_from_extension(full_hdu, args, filter=filter, for_offset=True)

    if dir_img is not None:
        segmentation_map = full_hdu['SEG'].data
        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)

        smoothing_kernel = Box2DKernel(5, mode=args.kernel_mode)
        dir_img_smoothed = convolve(dir_img, smoothing_kernel)
        smooth_shape = dir_img_smoothed.shape
        ncells = 5 # only searches within +/- 5 cells of the original center
        dir_img_smoothed_subarea = dir_img_smoothed[smooth_shape[0] // 2 - ncells: smooth_shape[0] // 2 + ncells, smooth_shape[1] // 2 - ncells: smooth_shape[1] // 2 + ncells] # only searching for brightest pixel in the vicinity of the original center, otherwise might pick up neighbouring galaxies
        brightest_coords = np.where(dir_img_smoothed == np.nanmax(dir_img_smoothed_subarea))
        brightest_x, brightest_y = brightest_coords[0][0], brightest_coords[1][0]
        cen_x, cen_y = int(np.shape(dir_img)[0] / 2), int(np.shape(dir_img)[1] / 2)
        ndelta_xpix = cen_x - brightest_x
        ndelta_ypix = cen_y - brightest_y
        #ndelta_xpix, ndelta_ypix = 0, 0 ##
        if not silent: print(f'For {args.field}:{args.id}: Determined x and y offsets from {filter} direct image = {ndelta_xpix}, {ndelta_ypix}')

        if args.debug_offset:
            print(f'Deb2934: original shape = {np.shape(dir_img)}') ##
            print(f'Deb2935: original center = {cen_x}, {cen_y}') ##
            print(f'Deb2934: brightest pixel = {brightest_x}, {brightest_y}') ##

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.07, wspace=0.3, hspace=0.1)
            cmap = 'viridis'
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: Centering offset diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

            axes[0].imshow(dir_img, origin='lower', cmap=cmap)
            axes[0].scatter(cen_y, cen_x, marker='x', c='k')
            axes[0].scatter(brightest_y, brightest_x, marker='x', c='r')
            axes[0].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[0].text(0.9, 0.9, 'Original', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[0].transAxes)

            axes[1].imshow(dir_img_smoothed, origin='lower', cmap=cmap)
            axes[1].scatter(cen_y, cen_x, marker='x', c='k')
            axes[1].scatter(brightest_y, brightest_x, marker='x', c='r')
            axes[1].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[1].text(0.9, 0.9, 'Smoothed', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[1].transAxes)

            dir_img_shifted = np.roll(dir_img_smoothed, ndelta_xpix, axis=0)
            dir_img_shifted = np.roll(dir_img_shifted, ndelta_ypix, axis=1)
            
            segmentation_map = np.roll(segmentation_map, ndelta_xpix, axis=0)
            segmentation_map = np.roll(segmentation_map, ndelta_ypix, axis=1)

            new_cen_x, new_cen_y = int(np.shape(dir_img_shifted)[0] / 2), int(np.shape(dir_img_shifted)[1] / 2)
            
            axes[2].imshow(dir_img_shifted, origin='lower', cmap=cmap)
            axes[2].scatter(new_cen_y, new_cen_x, marker='x', c='k')
            axes[2].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[2].text(0.9, 0.9, 'Shifted', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[2].transAxes)

            plt.show(block=False)
            sys.exit(f'Exiting here because of --debug_offset mode; if you want to run the full code as usual then remove the --debug_offset option and re-run')
    
    else:
        ndelta_xpix = 0
        ndelta_ypix = 0
        print(f'Could not find image for filter {filter} in {args.field}:{args.id}, so forcing center offsets to be {ndelta_xpix}, {ndelta_ypix}')

    return ndelta_xpix, ndelta_ypix

# --------------------------------------------------------------------------------------------------------------------
def get_obs_params_from_drz_img(filter_name, args, keywords_to_extract=['EXPSTART', 'PA_APER']):
    '''
    Extracts observation parameters (e.g. observing time and PA_APER) keywords from the header of the drz image file 
    of corresponding filter
    Returns keyword values
    '''
    keyword_values_list = []
    drz_filename = full_filename.parent.parent / f'{args.field}_{filter_name.lower()}-gr150r_drz_sci.fits'
    with fits.open(drz_filename, memmap=True) as beams_hdu:
        drz_header = beams_hdu[0].header
        for kwd in keywords_to_extract:
            if 'BEG' in kwd or 'START' in kwd: # this is a Time value
                value = Time(drz_header[kwd], format='mjd')  # + 'T' + header['TIME-OBS'])
            else:
                value = drz_header[kwd]
            keyword_values_list.append(value)
    del drz_header

    return keyword_values_list

# --------------------------------------------------------------------------------------------------------------------
def find_filter(wavelength, filter_dict):
    """
    Returns the first filter name where the wavelength falls inside the range.
    Returns None if no match is found.
    """
    for filt, (w_min, w_max) in filter_dict.items():
        if w_min <= wavelength <= w_max:
            return filt
    return None

# --------------------------------------------------------------------------------------------------------------------
def match_pypher(obs_psf_hdu, target_psf_hdu, reg=0.01, angle_source=0, angle_target=0):
    """
    Calls the pypher CLI from within Python to generate the kernel.
    Returns PSF-matched kernel
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Define paths for the temp FITS files
        obs_psf_path = os.path.join(tmpdir, 'obs_psf.fits')
        target_psf_path = os.path.join(tmpdir, 'target_psf_hdu.fits')
        output_kernel_path = os.path.join(tmpdir, 'matched_kernel.fits')

        # Save arrays as FITS files
        obs_psf_hdu.writeto(obs_psf_path, overwrite=True)
        target_psf_hdu.writeto(target_psf_path, overwrite=True)

        # run pypher
        cmd = ['pypher', obs_psf_path, target_psf_path, output_kernel_path, '-r', str(reg), '-s', str(angle_source), '-t', str(angle_target)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    
        # loading the resulting FITS file
        kernel = fits.getdata(output_kernel_path)
    
    return kernel

# --------------------------------------------------------------------------------------------------------------------
def get_niriss_psf(filter_name, args, supersampling_factor=1, is_target_psf=False, axes=None, img_wcs=None):
    '''
    Computes NIRISS PSF in a given filter
    Returns 2D array
    Borrowed in part from PJW
    '''
    fov_pixels = 2 * args.arcsec_limit / args.pix_size_arcsec
    if fov_pixels % 2 == 0: fov_pixels += 1
    niriss = stpsf.NIRISS()

    if is_target_psf: 
        pa_angle = 0
    else:
        pa_angle = args.pa_aper
        niriss.load_wss_opd_by_date(args.obs_date, verbose=False, plot=False, choice='closest')
    
    if type(filter_name) == str:
        print(f'\t\tCreating PSF at wavelength {filter_name} filter with fov_pixels={fov_pixels}..\n')
        niriss.filter = filter_name    
        psf = niriss.calc_psf(fov_pixels=fov_pixels, oversample=supersampling_factor)
    elif type(filter_name) == np.float64 or type(filter_name) == float: # here 'filter' is actually the observed wavelength in microns, not the filter
        print(f'\t\tCreating PSF at wavelength {filter_name} microns with fov_pixels={fov_pixels}..\n')
        psf = niriss.calc_psf(monochromatic=filter_name * 1e-6, fov_pixels=fov_pixels, oversample=supersampling_factor)
    else:
        sys.exit(f'Unrecognised data type for filter_name={filter_name}; it can only be str or float')

    psf_hdu  = psf['DET_DIST'] # choose from 'OVERSAMP', 'DET_SAMP', 'DETDIST', 'DET_DIST'
    psf_data = psf_hdu.data
    psf_wcs = pywcs.WCS(psf_hdu)

    # plots for testing/debugging
    if axes is not None:
        axes[0] = myimshow(psf_data, axes[0], contour=args.segmentation_map != args.id, label='Base PSF', cmap='viridis', col='w')

    # -------rotating the PSF-------------
    if img_wcs is not None:
        '''
        psf_data = rotate(psf_data, pa_angle, reshape=False)
        '''
        psf_wcs.wcs.crpix = (np.asarray(psf_data.shape) + 1) / 2
        psf_wcs.wcs.crval = [args.ra, args.dec]
        rotation_angle_rad = np.radians(pa_angle - 360)
        psf_wcs.wcs.cd = (np.array([
                    [np.cos(rotation_angle_rad), -np.sin(rotation_angle_rad)],
                    [np.sin(rotation_angle_rad), np.cos(rotation_angle_rad)],])
            * (niriss.pixelscale * u.arcsec).to(u.deg).value)
        psf_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        psf_wcs.pscale = get_wcs_pscale(psf_wcs)
        img_wcs.pscale = get_wcs_pscale(img_wcs)
        if img_wcs._naxis[0] % 2 == 0:
            img_wcs = transform_wcs(img_wcs, translation=(0.5, 0.5))
            img_wcs._naxis = np.asarray(img_wcs._naxis) + 1

        out_sci = np.zeros(img_wcs._naxis, dtype=np.float32)
        out_wht = np.zeros(img_wcs._naxis, dtype=np.float32)
        out_ctx = np.zeros(img_wcs._naxis, dtype=np.int32)

        adrizzle.do_driz(
            psf_data.astype(np.float32),
            psf_wcs,
            np.ones_like(psf_data, dtype=np.float32),
            img_wcs,
            outsci=out_sci,
            outwht=out_wht,
            outcon=out_ctx,
            expin=1.0,
            in_units="cps",
            wt_scl=1,
            stepsize=1,
        )
        psf_data = out_sci
    
    psf_hdu = fits.PrimaryHDU(data=psf_data, header=psf_wcs.to_header() )
    psf_hdu.header['PIXASEC'] = float(args.pix_size_arcsec)
    psf_hdu.header['PRSCALE'] = float(args.pix_size_arcsec)

    # plots for testing/debugging
    if axes is not None:
        psf_data_to_plot = trim_image(psf_data, args, skip_re_trim=True)
        axes[1] = myimshow(psf_data_to_plot, axes[1], contour=args.segmentation_map != args.id, label=f'Rotated PSF by {pa_angle:.1f} deg', cmap='viridis', col='w')

    return psf_hdu

# --------------------------------------------------------------------------------------------------------------------
def match_to_psf(image, obs_wave, psf_cache, target_psf_hdu, args, supersampling_factor=2, label='', img_wcs=None):
    '''
    Function to match a given 2D image from a given source PSF to taget_PSF of F200W
    Returns convolved image fo same shape
    '''
    lines_to_plot = ['OIII', 'Ha']
    make_psf_plot = args.debug_psf or (args.plot_psf and label.endswith('LINE') and any(line in label for line in lines_to_plot))
    if type(obs_wave) == np.float64 or type(obs_wave) == float:
        obs_wave = np.round(obs_wave, 3)

    # plots for testing/debugging
    if make_psf_plot:
        image_to_plot = trim_image(image, args, skip_re_trim=True)
        target_psf_to_plot = trim_image(target_psf_hdu.data, args, skip_re_trim=True)
        
        fig, axes = plt.subplots(2, 3, figsize=(8, 6), layout='constrained')
        axes[0][0] = myimshow(image_to_plot, axes[0][0], contour=args.segmentation_map != args.id, label=f'Original {label}', cmap='cividis', col='w')
        axes[0][1] = myimshow(target_psf_to_plot, axes[0][1], contour=args.segmentation_map != args.id, label='Target PSF', cmap='viridis', col='w')

    if obs_wave in psf_cache:
        obs_psf_hdu = psf_cache[obs_wave]
    else:
        obs_psf_hdu = get_niriss_psf(obs_wave, args, supersampling_factor=supersampling_factor, axes=[axes[0][2], axes[1][0]] if make_psf_plot else None, img_wcs=img_wcs)
        psf_cache[obs_wave] = obs_psf_hdu

    psf_kernel = match_pypher(obs_psf_hdu, target_psf_hdu, angle_target=0, angle_source=args.pa_aper)
    psf_matched_image = convolve(image, psf_kernel, boundary='extend')

    # plots for testing/debugging
    if make_psf_plot:
        psf_kernel_to_plot = trim_image(psf_kernel, args, skip_re_trim=True)
        psf_matched_image_to_plot = trim_image(psf_matched_image, args, skip_re_trim=True)

        axes[1][1] = myimshow(psf_kernel_to_plot, axes[1][1], contour=args.segmentation_map != args.id, label='matched kernel', cmap='viridis', col='w')
        axes[1][2] = myimshow(psf_matched_image_to_plot, axes[1][2], contour=args.segmentation_map != args.id, label=f'PSF-matched {label}', cmap='cividis', col='w')
        fig.text(0.01, 0.98, f'{args.field}: {args.id}', ha='left', va='top', fontsize=args.fontsize, color='k')
        save_fig(fig, args.input_dir / args.field / 'figs', f'{args.field}_{args.id:05d}_{label.replace(" ", "_")}_PSF_matching.png', args)
        
        if args.debug_psf:
            sys.exit(f'Stopping here since you used --debug_psf; remove that option to go through')
        else:
            plt.close('all')

    return psf_matched_image, psf_cache
    
# --------------------------------------------------------------------------------------------------------------------
def process_fits_extensions(full_hdu, target_psf_hdu, args):
    '''
    Processes each extension of a given fits file and PSF-matches it to a given target PSF
    Returns hdu list of new modified fits file, which can then be written to disk
    '''
    match_list = ['LINE', 'LINEWHT', 'CONTAM', 'CONTINUUM']
    new_hdul = fits.HDUList()

    psf_cache = {}

    # -------looping over extensions---------------
    for index, ext in enumerate(full_hdu):
        if args.debug_psf and "EXTVER" in ext.header and not (ext.header["EXTVER"] == 'Ha' and 'LINE' in ext.name.upper()):
            print(f'\tSkipping {ext.header["EXTVER"]} {ext_name} ({index + 1}/{len(full_hdu)}) because of the --debug_psf mode')
            continue
        ext_name = ext.name.upper()
        if any(substring in ext_name for substring in match_list):
            print(f'\tPSF-matching extension: {ext.header["EXTVER"]} {ext_name} ({index + 1}/{len(full_hdu)})')
            data_to_match = ext.data.copy()
            data_wcs = pywcs.WCS(ext)

            # ----------re-center line map-----------
            data_to_match = shift(data_to_match, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)

            # ----------determining filter for this extension----------
            obs_line_wave = ext.header['WAVELEN'] / 1e4 # to get wavelength in microns
            filt = find_filter(obs_line_wave, filter_waverange_dict)
            if filt is None:
                sys.exit(f'Deb526: filter is determined to None, based on obswave={obs_line_wave} mu')
            else:
                args.obs_date, args.pa_aper = get_obs_params_from_drz_img(filt, args, keywords_to_extract=['EXPSTART', 'PA_APER'])                    

            # ------account for LINEWHT----------------
            if ext_name.endswith('WHT'):
                print(f'\t\tConsidering {ext_name} extension as weights, taking inverse root of data..')
                data_to_match = 1. / np.sqrt(data_to_match)

            # -----------PSF-matching-------------
            psf_mached_image, psf_cache = match_to_psf(data_to_match, obs_line_wave, psf_cache, target_psf_hdu, args, label=f'{ext.header["EXTVER"]} {ext_name}', img_wcs=data_wcs)

            # ------account for LINEWHT----------------
            if ext_name.endswith('WHT'):
                print(f'\t\tReverting to inverse square of data..')
                psf_mached_image = 1. / (psf_mached_image)**2

            new_ext = fits.ImageHDU(data=psf_mached_image, header=ext.header, name=ext_name)            
        else:
            print(f'\tCopying untouched: {ext_name} ({index + 1}/{len(full_hdu)})')
            if isinstance(ext, fits.PrimaryHDU):
                new_ext = fits.PrimaryHDU(data=ext.data, header=ext.header)
                new_ext.header['PSF_MATCHED_TO_MICRON'] = target_psf_wave
            elif isinstance(ext, fits.BinTableHDU):
                new_ext = fits.BinTableHDU(data=ext.data, header=ext.header, name=ext_name)
            else:
                new_ext = fits.ImageHDU(data=ext.data, header=ext.header, name=ext_name)
                new_ext.header['PSF_MATCHED_TO_MICRON'] = target_psf_wave

        new_hdul.append(new_ext)
    
    return new_hdul

# --------------------------------------------------------------------------------------------------------------------
# Below: NIRISS filter data taken from Table 1 in https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-filters#gsc.tab=0
filter_waverange_dict = {'F115W': [1.0, 1.3], 'F150W': [1.31, 1.7], 'F200W': [1.71, 2.3]}  # microns

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    filter_for_re = 'F150W'

    # ---------determining list of fields----------------
    output_dir = args.output_dir / 'catalogs'
    passage_catalog_filename = output_dir / 'SED_fits_v1.0.2_cosmosweb.fits'
    df_sed = read_passage_sed_catalog(passage_catalog_filename)

    if not args.do_all_fields:
        df_sed = df_sed[df_sed['field'].isin(args.field_arr)]

    field_list = list(pd.unique(df_sed['field']))
    field_list.sort(key=natural_keys)

    # -------------setup PSF calculations------------------
    target_psf_wave = 2.2 # corresponding to F200W reddest wavelength, in microns
    get_target_psf = True

    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'Starting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'

        if not os.path.exists(product_dir):
            print(f'Could not find {product_dir}, so skipping this field {args.field}.')
            continue

        # ---------read the photometric catalog file--------------------
        if args.do_all_obj:
            #df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
            df = df_sed[df_sed['field'] == args.field]
            id_arr = df['id'].values
        else:
            id_arr = args.id
 
        # ------------looping over the objects-----------------------
        for index2, this_id in enumerate(id_arr):
            args.id = this_id
            print(f'\t\tCommencing ({index2 + 1}/{len(id_arr)}) ID {args.id}..')

            # ------determining directories and filenames---------
            full_fits_file = product_dir / 'full' / f'{args.field}_{args.id:05d}.full.fits'
            maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'

            if os.path.exists(maps_fits_file): # if the fits files are in maps/
                full_filename = maps_fits_file
            elif os.path.exists(full_fits_file): # if the fits files are in full/
                full_filename = full_fits_file
            else:
                print(f'\t\tCould not find {full_fits_file} or {maps_fits_file}, so skipping it.')
                continue

            outfiledir = Path(str(full_filename.parent).replace('maps', 'maps_psf_matched').replace('full', 'full_psf_matched'))
            outfiledir.mkdir(parents=True, exist_ok=True)
            outfilename = outfiledir / full_filename.name
            if os.path.exists(outfilename) and not args.clobber:
                print(f'\t\tResult file for {args.field}, {args.id} already exists as {outfilename}, so skipping this object.')
                continue

            # ------------read in maps files--------------------------------
            with fits.open(full_filename, memmap=True) as full_hdu:
                # ----------determining object parameters------------
                args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
                args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
                args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
                args.z = full_hdu[0].header['REDSHIFT']
                args.ra = full_hdu[0].header['RA']
                args.dec = full_hdu[0].header['DEC']

                if args.plot_psf or args.debug_psf:
                    args.segmentation_map = full_hdu['SEG'].data
                    args.segmentation_map = trim_image(args.segmentation_map, args, skip_re_trim=True)

               # --------determining true center of object---------------------
                args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter=filter_for_re)

                # -----------determine PSFs for PSF matching-----------------
                if get_target_psf:
                    # ----------getting observation keywords from drz file----------
                    args.obs_date, args.pa_aper = get_obs_params_from_drz_img('F200W', args, keywords_to_extract=['EXPSTART', 'PA_APER'])
                    img_wcs = pywcs.WCS(full_hdu['LINE'])                  
                    target_psf_hdu = get_niriss_psf(target_psf_wave, args)#, img_wcs=img_wcs) # not providing img_wcs to create target_psf becaus etarget_psf should not depend on the object/extension
                    get_target_psf = False

                # --------compute PSF-matched maps---------------------
                try:
                    new_hdul = process_fits_extensions(full_hdu, target_psf_hdu, args)
                    new_hdul.writeto(outfilename, overwrite=True)
                    print(f'\n\t\tSuccessfully saved matched FITS to: {outfilename}')
                
                except Exception as e:
                    if args.debug_psf:
                        raise
                    else:
                        print(f'\t\tCould not produce PSF-matched fits file for object {args.id} due to {e}. Skipping this object.')
                        continue
                
                del full_hdu
                del new_hdul
                gc.collect()
            
        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')
        
        if not args.debug:
            plt.close('all')
    
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
