'''
    Filename: plot_footprints_standalone.py
    Notes: Plots all PASSAGE field footprints and overlays with existing large survey footprints; this is a standalone code, i.e. does not depend on other utility scripts from this repo
    Author : Ayan
    Created: 05-03-2025
    Example: run plot_footprints_standalone.py --bg_image_dir /Users/acharyya/Work/astro/passage/passage_data/footprints/survey_images --bg_file COSMOS-HST-ACS_mosaic_Shrink100.fits
             run plot_footprints_standalone.py --fg_image_dir /Users/acharyya/Work/astro/passage/passage_data/COSMOS/imaging --fg_file mosaic_nircam_f444w_COSMOS-Web_60mas_B*_v0_8_sci.fits
             run plot_footprints_standalone.py --bg_image_dir /Users/acharyya/Work/astro/passage/passage_data/COSMOS/imaging --bg_file mosaic_nircam_f444w_COSMOS-Web_60mas_B7_v0_8_sci.fits --only_passage_regions --plot_objects_from_file /Users/acharyya/Work/astro/passage/passage_output/v0.5/catalogs/Par028_v0.5_venn_OIII,Ha,OII,Hb,SII,SNR>2.0_df.txt
'''

from __future__ import print_function
import argparse, glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from pathlib import Path
from astropy import wcs as pywcs
from astropy.io import fits
import regions
from regions import Regions
from astropy.coordinates import SkyCoord
import pandas as pd
from reproject import reproject_interp

import warnings
warnings.filterwarnings("ignore")

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')

    # ---- common args used widely over the full codebase ------------
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=10, help='fontsize of plot labels, etc.; default is 15')

    # ------- args added for plot_footprints.py ------------------------------
    parser.add_argument('--bg_image_dir', metavar='bg_image_dir', type=str, action='store', default='/Users/acharyya/Work/astro/passage/passage_data/footprints/survey_images', help='Which folder to be used for looking for the background image?')
    parser.add_argument('--bg_file', metavar='bg_file', type=str, action='store', default='COSMOS-HST-ACS_mosaic_Shrink100.fits', help='Which file to be used for plotting the background image?')
    parser.add_argument('--reg_files_dir', metavar='reg_files_dir', type=str, action='store', default='/Users/acharyya/Work/astro/passage/passage_data/footprints/region_files', help='Which folder to be used for looking for the region files?')
    parser.add_argument('--fg_image_dir', metavar='fg_image_dir', type=str, action='store', default='/Users/acharyya/Work/astro/passage/passage_data/COSMOS/imaging', help='Which folder to be used for looking for the foreground image?')
    parser.add_argument('--fg_file', metavar='fg_file', type=str, action='store', default=None, help='Which file to be used for plotting the foreground images?')
    parser.add_argument('--only_passage_regions', dest='only_passage_regions', action='store_true', default=False, help='Overplot ONLY the PASSAGE Par regions? Default is no.')
    parser.add_argument('--plot_objects_from_file', metavar='plot_objects_from_file', type=str, action='store', default=None, help='Overplot object locations from a list of coordinatess? Is fo, put the filename here')
    
    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()

    return args

# -------------------------------------------------------------------------------------------------------
def overplot_regions(region_files, bg_img_hdu, fig, fontsize=15):
    '''
    Overplots the footprint/s for a given region file and existing figure and header of the background file plotted on the figure
    Returns fig handle, and list of regions
    '''
    col_arr = ['blue', 'green', 'yellow', 'cyan']
    pix_offset_forlabels = 20

    ax = fig.gca()
    wcs_header = pywcs.WCS(bg_img_hdu[0].header)
    region_files = np.atleast_1d(region_files)

    for index, region_file in enumerate(region_files):
        print(f'Reading in region file {region_file}..')

        sky_regions = Regions.read(region_file, format='ds9')

        if 'color' in sky_regions[0].visual: color = sky_regions[0].visual['color']
        else: color = col_arr[index]

        for index2, sky_region in enumerate(sky_regions):
            if type(sky_region) == regions.shapes.text.TextSkyRegion: # label if it is text
                label = sky_region.text
                ax.text(ax.get_xlim()[0] + 0.1 * np.diff(ax.get_xlim())[0], ax.get_ylim()[1] * 0.98 - (index + 1) * 0.05 * np.diff(ax.get_ylim())[0], label, c=color, ha='left', va='top', fontsize=fontsize)
            else: # otherwise plot it
                # plotting the region
                pixel_region = sky_region.to_pixel(wcs_header)
                pixel_region.plot(ax=ax, lw=1, color=color)

                # labeling the region
                if type(sky_region) == regions.shapes.rectangle.RectangleSkyRegion:
                    label_pixcoord_x = pixel_region.center.xy[0] + pixel_region.width/2 + pix_offset_forlabels
                    label_pixcoord_y = pixel_region.center.xy[1] + pixel_region.height/2 + pix_offset_forlabels
                    label_text = f'P{pixel_region.meta["text"]}'
                    ax.text(label_pixcoord_x, label_pixcoord_y, label_text, c=color, ha='left', va='top', fontsize=fontsize/1.5)

    return fig, sky_regions

# -------------------------------------------------------------------------------------------------------
def plot_skycoord_from_df(df, fig, bg_img_hdu, color='aqua', alpha=0.3, size=1):
    '''
    Plots the location of all RA/DEC from a given dataframe on an existing background image, given the header of the background image
    Returns fig handle
    '''
    sky_coords = SkyCoord(df['ra'], df['dec'], unit='deg')

    wcs_header = pywcs.WCS(bg_img_hdu[0].header)
    ax=fig.gca()

    if len(df) < 20:
        # use this option for a slower plot, but with scale-able scatter points
        pix_coords = np.transpose(sky_coords.to_pixel(wcs_header))
        for index in range(len(df)):
            circle = plt.Circle(xy=pix_coords[index], radius=1, color=color, alpha=alpha)
            ax.add_patch(circle)
            if 'objid' in df: ax.text(pix_coords[index][0], pix_coords[index][1], df['objid'].values[index], color=color, fontsize=15)
    else:
        # use this option for a quicker plot, but not with scale-able scatter points
        ra_coords, dec_coords = sky_coords.to_pixel(wcs_header)
        ax.scatter(ra_coords, dec_coords, color=color, s=size, alpha=alpha, zorder=-1)

    return fig

# -------------------------------------------------------------------------------------------------------
def plot_background(filename, cmap='Greys', fontsize=15):
    '''
    Plots the background image for a given input filename
    Returns fig handle
    '''

    print(f'Reading in background image file {filename}')
    bg_img_hdu = fits.open(filename)
    sci_ext = 0 if 'COSMOS-Web' in str(filename) else 0
    data = np.log10(bg_img_hdu[sci_ext].data)
    header = bg_img_hdu[sci_ext].header

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(right=0.99, top=0.95, bottom=0.1, left=0.01)

    ax.imshow(data, origin='lower', cmap=cmap, vmin=-4, vmax=-3) # full, relevant range of data, for COSMOS field, is roughly (-6,-1)

    CDELT1 = 'CDELT1' if 'CDELT1' in header else 'CD1_1'
    CDELT2 = 'CDELT2' if 'CDELT2' in header else 'CD2_2'
    ra_offset = header['CRVAL1'] - header['CRPIX1'] * header[CDELT1]
    ra_per_pix = header[CDELT1]
    dec_offset = header['CRVAL2'] - header['CRPIX2'] * header[CDELT2]
    dec_per_pix = header[CDELT2]

    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(['%.3F' % (item * ra_per_pix + ra_offset) for item in ax.get_xticks()], fontsize=fontsize)
    ax.set_yticklabels(['%.3F' % (item * dec_per_pix + dec_offset) for item in ax.get_yticks()], fontsize=fontsize)

    ax.set_xlabel('RA (deg)', fontsize=fontsize)
    ax.set_ylabel('Dec (deg)', fontsize=fontsize)

    return fig, bg_img_hdu

# -------------------------------------------------------------------------------------------------------
def overplot_skyregion_from_fits(filename, bg_img_hdu, ax, ext=0, color='green', label=None):
    '''
    Overplots the sky region corresponding to a given fits filename, on a given background hdu, on a given axis handle
    Returns axis handle
    '''
    target_wcs = pywcs.WCS(bg_img_hdu[0].header)
    filename = Path(filename)
    hdul = fits.open(filename)
    source_wcs = pywcs.WCS(hdul[ext].header)

    region_dir = filename.parent / 'regions'
    region_dir.mkdir(parents=True, exist_ok=True)
    region_file = region_dir / Path(filename.stem + '.reg')
    source_wcs.footprint_to_file(region_file, color=color, width=1)

    sky_region = Regions.read(region_file, format='ds9')[0]
    color = sky_region.visual['facecolor']
    color = color.replace(',', '') # remove any weird commas
    pixel_region = sky_region.to_pixel(target_wcs)
    pixel_region.plot(ax=ax, color=color)
    ax.text(pixel_region.vertices.x[0], pixel_region.vertices.y[0], filename.stem if label is None else label, color=color)

    return ax

# -------------------------------------------------------------------------------------------------------
def overplot_data_from_fits(filename, bg_img_hdu, ax, ext=0, cmap='Greens'):
    '''
    Overplots the data from a given fits filename, on a given background hdu, on a given axis handle
    Returns axis handle
    '''
    foreground_hdu = fits.open(filename)[ext]
    foreground_data, footprint = reproject_interp(foreground_hdu, bg_img_hdu[0].header)
    data_to_plot = np.log10(foreground_data)
    ax.imshow(data_to_plot, cmap=cmap)

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------directory structures and file names-------------
    args.bg_image_dir = Path(args.bg_image_dir)
    args.reg_files_dir = Path(args.reg_files_dir)
    args.fg_image_dir = Path(args.fg_image_dir)

    bg_filename = args.bg_image_dir / args.bg_file

    # ------plotting the background------------
    fig, bg_img_hdu = plot_background(bg_filename, cmap='Greys', fontsize=args.fontsize)

    # --------getting region file names-------------
    reg_filenames = list(args.reg_files_dir.glob('*PASSAGE*.reg'))
    if not args.only_passage_regions:
        if args.fg_file is not None and 'miri' in args.fg_file: reg_filenames += list(args.reg_files_dir.glob('*COSMOS*MIRI*.reg'))
        elif args.fg_file is not None and 'nircam' in args.fg_file: reg_filenames += list(args.reg_files_dir.glob('*COSMOS*NIRCam*.reg'))
        else: reg_filenames += list(args.reg_files_dir.glob('*COSMOS*.reg'))

    # ------plotting the region files------------
    if len(reg_filenames) > 0: fig, sky_regions = overplot_regions(reg_filenames, bg_img_hdu, fig, fontsize=args.fontsize)

    # ------overplotting data as foreground-----------------
    if args.fg_file is not None:
        fg_files = glob.glob(str(args.fg_image_dir) + '/' + str(args.fg_file))
        if len(fg_files) == 0: print(f'No foreground file found with the pattern {str(args.fg_image_dir) + "/" + str(args.fg_file)}')
        for index, fg_filename in enumerate(fg_files):
            thisfile = Path(fg_filename).stem
            print(f'Overplotting foreground {thisfile} which is {index + 1} of {len(fg_files)}..')
            ax = overplot_skyregion_from_fits(fg_filename, bg_img_hdu, fig.axes[0], ext=0, color='magenta', label=thisfile) # overplots just the outer region of the given fits file
            ax = overplot_data_from_fits(fg_filename, bg_img_hdu, fig.axes[0], ext=0, cmap='Greens' if 'nircam' in args.fg_file else 'Reds') # overplots actual data fo the given fits file

    # -------overplotting specific coordinates from a dataframe----------------
    if args.plot_objects_from_file is not None:
        df = pd.read_csv(args.plot_objects_from_file)
        radius = 10 # arcsec
        print(f'Trying to overplot {radius}" circles around {len(df)} objects..')
        fig = plot_skycoord_from_df(df, fig, bg_img_hdu, color='red', alpha=1, size=radius)

    figname = args.bg_image_dir.parent / f'COSMOS_with_footprints.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
