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
def plot_skycoord_from_df(df, ax, color='aqua', alpha=0.3, size=1, fontsize=15):
    '''
    Plots the location of all RA/DEC from a given dataframe on an existing background image, using a circle of radius=size in arcseconds,
    given the header of the background image
    Returns fig handle
    '''
    print(f'\nTrying to overplot {size}" circles around {len(df)} objects..')
    size = size / 3600. #conveting from arcseconds to degrees

    for index, row in df.iterrows():
        circle = plt.Circle(xy=(row['ra'], row['dec']), radius=size, color=color, alpha=alpha, transform=ax.get_transform('fk5'))
        ax.add_patch(circle)
        if len(df) < 20 and 'objid' in df: ax.text(row['ra'], row['dec'], row['objid'], color=color, fontsize=fontsize, transform=ax.get_transform('fk5'))

    return fig

# -------------------------------------------------------------------------------------------------------
def plot_background(filename, cmap='Greys', fontsize=15):
    '''
    Plots the background image for a given input filename
    Returns fig handle, and the image HDU
    '''

    print(f'Reading in background image file {filename}')
    bg_img_hdu = fits.open(filename)
    sci_ext = 0 if 'COSMOS-Web' in str(filename) else 0
    data = np.log10(bg_img_hdu[sci_ext].data)
    header = bg_img_hdu[sci_ext].header

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=pywcs.WCS(header))
    fig.subplots_adjust(right=0.99, top=0.95, bottom=0.1, left=0.01)

    ax.imshow(data, origin='lower', cmap=cmap, vmin=-4, vmax=-3) # full, relevant range of data, for COSMOS field, is roughly (-6,-1)

    ax.coords['ra'].set_major_formatter('d.ddd')
    ax.coords['dec'].set_major_formatter('d.ddd')
    ax.coords['ra'].set_axislabel('RA', fontsize=fontsize)
    ax.coords['dec'].set_axislabel('Dec', fontsize=fontsize)

    ax.coords['ra'].set_ticks(number=5)
    ax.coords['dec'].set_ticks(number=5)

    ax.coords['ra'].set_ticklabel(size=fontsize)
    ax.coords['dec'].set_ticklabel(size=fontsize)

    ax.coords.grid(color='k', alpha=0.5, linestyle='solid')

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
    reg_filenames += list(args.reg_files_dir.glob('*Cy2*.reg'))
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
        fig = plot_skycoord_from_df(df, fig.axes[0], color='red', alpha=1, size=2, fontsize=args.fontsize)

    figname = args.bg_image_dir.parent / f'COSMOS_with_footprints.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
