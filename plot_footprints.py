'''
    Filename: plot_footprints.py
    Notes: Plots all PASSAGE field footprints and overlays with existing large survey footprints from MAST
           This script is still W.I.P and does not yet work fully (for now, the workaround is plotting the footprints manually using DS9)
    Author : Ayan
    Created: 18-06-24
    Last modified: 18-06-24
    Example: run plot_footprints.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field COSMOS
             run plot_footprints.py --field COSMOS
'''

from header import *
from util import *
from make_region_files import get_passage_footprints

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def overlay_footprint(img_hdul_orig, img_hdul_new):
    '''
    Borrowed from Peter Watson
    Not yet used in this script
    '''
    xpix, ypix = img_hdul_orig[1].data.shape
    xx = np.asarray([0, 0, xpix, xpix,0])
    yy = np.asarray([0, ypix, ypix, 0, 0])
    orig_celestial = wcs.WCS(img_hdul_orig[0].header).celestial
    new_celestial = wcs.WCS(img_hdul_new[0].header).celestial
    x_p, y_p = wcs.utils.pixel_to_pixel(orig_celestial, new_celestial, yy, xx)

    # extent_in_sky_coords = wcs.WCS(img_hdul_orig[0].header).calc_footprint

    return x_p, y_p

# -------------------------------------------------------------------------------------------------------
def plot_footprints(region_files, bg_img_hdu, fig, args=None):
    '''
    Plots the footprint/s for a given region file and existing figure and header of the background file plotted on the figure
    Returns fig handle
    '''
    col_arr = ['blue', 'green', 'yellow', 'cyan']
    region_files = np.atleast_1d(region_files)

    for index, region_file in enumerate(region_files):
        print(f'Reading in region file {region_file}..')
        all_regions = Regions.read(region_file, format='ds9')
        wcs_header = wcs.WCS(bg_img_hdu[0].header)

        if type(all_regions[0]) == regions.shapes.text.TextSkyRegion: label = all_regions[0].text
        else: label = os.path.splitext(os.path.split(region_file)[1])[0]

        ax = fig.gca()
        ax.text(ax.get_xlim()[0] * 1.01, ax.get_ylim()[1] * 0.98 - index * 0.2, label, c=col_arr[index], ha='left', va='top', fontsize=args.fontsize)
        for region in all_regions:
            try: region.to_pixel(wcs_header).plot(ax=ax, lw=2, color=col_arr[index])
            except AttributeError: pass

    return fig, all_regions

# -------------------------------------------------------------------------------------------------------
def plot_background(filename, args):
    '''
    Plots the background image for a given input filename
    Returns fig handle
    '''
    cmap = 'Grays'
    fig, ax = plt.subplots(figsize=(8, 6))

    print(f'Reading in background image file {filename}')
    bg_img_hdu = fits.open(filename)
    data = np.log10(bg_img_hdu[0].data)
    header = bg_img_hdu[0].header

    ax.imshow(data, origin='lower', cmap=cmap)

    ra_offset = header['CRVAL1'] - header['CRPIX1'] * header['CDELT1']
    ra_per_pix = header['CDELT1']
    dec_offset = header['CRVAL2'] - header['CRPIX2'] * header['CDELT2']
    dec_per_pix = header['CDELT2']

    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(['%.1F' % (item * ra_per_pix + ra_offset) for item in ax.get_xticks()], fontsize=args.fontsize)
    ax.set_yticklabels(['%.1F' % (item * dec_per_pix + dec_offset) for item in ax.get_yticks()], fontsize=args.fontsize)

    ax.set_xlabel('RA (deg)', fontsize=args.fontsize)
    ax.set_ylabel('Dec (deg)', fontsize=args.fontsize)

    return fig, bg_img_hdu

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------directory structures and file names-------------
    bg_image_dir = args.input_dir / 'footprints/survey_images'
    reg_files_dir = args.input_dir / 'footprints/region_files'

    if 'Par' in args.field: args.field = 'COSMOS' # defaults to COSMOS field
    if args.bg_file is None:
        if args.field == 'COSMOS': bg_filename = bg_image_dir / 'COSMOS-HST-ACS_mosaic_Shrink100.fits' # if field is COSMOS, by default use this background image
        else: bg_filename = list(bg_image_dir.glob('*%s*.fits' %args.field))[0] # for other fields, look for available background image files
    else:
        bg_filename = bg_image_dir / args.bg_file

    reg_filenames = list(reg_files_dir.glob('*%s*.reg' %args.field))
    reg_filenames += list(reg_files_dir.glob('*PASSAGE*.reg'))
    if len(reg_filenames) == 0: sys.exit(f'No {args.field} reg file in {reg_files_dir}')

    # ------plotting the background------------
    fig, bg_img_hdu = plot_background(bg_filename, args)

    # ------plotting the footprints---------
    #reg_filenames = ['/Users/acharyya/Work/astro/passage/passage_data/footprints/region_files/COSMOS-Web_NIRCam.reg'] ## only for debugging
    fig, all_regions = plot_footprints(reg_filenames, bg_img_hdu, fig, args)

    # ------saving figure---------
    figname = args.input_dir / 'footprints' / f'{args.field}_with_footprints.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
