'''
    Filename: plot_footprints.py
    Notes: Plots all PASSAGE field footprints and overlays with existing large survey footprints from MAST
           This script is still W.I.P and does not yet work fully (for now, the workaround is plotting the footprints manually using DS9)
    Author : Ayan
    Created: 18-06-24
    Last modified: 18-06-24
    Example: run plot_footprints.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field COSMOS
             run plot_footprints.py --field COSMOS --plot_zcosmos
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def plot_footprints(region_files, bg_img_hdu, fig, args, table=None):
    '''
    Plots the footprint/s for a given region file and existing figure and header of the background file plotted on the figure
    Returns fig handle
    '''
    col_arr = ['blue', 'green', 'yellow', 'cyan']
    pix_offset_forlabels = 20

    ax = fig.gca()
    wcs_header = wcs.WCS(bg_img_hdu[0].header)
    region_files = np.atleast_1d(region_files)

    for index, region_file in enumerate(region_files):
        print(f'Reading in region file {region_file}..')

        sky_regions = Regions.read(region_file, format='ds9')

        if 'color' in sky_regions[0].visual: color = sky_regions[0].visual['color']
        else: color = col_arr[index]


        for sky_region in sky_regions:
            if type(sky_region) == regions.shapes.text.TextSkyRegion: # label if it is text
                label = sky_region.text
                ax.text(ax.get_xlim()[0] * 1.1, ax.get_ylim()[1] * 0.98 - index * 0.05 * np.diff(ax.get_ylim())[0], label, c=color, ha='left', va='top', fontsize=args.fontsize)
            else: # otherwise plot it
                # plotting the region
                pixel_region = sky_region.to_pixel(wcs_header)
                pixel_region.plot(ax=ax, lw=2, color=color)

                # labeling the region
                if type(sky_region) == regions.shapes.rectangle.RectangleSkyRegion:
                    label_pixcoord_x = pixel_region.center.xy[0] + pixel_region.width/2 + pix_offset_forlabels
                    label_pixcoord_y = pixel_region.center.xy[1] + pixel_region.height/2 + pix_offset_forlabels
                    label_text = f'P{pixel_region.meta["text"]}'
                    ax.text(label_pixcoord_x, label_pixcoord_y, label_text, c=color, ha='left', va='top', fontsize=args.fontsize/1.5)

                    # detecting sources from <table> that lie within this region, if <table> provided
                    if table is not None:
                        n_sources = np.sum(sky_region.contains(SkyCoord(table['RAJ2000'], table['DEJ2000'], unit='deg'), wcs_header))
                        if n_sources > 0:
                            print(f'Region {label_text} contains {n_sources} zCOSMOS sources') #
                            ax.text(label_pixcoord_x + 90, label_pixcoord_y, f'({n_sources})', c=color, ha='left', va='top', fontsize=args.fontsize/1.5)


    return fig, sky_regions

# -------------------------------------------------------------------------------------------------------
def plot_zCOSMOS(fig, bg_img_hdu, color='aqua'):
    '''
    Plots the location of all zCOSMOS galaxies (with spectra) given on an existing background image, given the header of the background image
    Returns fig handle
    '''
    zCOSMOS_catalog_file = '/Users/acharyya/Work/astro/passage/passage_data/zCOSMOS-DR3/zCOSMOS_VIMOS_BRIGHT_DR3_CATALOGUE.fits'

    data = fits.open(zCOSMOS_catalog_file)
    table = Table(data[1].data)
    sky_coords = SkyCoord(table['RAJ2000'], table['DEJ2000'], unit='deg')

    wcs_header = wcs.WCS(bg_img_hdu[0].header)
    ax=fig.gca()
    #'''
    # use this option for a quicker plot, but not with scale-able scatter points
    ra_coords, dec_coords = sky_coords.to_pixel(wcs_header)
    ax.scatter(ra_coords, dec_coords, color=color, s=1, alpha=0.3, zorder=-1)
    '''
    # use this option for a slower plot, but with scale-able scatter points
    pix_coords = np.transpose(sky_coords.to_pixel(wcs_header))
    for index in range(len(table)):
        circle = plt.Circle(xy=pix_coords[index], radius=1, color=color, alpha=0.5)
        ax.add_patch(circle)
    '''
    return fig, table

# -------------------------------------------------------------------------------------------------------
def plot_background(filename, args):
    '''
    Plots the background image for a given input filename
    Returns fig handle
    '''
    cmap = 'Greys'
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(right=0.99, top=0.95, bottom=0.1, left=0.01)

    print(f'Reading in background image file {filename}')
    bg_img_hdu = fits.open(filename)
    data = np.log10(bg_img_hdu[0].data)
    header = bg_img_hdu[0].header

    ax.imshow(data, origin='lower', cmap=cmap, vmin=-4, vmax=-3) # full, relevant range of data, for COSMOS field, is roughly (-6,-1)

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
    if args.plot_zcosmos: fig, table = plot_zCOSMOS(fig, bg_img_hdu)
    fig, sky_regions = plot_footprints(reg_filenames, bg_img_hdu, fig, args, table=table if 'table' in locals() else None)

    # ------saving figure---------
    zCOSMOS_text = '_with_zCOSMOS' if args.plot_zcosmos else ''
    figname = args.input_dir / 'footprints' / f'{args.field}_with_footprints{zCOSMOS_text}.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
