'''
    Filename: plot_footprints.py
    Notes: Plots all PASSAGE field footprints and overlays with existing large survey footprints from MAST
    Author : Ayan
    Created: 18-06-24
    Last modified: 18-06-24
    Example: run plot_footprints.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field COSMOS
             run plot_footprints.py --bg COSMOS --plot_zcosmos
             run plot_footprints.py --bg_image_dir /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_data/COSMOS/imaging_orig/ --bg_file ACS_814_030mas_077_sci.fits
             run plot_footprints.py --fortalk
             run plot_footprints.py --fg_file mosaic_nircam_f115w_COSMOS-Web_60mas*_v0_5_i2d.fits
             run plot_footprints.py --fg_file mosaic_miri_f770w_COSMOS-Web_60mas_A*_v0_5_i2d.fits
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def plot_footprints(region_files, bg_img_hdu, fig, args, df=None):
    '''
    Plots the footprint/s for a given region file and existing figure and header of the background file plotted on the figure
    Returns fig handle
    '''
    col_arr = ['cornflowerblue', 'lightsalmon', 'gold', 'cyan'] if args.fortalk else ['blue', 'green', 'yellow', 'cyan']
    pix_offset_forlabels = 20

    ax = fig.gca()
    wcs_header = pywcs.WCS(bg_img_hdu[0].header)
    region_files = np.atleast_1d(region_files)

    for index, region_file in enumerate(region_files):
        print(f'Reading in region file {region_file}..')

        sky_regions = Regions.read(region_file, format='ds9')

        if 'color' in sky_regions[0].visual and not args.fortalk: color = sky_regions[0].visual['color']
        else: color = col_arr[index]

        for index2, sky_region in enumerate(sky_regions):
            if type(sky_region) == regions.shapes.text.TextSkyRegion: # label if it is text
                label = sky_region.text
                ax.text(ax.get_xlim()[0] + 0.1 * np.diff(ax.get_xlim())[0], ax.get_ylim()[1] * 0.98 - (index + 1) * 0.05 * np.diff(ax.get_ylim())[0], label, c=color, ha='left', va='top', fontsize=args.fontsize, bbox=dict(facecolor='k', edgecolor='black', alpha=0.5) if args.fortalk else None)
            else: # otherwise plot it
                # plotting the region
                pixel_region = sky_region.to_pixel(wcs_header)
                pixel_region.plot(ax=ax, lw=2 if args.fg_file is None else 1, color=color)

                # labeling the region
                if type(sky_region) == regions.shapes.rectangle.RectangleSkyRegion:
                    label_pixcoord_x = pixel_region.center.xy[0] + pixel_region.width/2 + pix_offset_forlabels
                    label_pixcoord_y = pixel_region.center.xy[1] + pixel_region.height/2 + pix_offset_forlabels
                    label_text = f'P{pixel_region.meta["text"]}'
                    ax.text(label_pixcoord_x, label_pixcoord_y, label_text, c=color, ha='left', va='top', fontsize=args.fontsize/1.5, bbox=dict(facecolor='k', edgecolor='black', alpha=0.5) if args.fortalk else None)

                    # detecting sources from <table> that lie within this region, if <table> provided
                    if df is not None:
                        print(f'Doing field {label_text}, which is {index2 + 1} out of {len(sky_regions)}..')  #
                        contained_ids = sky_region.contains(SkyCoord(df['ra'], df['dec'], unit='deg'), wcs_header)
                        n_sources = np.sum(contained_ids)
                        if n_sources > 0:
                            df_contained = df[contained_ids]
                            fig = plot_skycoord_from_df(df_contained, fig, bg_img_hdu, color='red', alpha=1, size=10) # plotting contained objects in a different color

                            if args.plot_zcosmos: survey = 'zCOSMOS'
                            elif args.plot_cosmos2020: survey = 'COSMOS2020'
                            print(f'Region {label_text} contains {n_sources} {survey} sources') #
                            ax.text(label_pixcoord_x + 40, label_pixcoord_y, f'({n_sources})', c=color, ha='left', va='top', fontsize=args.fontsize/1.5)


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
    #'''
    # use this option for a quicker plot, but not with scale-able scatter points
    ra_coords, dec_coords = sky_coords.to_pixel(wcs_header)
    ax.scatter(ra_coords, dec_coords, color=color, s=size, alpha=alpha, zorder=-1)
    '''
    # use this option for a slower plot, but with scale-able scatter points
    pix_coords = np.transpose(sky_coords.to_pixel(wcs_header))
    for index in range(len(table)):
        circle = plt.Circle(xy=pix_coords[index], radius=1, color=color, alpha=alpha)
        ax.add_patch(circle)
    '''
    return fig

# -------------------------------------------------------------------------------------------------------
def plot_background(filename, args):
    '''
    Plots the background image for a given input filename
    Returns fig handle
    '''
    cmap = 'Greys_r' if args.fortalk else 'Greys'

    print(f'Reading in background image file {filename}')
    bg_img_hdu = fits.open(filename)
    sci_ext = 1 if 'COSMOS-Web' in str(filename) else 0
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
    ax.set_xticklabels(['%.3F' % (item * ra_per_pix + ra_offset) for item in ax.get_xticks()], fontsize=args.fontsize)
    ax.set_yticklabels(['%.3F' % (item * dec_per_pix + dec_offset) for item in ax.get_yticks()], fontsize=args.fontsize)

    ax.set_xlabel('RA (deg)', fontsize=args.fontsize)
    ax.set_ylabel('Dec (deg)', fontsize=args.fontsize)

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
    reg_files_dir = args.input_dir / 'footprints/region_files'

    if args.bg_image_dir is None: args.bg_image_dir = args.input_dir / 'footprints/survey_images'
    else: args.bg_image_dir = Path(args.bg_image_dir)

    if args.fg_image_dir is None: args.fg_image_dir = args.input_dir / 'COSMOS/imaging'
    else: args.fg_image_dir = Path(args.fg_image_dir)

    if args.bg_file is None:
        if args.bg == 'COSMOS': bg_filename = args.bg_image_dir / 'COSMOS-HST-ACS_mosaic_Shrink100.fits' # if field is COSMOS, by default use this background image
        else: bg_filename = list(args.bg_image_dir.glob('*%s*.fits' %args.bg))[0] # for other fields, look for available background image files
    else:
        bg_filename = args.bg_image_dir / args.bg_file

    if args.only_passage_regions: reg_filenames = []
    elif 'miri' in args.fg_file: reg_filenames = list(reg_files_dir.glob('*%s*MIRI*.reg' %args.bg))
    elif 'nircam' in args.fg_file: reg_filenames = list(reg_files_dir.glob('*%s*NIRCam*.reg' %args.bg))
    else: reg_filenames = list(reg_files_dir.glob('*%s*.reg' %args.bg))
    reg_filenames += list(reg_files_dir.glob('*PASSAGE*.reg'))

    if len(reg_filenames) == 0: sys.exit(f'No {args.bg} reg file in {reg_files_dir}')

    # ------plotting the background------------
    fig, bg_img_hdu = plot_background(bg_filename, args)
    fig, sky_regions = plot_footprints(reg_filenames, bg_img_hdu, fig, args, df=df if 'df' in locals() else None)

    # ------overplotting COSMOS-Web mosaics-----------------
    if args.fg_file is not None:
        fg_files = glob.glob(str(args.fg_image_dir) + '/' + str(args.fg_file))
        for index, fg_filename in enumerate(fg_files):
            thisfile = Path(fg_filename).stem
            print(f'Overplotting foreground {thisfile} which is {index + 1} of {len(fg_files)}..')
            ax = overplot_skyregion_from_fits(fg_filename, bg_img_hdu, fig.axes[0], ext=0, color='magenta', label=thisfile)
            ax = overplot_data_from_fits(fg_filename, bg_img_hdu, fig.axes[0], ext=0, cmap='Greens' if 'nircam' in args.fg_file else 'Reds')

    # ------plotting the footprints---------
    if args.plot_zcosmos or args.plot_cosmos2020:
        if args.plot_zcosmos: df = read_zCOSMOS_catalog(args=args)
        elif args.plot_cosmos2020: df = read_COSMOS2020_catalog(args=args)
        fig = plot_skycoord_from_df(df, fig, bg_img_hdu, color='aqua', alpha=0.3, size=1)

    # ------saving figure---------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    zCOSMOS_text = '_with_zCOSMOS' if args.plot_zcosmos else '_with_COSMOS2020' if args.plot_cosmos2020 else ''
    figname = args.input_dir / 'footprints' / f'{args.bg}_with_footprints{zCOSMOS_text}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
