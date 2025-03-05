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
             run plot_footprints.py --fg_file mosaic_miri_f770w_COSMOS-Web_60mas_A*_v0_5_i2d.fits
             run plot_footprints.py --bg_image_dir /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_data/COSMOS/imaging/ --bg_file mosaic_nircam_f444w_COSMOS-Web_60mas_B7_v0_8_sci.fits --only_passage_regions --do_field Par028 --drv 0.5 --plot_conditions SNR --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --use_only_good
             run plot_footprints.py --fg_file mosaic_nircam_f444w_COSMOS-Web_60mas_B7_v0_8_sci.fits  --only_passage_regions --do_field Par028 --drv 0.5 --plot_conditions SNR --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --use_only_good --plot_cosmoswebb --plot_cosmos2020
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
                            fig = plot_skycoord_from_df(df_contained, ax, color='red', alpha=1, size=2, fontsize=args.fontsize) # plotting contained objects in a different color

                            if args.plot_zcosmos_objects: survey = 'zCOSMOS'
                            elif args.plot_cosmos2020_objects: survey = 'COSMOS2020'
                            print(f'Region {label_text} contains {n_sources} {survey} sources') #
                            ax.text(label_pixcoord_x + 40, label_pixcoord_y, f'({n_sources})', c=color, ha='left', va='top', fontsize=args.fontsize/1.5)


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

    # use this option for a slower plot, but with scale-able scatter points
    for index, row in df.iterrows():
        circle = plt.Circle(xy=(row['ra'], row['dec']), radius=size, color=color, alpha=alpha, transform=ax.get_transform('fk5'))
        ax.add_patch(circle)
        if len(df) < 20 and 'objid' in df: ax.text(row['ra'], row['dec'], row['objid'], color=color, fontsize=fontsize, transform=ax.get_transform('fk5'))

    '''
    # use this option for a quicker plot, but not with scale-able scatter points
    ra_coords, dec_coords = sky_coords.to_pixel(wcs_header)
    ax.scatter(df['ra'], df['dec'], color=color, s=size, alpha=alpha, zorder=-1, transform=ax.get_transform('fk5'))
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
    sci_ext = 0 if 'COSMOS-Web' in str(filename) else 0
    data = np.log10(bg_img_hdu[sci_ext].data)
    header = bg_img_hdu[sci_ext].header

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=pywcs.WCS(header))
    fig.subplots_adjust(right=0.99, top=0.95, bottom=0.1, left=0.01)

    ax.imshow(data, origin='lower', cmap=cmap, vmin=-4, vmax=-3) # full, relevant range of data, for COSMOS field, is roughly (-6,-1)

    ax.coords['ra'].set_major_formatter('d.ddd')
    ax.coords['dec'].set_major_formatter('d.ddd')
    ax.coords['ra'].set_axislabel('RA', fontsize=args.fontsize)
    ax.coords['dec'].set_axislabel('Dec', fontsize=args.fontsize)

    ax.coords['ra'].set_ticks(number=5)
    ax.coords['dec'].set_ticks(number=5)

    ax.coords['ra'].set_ticklabel(size=args.fontsize)
    ax.coords['dec'].set_ticklabel(size=args.fontsize)

    ax.coords.grid(color='k', alpha=0.5, linestyle='solid')

    #ax.coords['ra'].set_limit(150.15, 150.05) ##
    #ax.coords['dec'].set_limit(2.36, 2.45) ##

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
    elif args.fg_file is not None and 'miri' in args.fg_file: reg_filenames = list(reg_files_dir.glob('*%s*MIRI*.reg' %args.bg))
    elif args.fg_file is not None and 'nircam' in args.fg_file: reg_filenames = list(reg_files_dir.glob('*%s*NIRCam*.reg' %args.bg))
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
            if args.plot_fg_data: ax = overplot_data_from_fits(fg_filename, bg_img_hdu, fig.axes[0], ext=0, cmap='Greens' if 'nircam' in args.fg_file else 'Reds')

    # -------reading in dataframe produced by get_field_stats.py or by compute_stellar_masses.py----------------
    if args.do_field is not None:
        if args.plot_conditions == 'all_match':
            args.field = f'Par{int(args.do_field.split("Par")[1]):03d}'
            args.field_set_plot_conditions_text = f'{args.field}_{args.drv}_allmatch'
            df_infilename = args.output_dir / args.field / f'{args.field}_all_diag_results_withSED_{args.run}.csv'
        else:
            args.field = f'Par{int(args.do_field.split("Par")[1]):03d}'
            plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
            plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW',
                                                                                                         f'EW>{args.EW_thresh}').replace(
                'a_image', f'a>{args.a_thresh}')
            args.field_set_plot_conditions_text = f'{args.field}_{args.drv}_venn_{plot_conditions_text}'
            if len(args.run) > 0 and args.run[0] != '_': args.run = '_' + args.run
            df_infilename = args.output_dir / 'catalogs' / f'{args.field_set_plot_conditions_text}_df_withSED{args.run}.csv'

        if os.path.exists(df_infilename):
            print(f'Reading in main df from {df_infilename}')
        else:
            print(f'Could not find {df_infilename},')
            df_infilename = Path(str(df_infilename)[:str(df_infilename).find('withSED') - 1] + '.txt')
            if os.path.exists(df_infilename):
                print(f'Loading pre-SED df from {df_infilename} instead')
            else:
                print(f'Could not find {df_infilename},')
                df_infilename = Path(str(df_infilename).replace(f'_{args.drv}', ''))
                if os.path.exists(df_infilename):
                    print(f'Loading df from {df_infilename} instead')
                else:
                    sys.exit(f'{df_infilename} does not exist')

        df = pd.read_csv(df_infilename)
        if args.use_only_good and args.drv == 'v0.5' and set(args.plot_conditions) == set(['SNR']) and set(args.line_list) == set(['OIII', 'Ha', 'OII', 'Hb', 'SII']):
            df = df[df['objid'].isin([1303,1934,2734,2867,300,2903])].reset_index(drop=True) # only choosing the pre-determined good galaxies
            print(f'Using only the pre-determined good galaxies, and there are {len(df)} of them..')

        fig = plot_skycoord_from_df(df, fig.axes[0], color='red', alpha=1, size=1, fontsize=args.fontsize)

    # ------getting sky region of interest-----------------
    passage_reg_filename = list(reg_files_dir.glob('*PASSAGE*.reg'))[0]
    sky_regions = Regions.read(passage_reg_filename, format='ds9')
    Par28_sky_region = [item for item in sky_regions if 'text' in item.meta and item.meta['text'] == 28][0]

    # ------plotting the circles and objects from a list of catalog---------
    if args.plot_zcosmos_objects:
        df = read_zCOSMOS_catalog(args=args)
        df = df[Par28_sky_region.contains(SkyCoord(df['ra'], df['dec'], unit='deg'), pywcs.WCS(bg_img_hdu[0].header))]
        fig = plot_skycoord_from_df(df, fig.axes[0], color='aqua', alpha=0.3, size=1, fontsize=args.fontsize)
    elif args.plot_cosmos2020_objects:
        df = read_COSMOS2020_catalog(args=args, filename=args.input_dir / 'COSMOS/COSMOS2020_CLASSIC_R1_v2.2_p3_subsetcolumns_nofluxes.fits')
        df = df[Par28_sky_region.contains(SkyCoord(df['ra'], df['dec'], unit='deg'), pywcs.WCS(bg_img_hdu[0].header))]
        fig = plot_skycoord_from_df(df, fig.axes[0], color='pink', alpha=0.8, size=1, fontsize=args.fontsize)
    elif args.plot_cosmoswebb_objects:
        df = read_COSMOSWebb_catalog(args=args, filename=args.input_dir / 'COSMOS/COSMOS_Web_for_Ayan_Dec24_nofluxes.fits')
        df = df[Par28_sky_region.contains(SkyCoord(df['ra'], df['dec'], unit='deg'), pywcs.WCS(bg_img_hdu[0].header))]
        fig = plot_skycoord_from_df(df, fig.axes[0], color='lightblue', alpha=0.8, size=1, fontsize=args.fontsize)

    # ------saving figure---------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    zCOSMOS_text = '_with_zCOSMOS' if args.plot_zcosmos_objects else '_with_COSMOS2020' if args.plot_cosmos2020_objects else '_with_COSMOSWebb' if args.plot_cosmoswebb_objects else ''
    figname = args.input_dir / 'footprints' / f'{args.bg}_with_footprints{zCOSMOS_text}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
