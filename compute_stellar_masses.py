'''
    Filename: compute_stellar_masses.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020
    Author : Ayan
    Created: 19-08-24
    Example: run compute_stellar_masses.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run compute_stellar_masses.py --plot_transmission --plot_SED
             run compute_stellar_masses.py --plot_cutouts --plot_all --arcsec_limit 1 --only_seg
             run compute_stellar_masses.py --plot_cutouts --plot_cutout_errors
             run compute_stellar_masses.py --plot_cutouts
'''
from header import *
from util import *
import random

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def read_filter_transmission(filter_arr, args, verbose=False):
    '''
    Function to load transmission curves for filters in COSMOS2020 catalog
    Returns dataframe
    '''
    transmission_file_path = args.input_dir / 'COSMOS' / 'transmission_curves'
    df_master = pd.DataFrame()
    filter_arr = [filter[: filter.lower().find('_flux')].replace('SPLASH', 'IRAC') if '_flux' in filter.lower() else filter.replace('SPLASH', 'IRAC') for filter in filter_arr]
    filter_arr = list(dict.fromkeys(filter_arr)) # to remove duplicate filter entries

    for index, filter in enumerate(filter_arr):
        if verbose: print(f'Doing filter {filter} which is {index + 1} out of {len(filter_arr)}..')
        transmission_filename = transmission_file_path / f'{filter}.txt'

        if os.path.exists(transmission_filename):
            if 'NIRISS' in filter:
                df = pd.read_table(transmission_filename, comment='#', names=['wave', 'trans', 'dummy'], header=0, usecols=['wave', 'trans'], delim_whitespace=True)
                df['wave'] = df['wave'] * 1e4 # converting from microns to Angstroms
                df = df[df['trans'] > 1e-3]
            else:
                df = pd.read_table(transmission_filename, comment='#', names=['wave', 'trans'], delim_whitespace=True)
        elif os.path.exists(str(transmission_filename).replace('SC_IA', 'SC_IB')):
            df = pd.read_table(str(transmission_filename).replace('SC_IA', 'SC_IB'), comment='#', names=['wave', 'trans'], delim_whitespace=True)
        else:
            print(f'Transmission file for {filter} does not exist.')
            continue

        df['filter'] = filter
        df_master = pd.concat([df_master, df])

    return df_master

# -------------------------------------------------------------------------------------------------------
def plot_filter_transmission(df_master, args, x_scale='linear', color_by_wave=True, plot_fwhm=True):
    '''
    Function to plot transmission curves for filters in COSMOS2020 catalog
    Saves plot as png figure and
    Returns figure handle
    '''

    fig, ax = plt.subplots(figsize=(8, 6))
    jet = plt.get_cmap('rainbow')

    wave_min, wave_max = np.log10(df_master['wave'].min()), np.log10(df_master['wave'].max())
    cnorm = mplcolors.Normalize(vmin=wave_min, vmax=wave_max)
    scalarMap = mpl_cm.ScalarMappable(norm=cnorm, cmap=jet)
    filter_arr = pd.unique(df_master['filter'])

    for filter in filter_arr:
        df = df_master[df_master['filter'] == filter]
        if color_by_wave: color = scalarMap.to_rgba(np.log10(df['wave'].mean()))
        else: color = scalarMap.to_rgba(random.uniform(wave_min, wave_max))
        ax.plot(df['wave'], df['trans'], color=color, label=filter)

        if plot_fwhm:
            trans_max = df['trans'].max()
            wave_llim = df['wave'].iloc[np.where(df['trans'] >= 0.5 * trans_max)[0][0]]
            wave_ulim = df['wave'].iloc[np.where(df['trans'] >= 0.5 * trans_max)[0][-1]]
            ax.plot([wave_llim, wave_ulim], [0.5 * trans_max, 0.5 * trans_max], c=color, lw=1, ls='dashed')

    plt.legend(fontsize = args.fontsize/3)
    ax.set_xlabel(r'Wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel('Transmission', fontsize=args.fontsize)
    ax.set_xscale(x_scale)

    ax.set_xlim(1e3, 1.2e5)
    ax.set_ylim(-0.05, 1.05)

    ax.set_xticklabels(['%.0e' % item for item in ax.get_xticks()], fontsize=args.fontsize)
    ax.set_yticklabels(['%.1f' % item for item in ax.get_yticks()], fontsize=args.fontsize)

    figname = args.input_dir / 'COSMOS' / 'transmission_curves' / 'filter_responses.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# -------------------------------------------------------------------------------------------------------
def plot_SED(df_fluxes, df_trans, args, x_scale='linear'):
    '''
    Function to plot SED based on input dataframe of fluxes and a separate input dataframe with list of filter transmission curves
    Saves plot as png figure and
    Returns figure handle
    '''
    fluxcols = [item for item in df_fluxes.columns if 'flux' in item.lower() and 'fluxerr' not in item.lower()]
    filters, waves_cen, waves_width = [], [], []

    for fluxcol in fluxcols:
        filter = fluxcol[: fluxcol.lower().find('_flux')].replace('SPLASH', 'IRAC') if '_flux' in fluxcol.lower() else fluxcol.replace('SPLASH', 'IRAC')
        thisdf = df_trans[df_trans['filter'] == filter]
        thisdf = thisdf.sort_values(by='wave')

        trans_max = thisdf['trans'].max()
        wave_llim = thisdf['wave'].iloc[np.where(thisdf['trans'] >= 0.5 * trans_max)[0][0]]
        wave_ulim = thisdf['wave'].iloc[np.where(thisdf['trans'] >= 0.5 * trans_max)[0][-1]]
        wave_width = wave_ulim - wave_llim
        wave_cen = np.mean([wave_llim, wave_ulim])

        filters.append(filter)
        waves_cen.append(wave_cen)
        waves_width.append(wave_width)

    waves_cen = np.array(waves_cen)
    waves_width = np.array(waves_width)

    fig, ax = plt.subplots(figsize=(8, 6))
    color_arr = ['rebeccapurple', 'chocolate', 'darkgreen', 'darkblue', 'crimson', 'darkkhaki', 'salmon', 'cornflowerblue']

    for index, row in df_fluxes.iterrows():
        fluxes = np.array([row[thiscol] for thiscol in fluxcols]) # micro Jansky units
        fluxes_err = np.array([row[thiscol.replace('flux', 'fluxerr').replace('FLUX', 'FLUXERR')] for thiscol in fluxcols])

        sorted_indices = np.argsort(waves_cen)
        fluxes_sorted = fluxes[sorted_indices]
        fluxes_err_sorted = fluxes_err[sorted_indices]
        waves_width_sorted = waves_width[sorted_indices]
        waves_cen_sorted = waves_cen[sorted_indices]

        ax.plot(waves_cen_sorted, fluxes_sorted, 'o', markersize=5, color=color_arr[index], label=f'ID #{row["objid"]}')
        ax.errorbar(waves_cen_sorted, fluxes_sorted, xerr=waves_width_sorted, yerr=fluxes_err_sorted, c=color_arr[index], fmt='none', lw=1, alpha=0.4)

    plt.legend(fontsize = args.fontsize/2.)
    ax.set_xlabel(r'Wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel(r'Flux ($\mu$Jy)', fontsize=args.fontsize)
    ax.set_xscale(x_scale)

    ax.set_xlim(1e3, 1.2e5)
    ax.set_ylim(-0.5, 0.8)

    ax.set_xticklabels(['%.0e' % item for item in ax.get_xticks()], fontsize=args.fontsize)
    ax.set_yticklabels(['%.1f' % item for item in ax.get_yticks()], fontsize=args.fontsize)

    figname = args.output_dir / f'{args.intersection_conditions}_SED.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# -------------------------------------------------------------------------------------------------------
def annotate_axis(ax, col_index, row_index, row, filter, n_obj, args):
    '''
    To annotate one individual axis properly
    Returns axis handle
    '''
    if row_index == 0: ax.text(0.9, 0.9, filter, c='k', ha='right', va='top', fontsize=args.fontsize / 2 if len(filter) < 10 else args.fontsize / 5, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
    if col_index == 0:
        ax.set_ylabel('Dec (")', fontsize=args.fontsize / 2)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize / 2)
        ax.text(0.1, 0.1, f'{row["redshift"]:.2f}', c='k', ha='left', va='bottom', fontsize=args.fontsize / 2, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
    else:
        ax.set_yticklabels([])

    if row_index == n_obj - 1:
        ax.set_xlabel('RA (")', fontsize=args.fontsize / 2)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize / 2)
    else:
        ax.set_xticklabels([])

    return ax

# -------------------------------------------------------------------------------------------------------
def extract_header(source_data, source_header):
    '''
    Extracts a smaller, more digestible (for cutout2D) header from radio image header
    Returns extracted header
    Adapted from Alessandro Ignesti's code
    '''
    hdu = fits.PrimaryHDU(source_data)
    for kw in 'CTYPE', 'CRVAL', 'CRPIX', 'CDELT', 'CUNIT':
        for n in [1, 2]:
            if f'{kw}{n}'in source_header: hdu.header.append((f'{kw}{n}', source_header[f'{kw}{n}']))

    for kw in ['BMAJ', 'BMIN', 'BPA']:
        if kw in source_header: hdu.header.append((kw, source_header[kw]))

    return hdu.header

# -------------------------------------------------------------------------------------------------------
def get_error_map(filepath, shape, isweight=False):
    '''
    Function to read in the errormap for a given filename, if file exists, and returns dummy array if file des not exist
    If the provided file is a weightmap instead of erromap, then calculates errormap from the weightmap
    Returns uncertainty map as 2D image
    '''
    if os.path.isfile(filepath):
        data = fits.open(filepath)
        error_map = data[0].data
        if isweight: error_map = error_map ** (-0.5)
    else:
        print(f'Could not find error file {os.path.split(filepath)[-1]}')
        error_map = np.ones(shape)

    return error_map

# -------------------------------------------------------------------------------------------------------
def get_flux_error(filepath, shape):
    '''
    Function to determine the filename of the weightmap for a given fluxmap, and then compute the uncertainty
    Returns uncertainty map as 2D image
    '''
    filedir, filename = os.path.split(filepath)
    filedir = Path(filedir)

    if '-nd-' in filename or '-fd-' in filename:
        error_filename = filename.replace('int', 'skybg')
        error_map = get_error_map(filedir / error_filename, shape)
    elif '27T16_34_32' in filename:
        id = int(filename.split('.')[-2])
        weight_filename = filename.replace(f'{id}', f'{id+1}')
        error_map = get_error_map(filedir / weight_filename, shape, isweight=True)
    elif 'homogenized' in filename or 'vla' in filename.lower():
        error_filename = filename.replace('.fits', '.rms.fits')
        error_map = get_error_map(filedir / error_filename, shape)
    elif 'irac' in filename:
        error_filename = filename.replace('sci', 'unc')
        error_map = get_error_map(filedir / error_filename, shape)
    elif 'drz' in filename:
        error_filename = filename.replace('sci', 'wht')
        error_map = get_error_map(filedir / error_filename, shape, isweight=True)
    elif 'subaru' in filename or 'cfht' in filename:
        error_filename = filename.replace('sci', 'rms')
        error_map = get_error_map(filedir / error_filename, shape)
    elif 'xmm' in filename:
        error_filename = filename.replace('img', 'bg')
        error_map = get_error_map(filedir / error_filename, shape)
    else:
        print(f'No error file specified')
        error_map = np.ones(shape)

    return  error_map

# -------------------------------------------------------------------------------------------------------
def plot_cutouts(df_fluxes, args):
    '''
    Function to plot 2D image cutouts based on input dataframe of ra, dec and existing mosaic images
    Saves plot as png figure and
    Returns figure handle
    '''
    max_columns_per_page = 10
    cmap = 'viridis'
    image_dir = args.input_dir / 'COSMOS' / 'imaging'
    unc_files = ['rms', 'unc', 'skybg', 'wht', '522', '524', '526', '533', '536'] # removing any files that are actually uncertainty/weight maps
    files_to_not_plot = unc_files + ['xmm', '-int'] # removing x-ray and galex because of their extremely poor spatial res

    if args.fontsize == 10: args.fontsize = 15
    cutout_size = 2 * args.arcsec_limit # in arcsec

    fits_images = glob.glob(str(image_dir / '*.fits'))
    if not args.plot_all: fits_images = [os.path.split(item)[-1] for item in fits_images if not np.array([item2.lower() in item.lower() for item2 in files_to_not_plot]).any()]
    else: fits_images = [os.path.split(item)[-1] for item in fits_images if not np.array([item2.lower() in item.lower() for item2 in unc_files]).any()]

    fits_images.sort()
    n_filters = len(fits_images)

    if args.plot_cutout_errors:
        max_filters_per_page = int(max_columns_per_page / 2)
        n_columns = n_filters * 2
    else:
        max_filters_per_page = max_columns_per_page
        n_columns = n_filters

    n_obj = len(df_fluxes)
    n_figs = int(np.ceil(n_columns / max_columns_per_page))

    all_text = 'all' if args.plot_all else 'subset'
    if args.plot_cutout_errors: all_text += '_wunc'
    figname = args.output_dir / f'{args.intersection_conditions}_{all_text}_{cutout_size:.1f}"_cutouts.pdf'
    pdf = PdfPages(figname)

    # ----------getting the seg map----------------
    if args.only_seg:
        field = np.unique(df_fluxes['field'])[0]
        segmentation_file = args.input_dir /f'{field}' / 'Products' / f'{field}_comb_seg.fits'
        seg_data = fits.open(segmentation_file)
        seg_image = seg_data[0].data
        seg_header = seg_data[0].header
        wcs_seg_header = pywcs.WCS(seg_header)

    # -------setting up for plotting the cutouts-------------------
    print(f'\nTotal {n_obj} x {n_filters} x {2 if args.plot_cutout_errors else 1} = {n_obj * n_columns} cutouts to plot, of size {cutout_size}" each; hence splitting in to {n_figs} figures..')
    fig_arr = []
    for fig_index in range(n_figs):
        print(f'\nMaking figure {fig_index + 1} of {n_figs}..')
        these_fits_images = fits_images[fig_index * max_filters_per_page : min((fig_index + 1) * max_filters_per_page, n_filters)] # slicing the fits image array

        fig, axes = plt.subplots(nrows=n_obj, ncols=min(n_columns, max_columns_per_page), figsize=(12, 7))
        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.07, top=0.97, hspace=0.05, wspace=0.05)

        # ------looping over filters-------------
        for filter_index, thisfile in enumerate(these_fits_images):
            col_index = filter_index * 2 if args.plot_cutout_errors else filter_index
            thisfilename = os.path.splitext(thisfile)[0]
            print(f'Reading in file {thisfilename} which is {fig_index * max_filters_per_page + filter_index + 1} of {n_filters}..')
            thisfilename = thisfilename.replace('COSMOS', '').replace('_030mas_077_sci', '').replace('original', '').replace('psf', '').replace('v1', '').replace('v2', '').replace('v3', '').replace('v5', '').replace('_go2_sci_10', '').replace('img', '').replace('mosaic_Shrink10', '').replace('vla', '').replace('lg_sin_10', '').replace('msmf', '').replace(df_fluxes['field'].values[0], '').replace('drz_sci', '').replace('.', '').replace('_', '').replace('-', '')

            data = fits.open(image_dir / thisfile)
            image = data[0].data
            header = data[0].header
            if args.plot_cutout_errors: image_error = get_flux_error(image_dir / thisfile, np.shape(image))

            if 'CTYPE3' in header: # for radio images
                print(f'Modifying header because {thisfilename} is in radio data format..')
                if len(np.shape(image)) > 2: image = image[0][0]
                if args.plot_cutout_errors and len(np.shape(image_error)) > 2: image_error = image_error[0][0]
                header = extract_header(image, header)

            wcs_header = pywcs.WCS(header)
            filter = header['FILTER'] if 'FILTER' in header else thisfilename
            if filter == 'CLEAR' or  'Thin' in filter: filter = thisfilename

            # ------looping over objects-------------
            for row_index, row in df_fluxes.iterrows():
                coord = SkyCoord(row['ra'],row['dec'], unit = 'deg')
                cutout = Cutout2D(image, coord, cutout_size * u.arcsec, wcs = wcs_header)
                if args.plot_cutout_errors: cutout_error = Cutout2D(image_error, coord, cutout_size * u.arcsec, wcs = wcs_header)
                cutout_header = cutout.wcs.to_header()

                # ------for proper rebinning and applying segmentation map on cutout-------------
                if args.only_seg:
                    source_header = seg_header.copy()
                    seg_cutout = Cutout2D(seg_image, coord, cutout_size * u.arcsec, wcs = wcs_seg_header)
                    seg_cutout_header = seg_cutout.wcs.to_header()
                    source_header.update(seg_cutout_header)  # important step to update the header of the full mosaic to that of just the cut out

                    seg_cutout_data_hdu = fits.ImageHDU(seg_cutout.data, header=source_header)
                    seg_cutout_data_rebinned, _ = reproject_interp(seg_cutout_data_hdu, cutout_header)

                # ------now plotting the cutout flux-------------
                ax = np.atleast_1d(axes[row_index])[col_index]
                ax.imshow(np.log10(cutout.data), origin='lower', extent=args.extent, cmap=cmap)#, vmin=-8, vmax=1)
                #print(f'Deb280: min = {np.nanmin(np.log10(cutout.data))}, max={np.nanmax(np.log10(cutout.data))}') ##
                ax.scatter(0, 0, marker='x', c='r', s=30)
                ax = annotate_axis(ax, col_index, row_index, row, filter, n_obj, args)
                if args.only_seg: ax.contour(seg_cutout_data_rebinned != row['objid'], levels=0, colors='k', extent=args.extent, linewidths=0.5)

                # ------now plotting the cutout error-------------
                if args.plot_cutout_errors:
                    ax = np.atleast_1d(axes[row_index])[col_index + 1]
                    ax.imshow(np.log10(cutout_error.data), origin='lower', extent=args.extent,cmap=cmap)
                    ax.scatter(0, 0, marker='x', c='r', s=30)
                    ax = annotate_axis(ax, col_index + 1, row_index, row, filter + '_U', n_obj, args)
                    if args.only_seg: ax.contour(seg_cutout_data_rebinned != row['objid'], levels=0, colors='k', extent=args.extent, linewidths=0.5)

        # --------hiding excess axes frames--------------
        cols_used_up = len(these_fits_images) * 2 if args.plot_cutout_errors else len(these_fits_images)
        spare_columns = min(n_columns, max_columns_per_page) - cols_used_up
        for col_index in range(spare_columns):
            for row_index in range(len(df_fluxes)):
                ax = axes[row_index][col_index + cols_used_up]
                ax.set_visible(False)

        pdf.savefig(fig)
        fig_arr.append(fig)

    pdf.close()
    print(f'Saved {n_figs} figures in {figname}')
    plt.show(block=False)

    return fig_arr

# -------------------------------------------------------------------------------------------------------
def get_fluxcols(args):
    '''
    Function to load or generate the list of filters and correpsonding flux columns in COSMOS2020 catalog
    Returns list of columns, and optionally the full cosmos2020 dataframe
    '''
    filepath = args.input_dir / 'COSMOS' / 'cosmos_fluxcols.npy'

    if os.path.exists(filepath):
        print(f'Reading flux columns from existing {filepath}')
        fluxcols = np.load(filepath)
        df_cosmos = None
    else:
        print(f'{filepath} does not exist, so preparing the list..')
        df_cosmos = read_COSMOS2020_catalog(args=args, filename=args.cosmos2020_filename)

        all_flux_cols = [item for item in df_cosmos.columns if 'FLUX' in item and item != 'FLUX_RADIUS' and 'FLUXERR' not in item]
        filters = [item[:item.find('FLUX')] for item in all_flux_cols]
        fluxcols = [item + 'FLUX_AUTO' if item + 'FLUX_AUTO' in df_cosmos.columns else item + 'FLUX' for item in filters]
        fluxcols = list(dict.fromkeys(fluxcols)) # to remove duplicates
        np.save(filepath, fluxcols)

    return fluxcols, df_cosmos

# -------------------------------------------------------------------------------------------------------
def get_flux_catalog(df_int, args):
    '''
    Function to load or generate the catalog of flux values for galaxies of interest from the COSMOS2020 catalog, and add PASSAGE fluxes too
    Returns the dataframe containing all flux values
    '''
    fluxfilename = args.output_dir / 'passage_cosmos_fluxes.csv'
    if os.path.exists(fluxfilename) and not args.clobber:
        print(f'Reading flux values from existing {fluxfilename}')
        df_fluxes = pd.read_csv(fluxfilename)
    else:
        print(f'{fluxfilename} does not exist, so preparing the flux list..')
        filename = args.input_dir / 'COSMOS' / 'cosmos_fluxes_subset.csv'

        if os.path.exists(filename):
            print(f'Reading cosmos2020 flux values from existing {filename}')
            df_fluxes = pd.read_csv(filename)
        else:
            print(f'{filename} does not exist, so preparing the flux list..')

            # -------reading in flux columns names and df_cosmos-------
            fluxcols, df_cosmos = get_fluxcols(args)
            if df_cosmos is None: df_cosmos = read_COSMOS2020_catalog(args=args, filename=args.cosmos2020_filename)

            # -------making subset of df_cosmos-------
            cosmos_ids = df_int['cosmos_id'].values
            df_cosmos_subset = df_cosmos[df_cosmos['id'].isin(cosmos_ids)]
            df_cosmos_subset = df_cosmos_subset.rename(columns={'id': 'cosmos_id'})

            # -------determining other columns to extract from df_cosmos-------
            lp_cols_suffix = ['med', 'med_min68', 'med_max68', 'best']
            lp_cols = np.ravel([f'lp_{item}_{suffix}' for item in ['mass', 'SFR', 'sSFR'] for suffix in lp_cols_suffix])
            ez_cols_suffix = ['', '_p160', '_p500', '_p840']
            ez_cols = np.ravel([f'ez_{item}{suffix}' for item in ['mass', 'sfr', 'ssfr'] for suffix in ez_cols_suffix])
            flux_and_err_cols = np.ravel([[item, item.replace('FLUX', 'FLUXERR')] for item in fluxcols])
            cols_to_extract = np.hstack((['cosmos_id', 'ra', 'dec', 'ID_COSMOS2015', 'ez_z_phot', 'lp_MK', 'lp_zBEST'], lp_cols, ez_cols, flux_and_err_cols)).tolist()
            df_fluxes = pd.merge(df_int[['field', 'objid', 'redshift', 'cosmos_id']], df_cosmos_subset[cols_to_extract], on='cosmos_id', how='inner')
            df_fluxes = df_fluxes.dropna(axis=1, how='all')

            # -------writing cosmos fluxes df into file-------
            df_fluxes.to_csv(filename, index=None)
            print(f'Written cosmos2020 flux table as {filename}')

        # -------reading in passage photometric catalog-------
        args.field = df_int['field'].values[0]
        product_dir = args.input_dir / args.field / 'Products'
        photcat_file = product_dir / f'{args.field}_photcat.fits'
        df_photcat = Table(fits.open(photcat_file)[1].data).to_pandas()
        df_photcat = df_photcat.rename(columns={'id': 'objid'})
        df_photcat.columns = df_photcat.columns.str.replace('f115w', 'NIRISS_F115W', regex=True)
        df_photcat.columns = df_photcat.columns.str.replace('f150w', 'NIRISS_F150W', regex=True)
        df_photcat.columns = df_photcat.columns.str.replace('f200w', 'NIRISS_F200W', regex=True)

        # -------determining flux and fluxerr columns from passage-------
        passage_filters = ['NIRISS_F115W', 'NIRISS_F150W', 'NIRISS_F200W']
        aper_num = 4 # aper_num = 0, 1, 2, 3, 4, 5, 6 correspond to fluxes measured within apertures of sizes [0.36", 0.5", 0.7", 1", 1.2", 1.5", 3.0"]
        cols_to_extract = []
        for thisfilter in passage_filters:
            fluxcol = f'{thisfilter}_flux_aper_{aper_num:0d}'
            errcol = fluxcol.replace('flux', 'fluxerr')
            cols_to_extract.append([fluxcol, errcol])
        cols_to_extract = np.array(cols_to_extract).flatten()
        cols_to_extract = np.hstack((['objid'], cols_to_extract)).tolist()

        # -------merging photcat with fluxes dataframe-------
        df_fluxes = pd.merge(df_fluxes, df_photcat[cols_to_extract], on='objid', how='inner')

        # -------writing master fluxes df into file-------
        df_fluxes.to_csv(fluxfilename, index=None)
        print(f'Written passage+cosmos2020 flux table as {fluxfilename}')

    return df_fluxes

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.intersection_conditions = 'allpar_venn_EW,mass,PA'
    if args.fontsize == 10: args.fontsize = 15
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)

    args.cosmos2020_filename = args.input_dir / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3.fits'
    df_int_filename = args.output_dir / f'{args.intersection_conditions}_df.txt'
    df_int = pd.read_csv(df_int_filename)
    df_fluxes = get_flux_catalog(df_int, args)

    fluxcols = [item for item in df_fluxes.columns if 'flux' in item.lower() and 'fluxerr' not in item.lower()]
    df_trans = read_filter_transmission(fluxcols, args)
    filters = pd.unique(df_trans['filter'])

    # ----------plotting-----------------
    if args.plot_transmission: fig = plot_filter_transmission(df_trans, args, x_scale='log')
    if args.plot_SED: fig2 = plot_SED(df_fluxes, df_trans, args, x_scale='log')
    if args.plot_cutouts: fig3 = plot_cutouts(df_fluxes, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
