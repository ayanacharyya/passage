'''
    Filename: compute_stellar_masses.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020
    Author : Ayan
    Created: 19-08-24
    Example: run compute_stellar_masses.py  --plot_conditions EW,mass,PA --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run compute_stellar_masses.py  --plot_conditions EW,mass,PA --plot_transmission --plot_SED
             run compute_stellar_masses.py --plot_conditions EW,mass,PA --plot_filter_cutouts --plot_all --arcsec_limit 1 --only_seg
             run compute_stellar_masses.py --plot_conditions EW,mass,PA --plot_filter_cutouts --plot_cutout_errors
             run compute_stellar_masses.py --plot_conditions EW,mass,PA --plot_filter_cutouts
             run compute_stellar_masses.py --plot_conditions EW,mass,PA --fit_sed --run narrow_z
             run compute_stellar_masses.py --plot_conditions SNR,mass,F115W,F150W,F200W --plot_niriss_direct --filters F115W,F150W,F200W
             run compute_stellar_masses.py --plot_conditions SNR,mass,F115W,F150W,F200W --fit_sed --run narrow_z
             run compute_stellar_masses.py --plot_conditions SNR,mass,F115W,F150W,F200W --fit_sed --run narrow_z_narrow_mass
'''
from header import *
from util import *
import random

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def read_filter_transmission(filter_dir, filter_arr, args, verbose=False):
    '''
    Function to load transmission curves for filters in COSMOS2020 catalog
    Returns dataframe
    '''
    print(f'Attempting to read all {len(filter_arr)} filter transmission files from {filter_dir}..')
    df_master = pd.DataFrame()
    filter_arr = [filter[: filter.lower().find('_sci')].replace('SPLASH', 'IRAC') if '_sci' in filter.lower() else filter.replace('SPLASH', 'IRAC') for filter in filter_arr]
    filter_arr = list(dict.fromkeys(filter_arr)) # to remove duplicate filter entries

    for index, filter in enumerate(filter_arr):
        if verbose: print(f'Doing filter {filter} which is {index + 1} out of {len(filter_arr)}..')
        transmission_filename = filter_dir / f'{filter}.txt'

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
        filter = fluxcol[: fluxcol.lower().find('_sci')].replace('SPLASH', 'IRAC') if '_sci' in fluxcol.lower() else fluxcol.replace('SPLASH', 'IRAC')
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

    figname = args.output_dir / f'{args.plot_conditions}_SED.png'
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
        if args.plot_filter_cutouts: ax.text(0.1, 0.1, f'{row["redshift"]:.2f}', c='k', ha='left', va='bottom', fontsize=args.fontsize / 2, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
        else: ax.text(0.1, 0.1, f'{row["objid"]}', c='k', ha='left', va='bottom', fontsize=args.fontsize / 2, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
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
def plot_niriss_direct(df_fluxes, args):
    '''
    Function to plot direct PASSAGE images based on input dataframe of ra, dec
    Saves plot as png figure and
    Returns figure handle
    '''
    max_columns_per_page = 10
    max_rows_per_page = 7
    max_cutouts_per_page = max_columns_per_page * max_rows_per_page
    n_obj = len(df_fluxes)
    n_figs = int(np.ceil(n_obj / max_cutouts_per_page))
    cmap = 'viridis'
    image_dir = args.input_dir / args.field / 'Products'

    if args.fontsize == 10: args.fontsize = 15
    cutout_size = 2 * args.arcsec_limit # in arcsec

    for filter in args.filters:
        print(f'\nDoing filter {filter} of {len(args.filters)} filters..')

        field = pd.unique(df_fluxes['field'])[0]
        fits_image = args.input_dir / f'{field}' / 'Products' / f'{field}_228_{filter}-clear_drz_sci.fits'
        data = fits.open(fits_image)
        image = data[0].data
        header = data[0].header
        wcs_header = pywcs.WCS(header)

        figname = args.output_dir / f'{args.plot_conditions}_niriss_direct_{filter}_{cutout_size:.1f}"_cutouts.pdf'
        pdf = PdfPages(figname)

        # ----------getting the seg map----------------
        if args.only_seg:
            segmentation_file = args.input_dir /f'{field}' / 'Products' / f'{field}_comb_seg.fits'
            seg_data = fits.open(segmentation_file)
            seg_image = seg_data[0].data
            seg_header = seg_data[0].header
            wcs_seg_header = pywcs.WCS(seg_header)

        # -------setting up for plotting the cutouts-------------------
        print(f'\nTotal {n_obj} cutouts to plot, of size {cutout_size}" each; hence splitting in to {n_figs} figures..')
        fig_arr = []
        for fig_index in range(n_figs):
            print(f'\nMaking figure {fig_index + 1} of {n_figs}..')
            this_fig_df = df_fluxes[fig_index * max_cutouts_per_page : min((fig_index + 1) * max_cutouts_per_page, n_obj)].reset_index(drop=False) # slicing the df_fluxes dataframe

            fig, axes = plt.subplots(nrows=max_rows_per_page, ncols=max_columns_per_page, figsize=(12, 7))
            fig.subplots_adjust(left=0.05, right=0.98, bottom=0.07, top=0.97, hspace=0.05, wspace=0.05)

            # ------looping over filters-------------
            for index, row in this_fig_df.iterrows():
                print(f'Doing object {fig_index * max_cutouts_per_page + index + 1} of {n_obj}..')
                col_index = index % max_columns_per_page
                row_index = int(index / max_columns_per_page)

                # ------looping over objects-------------
                coord = SkyCoord(row['ra'],row['dec'], unit = 'deg')
                cutout = Cutout2D(image, coord, cutout_size * u.arcsec, wcs = wcs_header)
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
                ax = axes[row_index][col_index]
                p = ax.imshow(np.log10(cutout.data), origin='lower', extent=args.extent, cmap=cmap)
                ax.scatter(0, 0, marker='x', c='r', s=30)
                ax.text(0.9, 0.9, filter, c='k', ha='right', va='top', fontsize=args.fontsize / 2, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                ax.text(0.1, 0.1, f'{row["objid"]}', c='k', ha='left', va='bottom', fontsize=args.fontsize / 2, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                if col_index == 0:
                    ax.set_ylabel('Dec (")', fontsize=args.fontsize / 2)
                    ax.tick_params(axis='y', which='major', labelsize=args.fontsize / 2)
                else:
                    ax.set_yticklabels([])

                if row_index == int(len(this_fig_df) / max_columns_per_page):
                    ax.set_xlabel('RA (")', fontsize=args.fontsize / 2)
                    ax.tick_params(axis='x', which='major', labelsize=args.fontsize / 2)
                else:
                    ax.set_xticklabels([])
                if args.only_seg: ax.contour(seg_cutout_data_rebinned != row['objid'], levels=0, colors='k', extent=args.extent, linewidths=0.5)

            # --------hiding excess axes frames--------------
            for index in range(len(this_fig_df), max_cutouts_per_page):
                col_index = index % max_columns_per_page
                row_index = int(index / max_columns_per_page)
                ax = axes[row_index][col_index]
                ax.set_visible(False)

            pdf.savefig(fig)
            fig_arr.append(fig)

        pdf.close()
        print(f'Saved {n_figs} figures in {figname}')
        plt.show(block=False)

    return fig_arr

# -------------------------------------------------------------------------------------------------------
def plot_filter_cutouts(df_fluxes, args):
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
    figname = args.output_dir / f'{args.plot_conditions}_{all_text}_{cutout_size:.1f}"_cutouts.pdf'
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
                p = ax.imshow(np.log10(cutout.data), origin='lower', extent=args.extent, cmap=cmap)
                ax.scatter(0, 0, marker='x', c='r', s=30)
                ax = annotate_axis(ax, col_index, row_index, row, filter, n_obj, args)
                if args.only_seg: ax.contour(seg_cutout_data_rebinned != row['objid'], levels=0, colors='k', extent=args.extent, linewidths=0.5)

                # ------now plotting the cutout error-------------
                if args.plot_cutout_errors:
                    ax = np.atleast_1d(axes[row_index])[col_index + 1]
                    dummy_error = np.mean(cutout_error.data) == 1 # if error map was filled with only 1s
                    if not dummy_error:
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
def get_flux_catalog(photcat_filename, df_int, args):
    '''
    Function to load or generate the catalog of flux values for galaxies of interest from the COSMOS2020 catalog, and add PASSAGE fluxes too
    Returns the dataframe containing all flux values
    '''
    if os.path.exists(photcat_filename) and not args.clobber:
        print(f'Reading flux values from existing {photcat_filename}')
        df_fluxes = pd.read_csv(photcat_filename)
    else:
        print(f'{photcat_filename} does not exist, so preparing the flux list..')
        filename = args.input_dir / 'COSMOS' / f'{args.plot_conditions}_cosmos_fluxes_subset.csv'

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

        # -------determining flux and fluxerr columns from passage-------
        passage_filters = ['NIRISS_F115W', 'NIRISS_F150W', 'NIRISS_F200W']
        aper_num = 4  # aper_num = 0, 1, 2, 3, 4, 5, 6 correspond to fluxes measured within apertures of sizes [0.36", 0.5", 0.7", 1", 1.2", 1.5", 3.0"]
        cols_to_extract = []
        for thisfilter in passage_filters:
            fluxcol = f'{thisfilter}_flux_aper_{aper_num:0d}'
            errcol = fluxcol.replace('flux', 'fluxerr')
            cols_to_extract.append([fluxcol, errcol])
        cols_to_extract = np.array(cols_to_extract).flatten()

        for thiscol in cols_to_extract: df_fluxes[thiscol] = np.zeros(len(df_fluxes))

        # -------reading in passage photometric catalog-------
        for index, row in df_fluxes.iterrows():
            print(f'Getting NIRISS fluxes for object {index+1} out of {len(df_fluxes)}..')
            field = row['field']
            product_dir = args.input_dir / field / 'Products'
            photcat_file = product_dir / f'{field}_photcat.fits'
            df_photcat = Table(fits.open(photcat_file)[1].data).to_pandas()
            df_photcat.columns = df_photcat.columns.str.replace('f115w', 'NIRISS_F115W', regex=True)
            df_photcat.columns = df_photcat.columns.str.replace('f150w', 'NIRISS_F150W', regex=True)
            df_photcat.columns = df_photcat.columns.str.replace('f200w', 'NIRISS_F200W', regex=True)

            index_in_photcat = df_photcat[df_photcat['id'] == row['objid']].index[0]
            for thiscol in cols_to_extract:
                try: df_fluxes.loc[index, thiscol] = df_photcat[thiscol][index_in_photcat]
                except KeyError: df_fluxes.loc[index, thiscol] = np.nan

        # -------dropping rows with no NIRISS fluxes-------
        df_fluxes = df_fluxes.dropna(axis=0, subset=cols_to_extract).reset_index(drop=True) # drop all objects (rows) that have any NaNs in them

        # -------modiyfying all flux columns to have uniform nomenclature-------
        fluxcols = [item for item in df_fluxes.columns if 'flux' in item.lower() and 'fluxerr' not in item.lower()]
        for fluxcol in fluxcols:
            filter = fluxcol[:fluxcol.lower().find('_flux')]
            errcol = fluxcol.replace('FLUX', 'FLUXERR').replace('flux', 'fluxerr')
            new_flux_col = filter + '_sci'
            new_err_col = filter + '_err'
            df_fluxes = df_fluxes.rename(columns={fluxcol: new_flux_col.replace('SC_IA', 'SC_IB'), errcol: new_err_col.replace('SC_IA', 'SC_IB')})

        # -------writing master fluxes df into file-------
        df_fluxes.to_csv(photcat_filename, index=None)
        print(f'Written passage+cosmos2020 flux table as {photcat_filename}')


    return df_fluxes

# -------------------------------------------------------------------------------------------------------
def load_photom_bagpipes(str_id, phot_cat, id_colname = 'bin_id', zeropoint = 28.9, cat_hdu_index = 0, extra_frac_err = 0.1):
    '''
    Code written by P J Watson
    Load photometry from a catalogue to bagpipes-formatted data.

    The output fluxes and uncertainties are scaled to microJanskys.

    Parameters
    ----------
    str_id : str
        The ID of the object in the photometric catalogue to fit.
    phot_cat : os.PathLike
        The location of the photometric catalogue.
    id_colname : str, optional
        The name of the column containing ``str_id``, by default
        ``"bin_id"``.
    zeropoint : float, optional
        The AB magnitude zeropoint, by default ``28.9``.
    cat_hdu_index : int | str, optional
        The index or name of the HDU containing the photometric catalogue,
        by default ``0``.
    extra_frac_err : float, optional
        An additional fractional error to be added to the photometric
        uncertainties. By default ``extra_frac_err=0.1``, i.e. 10% of the
        measured flux will be added in quadrature to the estimated
        uncertainty.

    Returns
    -------
    ArrayLike
        An Nx2 array containing the fluxes and their associated
        uncertainties in all photometric bands.
    '''

    if not isinstance(phot_cat, Table):
        try: phot_cat = Table.read(phot_cat, hdu=cat_hdu_index)
        except TypeError: phot_cat = Table.read(phot_cat)
    phot_cat[id_colname] = phot_cat[id_colname].astype(str)

    row_idx = (phot_cat[id_colname] == str_id).nonzero()[0][0]
    fluxes = []
    errs = []

    for c in phot_cat.colnames:
        if '_sci' in c.lower(): fluxes.append(phot_cat[c][row_idx])
        elif '_var' in c.lower(): errs.append(np.sqrt(phot_cat[c][row_idx]))
        elif '_err' in c.lower(): errs.append(phot_cat[c][row_idx])

    #if zeropoint == 28.9: flux_scale = 1e-2
    #else: flux_scale = 10 ** ((8.9 - zeropoint) / 2.5 + 6)
    flux_scale = 1

    flux = np.asarray(fluxes) * flux_scale
    flux_err = np.asarray(errs) * flux_scale
    flux_err = np.sqrt(flux_err**2 + (0.1 * flux) ** 2)
    flux = flux.copy()
    flux_err = flux_err.copy()
    bad_values = (~np.isfinite(flux) | (flux <= 0) | ~np.isfinite(flux_err) | (flux_err <= 0))
    flux[bad_values] = 0.0
    flux_err[bad_values] = 1e30
    return np.c_[flux, flux_err]

# -------------------------------------------------------------------------------------------------------
def generate_fit_params(obj_z, z_range = 0.01, num_age_bins = 5, min_age_bin = 30):
    '''
    Code written by P J Watson
    Generate a default set of fit parameters for ``bagpipes``.

    Parameters
    ----------
    obj_z : float | ArrayLike
        The redshift of the object to fit. If ``ArrayLike``, this
        indicates the maximum range of redshifts to fit to.
    z_range : float, optional
        The maximum redshift range to search over, by default 0.01. To fit
        to a single redshift, pass a single value for ``obj_z``, and set
        ``z_range=0.0``. If ``obj_z`` is ``ArrayLike``, this parameter is
        ignored.
    num_age_bins : int, optional
        The number of age bins to fit, each of which will have a constant
        star formation rate following Leja+19. By default ``5`` bins are
        generated.
    min_age_bin : float, optional
        The minimum age of any bin in Myr, by default 30.

    Returns
    -------
    dict
        A dictionary containing the necessary fit parameters for
        ``bagpipes``.
    '''
    fit_params = {}
    if (z_range == 0.0) or (type(obj_z) is ArrayLike): fit_params['redshift'] = obj_z
    else: fit_params['redshift'] = (obj_z - z_range / 2, obj_z + z_range / 2)

    # Set up necessary variables for cosmological calculations.
    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
    age_at_z = cosmo.age(np.nanmax(fit_params['redshift'])).value

    age_bins = np.geomspace(min_age_bin, age_at_z * 1e3, num=num_age_bins)
    age_bins = np.insert(age_bins, 0, 0.0)

    continuity = {
        'massformed': (6.0, 11.0),
        'metallicity': (0.0, 3.0),
        'metallicity_prior_mu': 1.0,
        'metallicity_prior_sigma': 0.5,
        'bin_edges': age_bins.tolist(),
    }

    for i in range(1, len(continuity['bin_edges']) - 1):
        continuity['dsfr' + str(i)] = (-10.0, 10.0)
        continuity['dsfr' + str(i) + '_prior'] = 'student_t'
        continuity['dsfr' + str(i) + '_prior_scale'] = (
            0.5  # Defaults to 0.3 (Leja19), we aim for a broader sample
        )
        continuity['dsfr' + str(i) + '_prior_df'] = (
            2  # Defaults to this value as in Leja19, but can be set
        )

    fit_params['continuity'] = continuity

    fit_params['dust'] = {
        'type': 'Cardelli',
        'Av': (0.0, 2.0),
        'eta': 2.0,
    }
    fit_params['nebular'] = {'logU': (-3.5, -2.0)}
    fit_params['t_bc'] = 0.02

    return fit_params

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.plot_conditions == ['detected']: args.plot_conditions = ['EW', 'mass', 'PA']
    args.plot_conditions = f'allpar_venn_{",".join(args.plot_conditions)}'
    if args.fontsize == 10: args.fontsize = 15
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)

    # --------------declaring all paths-----------------
    args.cosmos2020_filename = args.input_dir / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3.fits'
    photcat_filename = args.output_dir / f'{args.plot_conditions}_passage_cosmos_fluxes.csv'
    filter_dir = args.input_dir / 'COSMOS' / 'transmission_curves'
    pipes_dir = args.output_dir / 'pipes'
    pipes_dir.mkdir(exist_ok=True, parents=True)

    # -----------reading in flux and transmission files------------------
    df_int_filename = args.output_dir / f'{args.plot_conditions}_df.txt'
    df_int = pd.read_csv(df_int_filename)

    df_fluxes = get_flux_catalog(photcat_filename, df_int, args)
    df_fluxes = df_fluxes.sort_values(by='objid')
    fluxcols = [item for item in df_fluxes.columns if '_sci' in item]

    df_trans = read_filter_transmission(filter_dir, fluxcols, args)

    # ----------plotting (for tests)-----------------
    if args.plot_transmission: fig = plot_filter_transmission(df_trans, args, x_scale='log')
    if args.plot_SED: fig2 = plot_SED(df_fluxes, df_trans, args, x_scale='log')
    if args.plot_filter_cutouts: fig3 = plot_filter_cutouts(df_fluxes, args)
    if args.plot_niriss_direct: fig3 = plot_niriss_direct(df_fluxes, args)

    # ---------discarding some unusable flux columns-----------------------
    photcat_filename_sed = Path(str(photcat_filename).replace('.csv', '_for_bagpipe.csv'))
    if not os.path.isfile(photcat_filename_sed) or args.clobber:
        snr_thresh = 10
        flux_thresh = 0.05 # in uJy

        for fluxcol in fluxcols:
            snr = df_fluxes[fluxcol] / df_fluxes[fluxcol.replace('_sci', '_err')]
            flux = df_fluxes[fluxcol]
            snr = snr[np.isfinite(snr)]
            flux = flux[np.isfinite(flux)]
            if np.array(snr < snr_thresh).all() or np.array(flux < flux_thresh).all():
                print(f'Dropping {fluxcol}..')
                df_fluxes.drop(fluxcol, axis=1, inplace=True)
                df_fluxes.drop(fluxcol.replace('_sci', '_err'), axis=1, inplace=True)

        df_fluxes.to_csv(photcat_filename_sed, index=None)
        print(f'Written {photcat_filename_sed} with only the reliable flux columns.')
    else:
        print(f'Reading from existing {photcat_filename_sed}')

    df_sed = pd.read_csv(photcat_filename_sed)
    df_sed = df_sed.sort_values(by='objid')
    filter_list = [str(filter_dir) + '/' + item[:item.lower().find('_sci')] + '.txt' for item in df_sed.columns if '_sci' in item]
    print(f'Resultant photcat has {len(filter_list)} filters')

    # ----------SED fitting: the following part of the code is heavily borrowed from P J Watson-----------------
    if args.fit_sed:
        os.chdir(args.output_dir)
        load_fn = partial(load_photom_bagpipes, phot_cat=photcat_filename_sed, id_colname='objid', zeropoint=28.9)

        # --------create columns to store stellar masses---------------
        new_columns_dict = {'log_mass_bgp':('stellar_mass', False), 'z_bgp': ('redshift', False), 'log_sfr_bgp':('sfr', True)} # dict of new_label:(quantity_in_bagpipe, needs_to_be_log)
        new_columns = np.hstack([[item, item + '_u'] for item in list(new_columns_dict.keys())])
        for thiscol in new_columns: df_int[thiscol] = np.zeros(len(df_int))

        # ---------Loop over the objects-------------
        for index, obj in df_sed.iterrows():
            print(f'\nLooping over object {index + 1} of {len(df_sed)}..')
            fit_params = generate_fit_params(obj_z=obj['redshift'], z_range=0.01, num_age_bins=5, min_age_bin=30) # Generate the fit parameters

            galaxy = bagpipes.galaxy(ID=obj['objid'], load_data=load_fn, filt_list=filter_list, spectrum_exists=False) # Load the data for this object

            fit = bagpipes.fit(galaxy=galaxy, fit_instructions=fit_params, run=args.run) # Fit this galaxy
            fit.fit(verbose=True, sampler='nautilus', pool=4)

            # ---------Make some plots---------
            fig = fit.plot_spectrum_posterior(save=True, show=True)
            fig = fit.plot_sfh_posterior(save=True, show=True)
            fig = fit.plot_corner(save=True, show=True)

            # --------Save the stellar masses---------------
            index_in_df_int = df_int[df_int['objid'] == obj['objid']].index[0]
            for thisquant in list(new_columns_dict.keys()):
                peak_value = sci_mode(fit.posterior.samples[new_columns_dict[thisquant][0]]).mode
                low_value, up_value = np.percentile(fit.posterior.samples[new_columns_dict[thisquant][0]], (16, 84))
                err_value = (up_value - low_value) / 2
                quant = ufloat(peak_value, err_value)
                if new_columns_dict[thisquant][1]: quant = unp.log10(quant)
                df_int.loc[index_in_df_int, thisquant] = unp.nominal_values(quant)
                df_int.loc[index_in_df_int, thisquant + '_u'] = unp.std_devs(quant)

        os.chdir(args.code_dir)

        # ------writing modified df with stellar masses etc-------------------
        df_int_filename_sed = Path(str(df_int_filename).replace('.txt', f'_withSED_{args.run}.csv'))
        df_int.to_csv(df_int_filename_sed, index=None)
        print(f'Added SED results to df and saved in {df_int_filename_sed}.')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
