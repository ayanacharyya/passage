'''
    Filename: compute_stellar_masses.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020
    Author : Ayan
    Created: 19-08-24
    Example: run compute_stellar_masses.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run compute_stellar_masses.py --only_seg --arcsec_limit 1
             run compute_stellar_masses.py
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
def plot_cutouts(df_fluxes, args):
    '''
    Function to plot 2D image cutouts based on input dataframe of ra, dec and existing mosaic images
    Saves plot as png figure and
    Returns figure handle
    '''
    max_filters_per_page = 6
    cmap = 'viridis'
    image_dir = args.input_dir / 'COSMOS' / 'imaging'

    if args.fontsize == 10: args.fontsize = 15
    cutout_size = 2 * args.arcsec_limit # in arcsec
    fits_images = glob.glob(str(image_dir / '*.fits'))
    n_filters = len(fits_images)
    n_obj = len(df_fluxes)
    n_figs = int(np.ceil(n_filters / max_filters_per_page))

    figname = args.output_dir / f'{args.intersection_conditions}_all_cutouts.pdf'
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
    print(f'\nTotal {n_obj} x {n_filters} = {n_obj * n_filters} cutouts to plot, of size {cutout_size}" each; hence splitting in to {n_figs} figures..')
    fig_arr = []
    for fig_index in range(n_figs):
        print(f'\nMaking figure {fig_index + 1} of {n_figs}..')
        these_fits_images = fits_images[fig_index * max_filters_per_page : min((fig_index + 1) * max_filters_per_page, n_filters)] # slicing the fits image array

        fig, axes = plt.subplots(nrows=n_obj, ncols=min(n_filters, max_filters_per_page), figsize=(10, 8))
        fig.subplots_adjust(left=0.05, right=0.98, bottom=0.07, top=0.97, hspace=0.05, wspace=0.05)

        # ------looping over filters-------------
        for col_index, thisfile in enumerate(these_fits_images):
            print(f'Reading in file {os.path.split(thisfile)[-1]} which is {fig_index * max_filters_per_page + col_index + 1} of {n_filters}..')
            try:
                data = fits.open(thisfile)
                image = data[0].data
                header = data[0].header
            except Exception as e:
                print(f'Skipping {os.path.split(thisfile)[-1]} due to {e}')

                # --------making sure the tick labels still get marked properly--------
                for row_index, row in df_fluxes.iterrows():
                    ax = axes[row_index][col_index]
                    if row_index == 0: ax.text(0.98, 0.98, filter, c='k', ha='right', va='top', fontsize=args.fontsize / 2,transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                    if col_index == 0:
                        ax.set_ylabel('Dec (")', fontsize=args.fontsize / 2)
                        ax.tick_params(axis='y', which='major', labelsize=args.fontsize / 2)
                        ax.text(0.05, 0.05,  f'{row["redshift"]:.2f}', c='k', ha='left', va='bottom', fontsize=args.fontsize / 2, transform = ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                    else:
                        ax.set_yticklabels([])

                    if row_index == n_obj - 1:
                        ax.set_xlabel('RA (")', fontsize=args.fontsize / 2)
                        ax.tick_params(axis='x', which='major', labelsize=args.fontsize / 2)
                    else:
                        ax.set_xticklabels([])

                continue

            wcs_header = pywcs.WCS(header)
            filter = header['FILTER'] if 'FILTER' in header else ''

            # ------looping over objects-------------
            for row_index, row in df_fluxes.iterrows():
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

                # ------now plotting the cutout-------------
                ax = axes[row_index][col_index]
                ax.imshow(np.log10(cutout.data), origin='lower', extent=args.extent, cmap=cmap)
                if args.only_seg: ax.contour(seg_cutout_data_rebinned != row['objid'], levels=0, colors='k', extent=args.extent, linewidths=0.5)

                # ------annotations and tick labels-------------
                if row_index == 0: ax.text(0.98, 0.98, filter, c='k', ha='right', va='top', fontsize=args.fontsize/2, transform = ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                if col_index == 0:
                    ax.set_ylabel('Dec (")', fontsize=args.fontsize/2)
                    ax.tick_params(axis='y', which='major', labelsize=args.fontsize/2)
                    ax.text(0.05, 0.05, f'{row["redshift"]:.2f}', c='k', ha='left', va='bottom', fontsize=args.fontsize/2, transform = ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
                else:
                    ax.set_yticklabels([])

                if row_index == n_obj - 1:
                    ax.set_xlabel('RA (")', fontsize=args.fontsize/2)
                    ax.tick_params(axis='x', which='major', labelsize=args.fontsize/2)
                else:
                    ax.set_xticklabels([])

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
    #fig = plot_filter_transmission(df_trans, args, x_scale='log')
    #fig2 = plot_SED(df_fluxes, df_trans, args, x_scale='log')
    fig3 = plot_cutouts(df_fluxes, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
