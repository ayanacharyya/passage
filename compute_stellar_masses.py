'''
    Filename: compute_stellar_masses.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020
    Author : Ayan
    Created: 19-08-24
    Example: run compute_stellar_masses.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
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
def plot_filter_transmission(df_master, args, x_scale='linear', color_by_wave=True):
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

    plt.legend(fontsize = 5)
    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel('Transmission')
    ax.set_xscale(x_scale)

    figname = args.input_dir / 'COSMOS' / 'transmission_curves' / 'filter_responses.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

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

    args.cosmos2020_filename = args.input_dir / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3.fits'
    df_int_filename = args.output_dir / f'allpar_venn_EW,mass,PA_df.txt'
    df_int = pd.read_csv(df_int_filename)
    df_fluxes = get_flux_catalog(df_int, args)

    fluxcols = [item for item in df_fluxes.columns if 'flux' in item.lower() and 'fluxerr' not in item.lower()]
    df_trans = read_filter_transmission(fluxcols, args)
    filters = pd.unique(df_trans['filter'])
    fig = plot_filter_transmission(df_trans, args, x_scale='log')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
