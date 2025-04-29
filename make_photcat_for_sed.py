'''
    Filename: compute_stellar_masses.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020
    Author : Ayan
    Created: 19-08-24
    Example: run make_photcat_for_sed.py --field Par028 --use_only_bands acs,niriss,nircam,miri
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_cosmos2020_fluxes(df):
    '''
    Merges all available flux columns from the COSMOS2020 catalog to a given dataframe
    Returns merged dataframe
    '''

    for thisfield in np.unique(df['field']):
        print(f'Merging COSMOS 2020 catalog for field {thisfield}..')
        cosmos2020_photcat_for_thisfield_filename = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5') / 'COSMOS' / f'cosmos2020_objects_in_{thisfield}.fits'
        cosmos2020_photcat_for_thisfield = read_COSMOS2020_catalog(args=None, filename=cosmos2020_photcat_for_thisfield_filename)
        cosmos2020_photcat_for_thisfield = cosmos2020_photcat_for_thisfield.rename(columns={'id': 'COSMOS2020_ID'})
        cosmos2020_photcat_for_thisfield = cosmos2020_photcat_for_thisfield.drop(['ra', 'dec'], axis=1)
        df = pd.merge(df, cosmos2020_photcat_for_thisfield, on='passage_id', how='left')

    df = remove_excess_columns(df)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_cosmosweb_fluxes(df, aperture=1):
    '''
    Merges all available flux columns for a given apeture from the COSMOSWeb catalog to a given dataframe
    Returns merged dataframe
    '''

    filter_intrument_dict = {'F115W':'NIRCAM', 'F150W':'NIRCAM', 'F277W':'NIRCAM', 'F444W':'NIRCAM', 'F770W':'MIRI'}
    ergs_s_cm2_hz_to_ujy_factor = 1e29 # 1 ergs/s/cm^2/Hz = 10^29 uJy

    for thisfield in np.unique(df['field']):
        print(f'Merging COSMOS Webb catalog for field {thisfield}..') 
        cosmoswebb_photcat_for_thisfield_filename = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5') / 'COSMOS' / f'cosmoswebb_objects_in_{thisfield}.fits'
        cosmoswebb_photcat_for_thisfield = read_COSMOSWebb_catalog(args=None, filename=cosmoswebb_photcat_for_thisfield_filename, aperture=aperture)
        nircam_fluxcols = [item for item in cosmoswebb_photcat_for_thisfield.columns if 'FLUX_APER' in item]
        nircam_errcols = [item for item in cosmoswebb_photcat_for_thisfield.columns if 'FLUX_ERR_APER' in item]
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield[np.hstack([nircam_fluxcols, nircam_errcols, ['passage_id', 'id']])]
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.replace({-998: np.nan})
        for thiscol in nircam_fluxcols:
            cosmoswebb_photcat_for_thisfield[thiscol] = cosmoswebb_photcat_for_thisfield[thiscol] * ergs_s_cm2_hz_to_ujy_factor # converting from ergs/s/cm^2/Hz to micro Jansky
            cosmoswebb_photcat_for_thisfield[thiscol.replace('FLUX', 'FLUX_ERR')] = cosmoswebb_photcat_for_thisfield[thiscol.replace('FLUX', 'FLUX_ERR')] * ergs_s_cm2_hz_to_ujy_factor  # converting from ergs/s/cm^2/Hz to micro Jansky

        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns={'id': 'COSMOSWeb_ID'})
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:] + '_FLUXERR') for item in cosmoswebb_photcat_for_thisfield.columns if'FLUX_ERR_APER' in item]))
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:] + '_FLUX') for item in cosmoswebb_photcat_for_thisfield.columns if'FLUX_APER' in item]))
        df = pd.merge(df, cosmoswebb_photcat_for_thisfield, on='passage_id', how='inner')
    
    df = remove_excess_columns(df)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_passage_photcat(photcat_filename, aperture=1.0):
    '''
    Reads in the PASSAGE photcat, with flux columns corresponding to a given aperture, in to a dataframe
    Returns merged dataframe
    '''
    i = str(photcat_filename).find('Par')
    field = str(photcat_filename)[i:i+6]

    print(f'Reading in photcat from {photcat_filename} for {field}..') 
    df = Table(fits.open(photcat_filename)[1].data).to_pandas().rename(columns={'id': 'PASSAGE_ID'})
    
    # -------determining flux and fluxerr columns from passage-------
    fluxcols = [item for item in df.columns if '_flux' in item and '_fluxerr' not in item]
    passage_filters = np.unique([item[: item.find('_flux')] for item in fluxcols])
    aper_ind = np.where(np.array([0.36, 0.5, 0.7, 1, 1.2, 1.5, 3.0]) == aperture)[0][0]

    cols_to_extract = np.hstack([[f'{item}_flux_aper_{aper_ind:0d}', f'{item}_fluxerr_aper_{aper_ind:0d}'] for item in passage_filters])
    df = df[np.hstack([['PASSAGE_ID', 'ra', 'dec'], cols_to_extract])]
    df['field'] = field
    df['passage_id'] = df['field'].astype(str) + '-' + df['PASSAGE_ID'].astype(str)
    for thisfilter in passage_filters:
        df= df.rename(columns={f'{thisfilter}_flux_aper_{aper_ind:0d}': f'NIRISS_{thisfilter.upper()}_FLUX', f'{thisfilter}_fluxerr_aper_{aper_ind:0d}': f'NIRISS_{thisfilter.upper()}_FLUXERR'})

    df = remove_excess_columns(df)
    return df

# --------------------------------------------------------------------------------------------------------------------
def remove_excess_columns(df):
    '''
    Removes columns other than those to retainfrom a dataframe
    Returns merged dataframe
    '''
    cols_to_retain = ['_FLUX', '_FLUXERR', 'ID', 'id', 'field', 'redshift', 'ra', 'dec']
    cols_to_drop = [item for item in df.columns if not np.array([item2 in item for item2 in cols_to_retain]).any()]
    df = df.drop(cols_to_drop, axis=1)

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    passage_photcat_filename = args.input_dir / 'data/Par028/DATA/DIRECT_GRISM/Par028_photcat.fits'
    aperture = 1.

    df = get_passage_photcat(passage_photcat_filename, aperture=aperture)
    df = get_cosmosweb_fluxes(df, aperture=aperture)
    df = get_cosmos2020_fluxes(df)

    # ---using only specific bands-------------
    if args.use_only_bands is not None:
        print(f'\nSelecting only bands with {args.use_only_bands} in the filter name..')
        use_bands = args.use_only_bands.split(',')
        fluxcols = [item for item in df.columns if '_FLUX' in item and '_FLUXERR' not in item]

        useless_fluxcols = [col for col in fluxcols if np.array([band.lower() not in col.lower() for band in use_bands]).all()]
        print(f'Therefore keeping only {list(set(fluxcols) - set(useless_fluxcols))} bands\n')
        useless_errcols = [item.replace('_FLUX', '_FLUXERR') for item in useless_fluxcols]
        df.drop(useless_fluxcols, axis=1, inplace=True)
        df.drop(useless_errcols, axis=1, inplace=True)

    # ------writing out the catalog for SED fitting----------------
    df = df.drop('passage_id', axis=1)
    photcat_filename_sed = passage_photcat_filename.parent / f'{args.field}_PASSAGE_COSMOS2020_COSMOSWeb_fluxes.csv'
    df.to_csv(photcat_filename_sed, index=None)
    print(f'Saved SED photcat in {photcat_filename_sed}')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

