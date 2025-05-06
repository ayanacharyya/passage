'''
    Filename: run_glass_sed.py
    Notes: Computes stellar masses for a given list of GLASS galaxies that have fluxes from UNCOVER, and run SED fitting
    Author : Ayan
    Created: 05-05-25
    Example: run run_glass_sed.py --use_only_bands acs,wfc3,niriss,nircam,miri --run for_paper_only_st --clobber_sed_photcat
             run run_glass_sed.py --use_only_bands acs,wfc3,niriss,nircam,miri --run for_paper_only_st --plot_restframe --fit_sed --test_sed 1721
'''

from header import *
from util import *
from run_passage_sed import run_bagpipes   
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_uncover_photometry(df, photcat_dir, aperture=1, idcol='id'):
    '''
    Merges all available flux columns for a given apeture from the UNCOVER catalog to a given dataframe
    Returns merged dataframe
    '''
    print(f'Merging UNCOVER catalog for aperture {aperture} arcsec..')
    photcat_filename = photcat_dir / f'UNCOVER_DR3_LW_D{int(aperture * 100):3d}_catalog.fits'
    data = fits.open(photcat_filename)
    dfu = Table(data[1].data).to_pandas()

    ten_nJy_to_ujy_factor = 1e-2 # 10 nJy = 1e-8 Jy = 0.01 uJy
    filter_intrument_dict = {'f435w':'ACS', 'f606w':'ACS', 'f814w':'ACS', 'f070w':'NIRCAM', 'f090w':'NIRCAM', 'f105w':'WFC3', 'f115w':'NIRCAM', \
                             'f125w':'WFC3', 'f140w':'WFC3', 'f140m':'NIRCAM', 'f150w':'NIRCAM', 'f160w':'WFC3', 'f162m':'XXX', 'f182m':'NIRCAM', \
                             'f200w':'NIRCAM', 'f210m':'NIRCAM', 'f250m':'NIRCAM', 'f277w':'NIRCAM', 'f300m':'NIRCAM', 'f335m':'NIRCAM', 'f356w':'NIRCAM', \
                             'f360m':'NIRCAM', 'f410m':'NIRCAM', 'f430m':'NIRCAM', 'f444w':'NIRCAM', 'f460m':'NIRCAM', 'f480m':'NIRCAM'}

    filter_list = [item for item in dfu.columns if 'f_f' in item]
    cols_to_extract = np.hstack([['id', 'ra', 'dec'], np.hstack([[item, item.replace('f_', 'e_')] for item in filter_list])])
    dfu = dfu[cols_to_extract]

    print(f'\nDoing cross-matching between GLASS-NIRISS and UNCOVER catalogs..')
    df_crossmatch = get_crossmatch(df, dfu, sep_threshold=1.0, df1_idcol=idcol, df2_idcol='id')
    df_crossmatch = df_crossmatch.rename(columns={'df1_id':idcol, 'df2_id':'ID_UNCOVER'})
    df_crossmatch = pd.merge(df_crossmatch[[idcol,  'ID_UNCOVER']], dfu, left_on='ID_UNCOVER', right_on='id', how='inner').drop(['id', 'ra', 'dec'], axis=1)
    df = pd.merge(df, df_crossmatch)
    
    for thiscol in filter_list:
        df[thiscol] = df[thiscol] * ten_nJy_to_ujy_factor
        df[thiscol.replace('f_', 'e_')] = df[thiscol.replace('f_', 'e_')] * ten_nJy_to_ujy_factor

    df=df.dropna(axis=1, how='all') # dropping columns that have NaN values for all objects

    df = df.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:].upper() + '_err') for item in df.columns if'e_' in item]))
    df = df.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:].upper() + '_sci') for item in df.columns if'f_' in item]))

    return df

# --------------------------------------------------------------------------------------------------------------------
def merge_glass_speccat(speccat_filename, idcol='id', objids=None):
    '''
    Reads in the GLASS-NIRISS speccat, and merges redshift with given dataframe
    Returns merged dataframe
    '''
    print(f'Reading in speccat from {speccat_filename}..')
    tab = Table(fits.open(speccat_filename)[1].data)
    if 'cdf_z' in tab.columns: tab.remove_column('cdf_z') # because multi dimension column does not agree well with pandas
    df = tab.to_pandas()

    if objids is not None: df = df[df[idcol].isin(objids)]

    cols_to_extract = [idcol, 'Z_NIRISS', 'Z_FLAG']
    df = df[cols_to_extract]
    df = df.rename(columns={'Z_NIRISS':'redshift'})

    return df

# --------------------------------------------------------------------------------------------------------------------
def merge_glass_photcat(df, photcat_filename, aperture=1.0, idcol='id'):
    '''
    Reads in the GLASS photcat, with flux columns corresponding to a given aperture, in to a dataframe
    Returns merged dataframe
    '''
    print(f'Reading in photcat from {photcat_filename}..')
    hdu = fits.open(photcat_filename)[1]
    df_phot = Table(hdu.data).to_pandas().rename(columns={'id': idcol})

    # ---------determining the aperture id------------------------
    index = 0
    while True:
        try:
            this_ap = np.round(hdu.header[f'ASEC_{index}'], 1)
        except KeyError:
            print(f'{aperture} size not found in any of ASEC_0 to ASEC_{index - 1}')
            break
        if this_ap == aperture:
            aper_ind = index
            print(f'Found APER_{aper_ind} to match {aperture}')
            break
        index += 1

    # -------determining flux and fluxerr columns from passage-------
    fluxcols = [item for item in df_phot.columns if '_flux' in item and '_fluxerr' not in item]
    filters = np.unique([item[: item.find('_flux')] for item in fluxcols])
    cols_to_extract = np.hstack([[f'{item}_flux_aper_{aper_ind:0d}', f'{item}_fluxerr_aper_{aper_ind:0d}'] for item in filters])
    
    df_phot = df_phot[np.hstack([[idcol, 'ra', 'dec'], cols_to_extract])]
    df = pd.merge(df, df_phot, on=idcol, how='left')
    
    for thisfilter in filters:
        df= df.rename(columns={f'{thisfilter}_flux_aper_{aper_ind:0d}': f'NIRISS_{thisfilter.upper()[:-1]}_sci', f'{thisfilter}_fluxerr_aper_{aper_ind:0d}': f'NIRISS_{thisfilter.upper()[:-1]}_err'})

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    aperture = 1.
    idcol = 'ID_NIRISS'
    glass_object_ids = [1333, 1721, 1983, 2128]

    filter_dir = args.input_dir / 'COSMOS' / 'transmission_curves'
    uncover_catalog_dir = args.root_dir / 'glass_data'
    glass_photcat_filename = args.root_dir / 'glass_data' / 'glass-a2744_phot.fits'
    glass_speccat_filename = args.root_dir / 'glass_data' / 'a2744_spec_cat_niriss_20250401.fits'
    photcat_filename_sed = uncover_catalog_dir / f'GLASS_UNCOVER_photometry_{args.run}.csv'

    # ----------making the photometric catalog for SED fitting-----------------
    if not os.path.exists(photcat_filename_sed) or args.clobber_sed_photcat:
        print(f'Not found {photcat_filename_sed}, so creating new..')

        # -----------gathering the catalog for SED fitting------------------
        df = merge_glass_speccat(glass_speccat_filename, idcol=idcol, objids=glass_object_ids)
        df = merge_glass_photcat(df, glass_photcat_filename, aperture=aperture, idcol=idcol)
        df = get_uncover_photometry(df, uncover_catalog_dir, aperture=aperture, idcol=idcol)

        # ---using only specific bands-------------
        if args.use_only_bands is not None:
            print(f'\nSelecting only bands with {args.use_only_bands} in the filter name..')
            use_bands = args.use_only_bands.split(',')
            fluxcols = [item for item in df.columns if '_sci' in item]

            useless_fluxcols = [col for col in fluxcols if np.array([band.lower() not in col.lower() for band in use_bands]).all()]
            print(f'Therefore keeping only {list(set(fluxcols) - set(useless_fluxcols))} bands\n')
            useless_errcols = [item.replace('_sci', '_err') for item in useless_fluxcols]
            df.drop(useless_fluxcols, axis=1, inplace=True)
            df.drop(useless_errcols, axis=1, inplace=True)

        # ------writing out the catalog for SED fitting----------------
        df[idcol] = df[idcol].astype(int)
        df.to_csv(photcat_filename_sed, index=None)
        print(f'Saved SED photcat in {photcat_filename_sed}')
    else:
        print(f'Reading from existing {photcat_filename_sed}..')

    df = pd.read_csv(photcat_filename_sed)
    df = df[df[idcol].isin(glass_object_ids)]

    # --------running SED fitting---------------------
    if args.fit_sed:
        dfm = run_bagpipes(photcat_filename_sed, filter_dir, args, idcol=idcol)
    else:
        print(f'Photcat for SED produced; use --fit_sed to actually run SED fitting')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

