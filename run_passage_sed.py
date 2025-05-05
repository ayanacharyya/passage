'''
    Filename: run_passage_sed.py
    Notes: Computes stellar masses for a given list of PASSAGE galaxies that also have fluxes from COSMOS2020, and run SED fitting
    Author : Ayan
    Created: 19-08-24
    Example: run run_passage_sed.py --field Par028 --use_only_bands acs,niriss,nircam,miri --run only_st
             run run_passage_sed.py --field Par028 --use_only_bands acs,niriss,nircam,miri --run only_st --plot_restframe --fit_sed --test_sed 1303
'''

from header import *
from util import *
from compute_stellar_masses import load_photom_bagpipes, generate_fit_params    
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_number(item):
    '''
    Extracts numbers (floats or ints) from strings
    Returns list of numbers
    '''
    if item == 0 or item == 1 or item == -1: return int(item)
    else: return float(item)

# --------------------------------------------------------------------------------------------------------------------
def get_line_finding_zcat(df, zcat_filename, idcol='objid'):
    '''
    Read in redshift catalog produced by Kalina's line finding code (for a given user, before reconciliation) and
    Merges all available objects in the given redshift catalog with the given phot catalog
    Returns merged dataframe
    '''
    outfilename = zcat_filename.parent / Path(str(zcat_filename.stem) + '.csv')

    if not os.path.exists(outfilename) or args.clobber:
        print(f'Not found {outfilename}, so making a new one..')
        f = open(zcat_filename, 'r')
        lines = f.readlines()

        # ----------to get all column labels-----------
        column_labels, counter = [], 0
        for line in lines:
            if line[0] == '#':
                start_index = line.find(next(filter(str.isalpha, line)))
                end_index = line.find('\n') - 1
                label = line[start_index: end_index]
                column_labels.append(label)
                counter += 1
            else:
                break
        
        # ----------to insert a new column label at a specific location-----------
        insert_idx = np.where(np.array(column_labels) == 'fwhm_g141_error')[0][0]
        column_labels.insert(insert_idx, 'comp_fit_flag')

        # ----------making new dataframe-----------
        dfz = pd.DataFrame()
        for line in lines[counter:]: # looping over each line in file
            if line[0] == '#': continue
            this_row = []
            for item in line.split(): # looping over each item in each line
                if 'e+' in item or 'e-' in item: # exponential notations require special treatment...
                    if len(item.split('e')[1]) > 3: # ...if they have an additional digit at the end that was supposed to be a different column
                        this_row.append(get_number(item[:-1]))
                        this_row.append(get_number(item[-1]))
                    else:
                        this_row.append(get_number(item))
                elif 'False' in item: # items with False or True combined with digits in the same string require special handling to split them
                    items = item.split('False')
                    this_row.append(get_number(items[0]))
                    this_row.append(False)
                elif 'True' in item:
                    items = item.split('True')
                    this_row.append(get_number(items[0]))
                    this_row.append(True)
                else:
                    this_row.append(get_number(item))
            this_df = pd.DataFrame(dict(zip(column_labels, this_row)), index=[0])
            dfz = pd.concat([dfz, this_df])

        # ------writing out the dataframe---------
        dfz[idcol] = dfz[idcol].astype(int)
        dfz.to_csv(outfilename, index=None)
        print(f'Saved as new table in {outfilename}')
    else:
        print(f'Reading existing csv linelist from {outfilename}')

    # ------reading z catalog and merging---------------   
    dfz = pd.read_csv(outfilename).rename(columns={idcol: 'PASSAGE_ID'})
    dfz['PASSAGE_ID'] = dfz['PASSAGE_ID'].astype(int)
    df = df.merge(dfz[['PASSAGE_ID', 'redshift', 'redshift_error']], on='PASSAGE_ID', how='inner')
    df = df.rename(columns = {'redshift_error': 'redshift_u'})

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_cosmos2020_fluxes(df, cosmos_photcat_dir, idcol='id'):
    '''
    Merges all available flux columns from the COSMOS2020 catalog to a given dataframe
    Returns merged dataframe
    '''

    for thisfield in np.unique(df['field']):
        print(f'Merging COSMOS 2020 catalog for field {thisfield}..')
        cosmos2020_photcat_for_thisfield_filename = cosmos_photcat_dir / f'cosmos2020_objects_in_{thisfield}.fits'
        cosmos2020_photcat_for_thisfield = read_COSMOS2020_catalog(args=None, filename=cosmos2020_photcat_for_thisfield_filename)
        cosmos2020_photcat_for_thisfield = cosmos2020_photcat_for_thisfield.rename(columns={idcol: 'COSMOS2020_ID'})
        cosmos2020_photcat_for_thisfield = cosmos2020_photcat_for_thisfield.drop(['ra', 'dec'], axis=1)
        df = pd.merge(df, cosmos2020_photcat_for_thisfield, on='passage_id', how='left')

    df = remove_excess_columns(df)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_cosmosweb_fluxes(df, cosmos_photcat_dir, aperture=1, idcol='id'):
    '''
    Merges all available flux columns for a given apeture from the COSMOSWeb catalog to a given dataframe
    Returns merged dataframe
    '''

    filter_intrument_dict = {'F115W':'NIRCAM', 'F150W':'NIRCAM', 'F277W':'NIRCAM', 'F444W':'NIRCAM', 'F770W':'MIRI'}
    ergs_s_cm2_hz_to_ujy_factor = 1e29 # 1 ergs/s/cm^2/Hz = 10^29 uJy

    for thisfield in np.unique(df['field']):
        print(f'Merging COSMOS Webb catalog for field {thisfield}..') 
        cosmoswebb_photcat_for_thisfield_filename = cosmos_photcat_dir / f'cosmoswebb_objects_in_{thisfield}.fits'
        cosmoswebb_photcat_for_thisfield = read_COSMOSWebb_catalog(args=None, filename=cosmoswebb_photcat_for_thisfield_filename, aperture=aperture)
        nircam_fluxcols = [item for item in cosmoswebb_photcat_for_thisfield.columns if 'FLUX_APER' in item]
        nircam_errcols = [item for item in cosmoswebb_photcat_for_thisfield.columns if 'FLUX_ERR_APER' in item]
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield[np.hstack([nircam_fluxcols, nircam_errcols, ['passage_id', idcol]])]
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.replace({-998: np.nan})
        for thiscol in nircam_fluxcols:
            cosmoswebb_photcat_for_thisfield[thiscol] = cosmoswebb_photcat_for_thisfield[thiscol] * ergs_s_cm2_hz_to_ujy_factor # converting from ergs/s/cm^2/Hz to micro Jansky
            cosmoswebb_photcat_for_thisfield[thiscol.replace('FLUX', 'FLUX_ERR')] = cosmoswebb_photcat_for_thisfield[thiscol.replace('FLUX', 'FLUX_ERR')] * ergs_s_cm2_hz_to_ujy_factor  # converting from ergs/s/cm^2/Hz to micro Jansky

        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns={idcol: 'COSMOSWeb_ID'})
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:] + '_FLUXERR') for item in cosmoswebb_photcat_for_thisfield.columns if'FLUX_ERR_APER' in item]))
        cosmoswebb_photcat_for_thisfield = cosmoswebb_photcat_for_thisfield.rename(columns=dict([(item, filter_intrument_dict[item.split('_')[-1]] + '_' + item[-5:] + '_FLUX') for item in cosmoswebb_photcat_for_thisfield.columns if'FLUX_APER' in item]))
        df = pd.merge(df, cosmoswebb_photcat_for_thisfield, on='passage_id', how='inner')
    
    df = remove_excess_columns(df)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_pjw_zcat(df, zcat_filename, idcol='objid'):
    '''
    Merges all available objects in the given redshift catalog (made from KVN's line finding catalog by PJW) with the given phot catalog
    Returns merged dataframe
    '''

    i = str(zcat_filename).find('Par')
    field = str(zcat_filename)[i:i+6]


    print(f'Reading in zcat from {zcat_filename} for {field}..') 
    dfz = Table(fits.open(zcat_filename)[1].data).to_pandas().rename(columns={idcol: 'PASSAGE_ID'})
    dfz['PASSAGE_ID'] = dfz['PASSAGE_ID'].astype(int)

    df = df.merge(dfz[['PASSAGE_ID', 'redshift', 'redshift_error']], on='PASSAGE_ID', how='inner')

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_passage_photcat(photcat_filename, aperture=1.0, idcol='id'):
    '''
    Reads in the PASSAGE photcat, with flux columns corresponding to a given aperture, in to a dataframe
    Returns merged dataframe
    '''


    i = str(photcat_filename).find('Par')
    field = str(photcat_filename)[i:i+6]

    print(f'Reading in photcat from {photcat_filename} for {field}..')
    hdu = fits.open(photcat_filename)[1]
    df = Table(hdu.data).to_pandas().rename(columns={idcol: 'PASSAGE_ID'})
    
    # ---------determining the aperture id------------------------
    for index in range (10):
        if hdu.header[f'ASEC_{index}'] == aperture:
            aper_ind = index
            print(f'Found APER_{aper_ind} to match {aperture}')
            break

    # -------determining flux and fluxerr columns from passage-------
    fluxcols = [item for item in df.columns if '_flux' in item and '_fluxerr' not in item]
    passage_filters = np.unique([item[: item.find('_flux')] for item in fluxcols])

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
def run_bagpipes(photcat_filename_sed, filter_dir, args, idcol='PASSAGE_ID'):
    '''
    Run bagpipes on each object in given phot catalog photcat_filename, save stellar mass in the catalog and write out the resulting catalog
    Returns resulting catalog, containing stellar mass
    '''
    # -------reading and polishing phot catalog------------------
    df = pd.read_csv(photcat_filename_sed)
    filter_list = [str(filter_dir) + '/' + item[:item.lower().find('_sci')] + '.txt' for item in df.columns if '_sci' in item]
    print(f'Resultant photcat has {len(filter_list)} filters')

    # ---------loading bagpipes-related stuff------------------------
    os.chdir(args.output_dir)
    load_fn = partial(load_photom_bagpipes, phot_cat=photcat_filename_sed, id_colname=idcol, zeropoint=28.9)

    # --------create columns to store stellar masses---------------
    new_columns_dict = {'log_mass_bgp':('stellar_mass', False), 'z_bgp': ('redshift', False), 'log_sfr_bgp':('sfr', True)} # dict of new_label:(quantity_in_bagpipe, needs_to_be_log)
    new_columns = np.hstack([[item, item + '_u'] for item in list(new_columns_dict.keys())])
    for thiscol in new_columns: df[thiscol] = np.zeros(len(df))

    # ---------Loop over the objects-------------
    if args.test_sed is not None:
        df = df[df['PASSAGE_ID'] == args.test_sed].reset_index(drop=True)
        print(f'Only runing on object {args.test_sed} as a test; for doing SED all objects, remove --test_sed and re-run')
    
    for index, obj in df.iterrows():
        print(f'\nLooping over object {index + 1} of {len(df)}..')
        fit_params = generate_fit_params(obj_z=obj['redshift'], z_range=0.01, num_age_bins=5, min_age_bin=30) # Generate the fit parameters

        galaxy = bagpipes.galaxy(ID=int(obj[idcol]), load_data=load_fn, filt_list=filter_list, spectrum_exists=False) # Load the data for this object
        fit = bagpipes.fit(galaxy=galaxy, fit_instructions=fit_params, run=args.run) # Fit this galaxy
        fit.fit(verbose=True, sampler='nautilus', pool=args.ncpus)

        # --------converting everything to restframe----------------
        if args.plot_restframe:
            fit.posterior.get_advanced_quantities()
            redshift = np.median(fit.posterior.samples['redshift'])
            fit.galaxy.photometry[:, 0] = fit.galaxy.photometry[:, 0] / (1 + redshift)
            if fit.galaxy.spectrum_exists: fit.galaxy.spectrum[:, 0] = fit.galaxy.spectrum[:, 0] / (1 + redshift)
            fit.galaxy.filter_set.eff_wavs = fit.galaxy.filter_set.eff_wavs / (1 + redshift)
            fit.posterior.model_galaxy.wavelengths = fit.posterior.model_galaxy.wavelengths / (1 + redshift)

        # ---------Make some plots---------
        fig, ax = fit.plot_spectrum_posterior(save=True, show=True, log_x=True, xlim=[2.7, 4.5], ylim=[0, 6])
        fig, ax = fit.plot_spectrum_posterior(save=True, show=True, log_x=False, xlim=[500, 30000], ylim=[0, 6])
        fig = fit.plot_sfh_posterior(save=True, show=True, xlim=None, ylim=[0, 10])
        fig = fit.plot_corner(save=True, show=True)

        # --------Save the stellar masses---------------
        for thisquant in list(new_columns_dict.keys()):
            peak_value = sci_mode(fit.posterior.samples[new_columns_dict[thisquant][0]]).mode
            low_value, up_value = np.percentile(fit.posterior.samples[new_columns_dict[thisquant][0]], (16, 84))
            err_value = (up_value - low_value) / 2
            quant = ufloat(peak_value, err_value)
            if new_columns_dict[thisquant][1]: quant = unp.log10(quant)
            df.loc[index, thisquant] = unp.nominal_values(quant)
            df.loc[index, thisquant + '_u'] = unp.std_devs(quant)

    os.chdir(args.code_dir)

    # ------writing modified df with stellar masses etc-------------------
    df_filename_sed = Path(str(photcat_filename_sed).replace('.csv', f'_withSED_{args.run}.csv'))
    df.to_csv(df_filename_sed, index=None)
    print(f'Added SED results to df and saved in {df_filename_sed}.')

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    filter_dir = args.input_dir / 'COSMOS' / 'transmission_curves'
    passage_photcat_dir = args.input_dir / 'data' / args.field / 'DATA' / 'DIRECT_GRISM' 
    passage_photcat_filename = passage_photcat_dir / f'{args.field}_photcat.fits' # from Vihang's upload on Box
    passage_zcat_filename = args.input_dir / 'data' / args.field / 'line_finding' / f'{args.field}lines_catalog_scarlata.csv' # improved version of what was from KVN
    aperture = 1.

    # -----------gathering the catalog for SED fitting------------------
    df = get_passage_photcat(passage_photcat_filename, aperture=aperture, idcol='id')
    df = get_line_finding_zcat(df, passage_zcat_filename, idcol='objid')
    df = get_cosmosweb_fluxes(df, passage_photcat_dir, aperture=aperture, idcol='id')
    df = get_cosmos2020_fluxes(df, passage_photcat_dir, idcol='id')

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
    df.columns = df.columns.str.replace('_FLUXERR', '_err', regex=True)
    df.columns = df.columns.str.replace('_FLUX', '_sci', regex=True)
    photcat_filename_sed = passage_photcat_dir / f'{args.field}_PASSAGE_COSMOS2020_COSMOSWeb_fluxes_{args.run}.csv'
    df.to_csv(photcat_filename_sed, index=None)
    print(f'Saved SED photcat in {photcat_filename_sed}')

    # --------running SED fitting---------------------
    if args.fit_sed:
        dfm = run_bagpipes(photcat_filename_sed, filter_dir, args, idcol='PASSAGE_ID')
    else:
        print(f'Photcat for SED produced; use --fit_sed to actually run SED fitting')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

