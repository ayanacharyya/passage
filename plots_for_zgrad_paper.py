'''
    Filename: plots_for_zgrad_paper.py
    Notes: Makes ALL plots related to PASSAGE metallicity gradient paper and saves them separately; some intentional hard-coding and reduction of flexibility in this script, for posterity
    Author : Ayan
    Created: 09-04-25
    Example: run plots_for_zgrad_paper.py --phot_models nb --debug_Zsfr
             run plots_for_zgrad_paper.py --histbycol SNR
'''
from header import *
from util import *

from get_field_stats import get_crossmatch_with_cosmos, plot_venn, read_stats_df, make_set
from plot_mappings_grid import plot_ratio_grid, plot_ratio_model
from make_passage_plots import break_column_into_uncertainty, plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14, plot_MZR_literature
from make_diagnostic_maps import get_re, get_cutout, get_emission_line_map, get_offsets_from_center, annotate_PAs, get_linelist, get_EB_V, get_voronoi_bin_IDs, get_radial_bin_IDs, get_voronoi_bin_distances, get_AGN_func_methods, AGN_func, take_safe_log_ratio, overplot_AGN_line_on_BPT, get_distance_map, compute_SFR

plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.width'] = 1

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, figname, args):
    '''
    Saves a given figure handle as a given output filename
    '''

    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass
        plot_output_dir = args.root_dir / 'zgrad_paper_plots_fortalk'
    else:
        plot_output_dir = args.root_dir / 'zgrad_paper_plots'

    plot_output_dir.mkdir(exist_ok=True, parents=True)
    figname = plot_output_dir / figname
    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def load_full_df(fields, args, cosmos_name='web'):
    '''
    Loads and returns a dataframe with ALL objects in the fields listed within objlist
    Returns dataframe
    '''
    print(f'\nLoading full data frame for fields" {fields}..')

    df = pd.DataFrame()
    # ---------looping over fields-----------
    for index, args.field in enumerate(fields):
        args.filters = args.available_filters_for_field_dict[args.field]
        output_dir = get_output_path(args.field, args)
        df_filename = output_dir / f'{args.field}_all_diag_results.csv'
        thisdf = read_stats_df(df_filename, args)  # read in the stats dataframe
        thisdf['filters'] = ', '.join(args.filters)
        df = pd.concat([df, thisdf], ignore_index=True)

    df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype( str)  # making a unique combination of field and object id
    df = df.drop_duplicates('par_obj', keep='last')
    df['n_filters'] = df['filters'].map(lambda x: len(x.split(',')))

    conditions_from_cosmos = ['mass', 'sfr', 'sSFR', 'photometry']
    if len(set(conditions_from_cosmos).intersection(set(args.plot_conditions))) > 0:
        df = get_crossmatch_with_cosmos(df, args, cosmos_name=cosmos_name)

    return df

# --------------------------------------------------------------------------------------------------------------------
def plot_glass_venn(args, fontsize=10):
    '''
    Plots and saves the Venn diagram for GLASS-NIRISS given a set of conditions
    '''
    args.fontsize = fontsize
    print(f'Plotting Venn diagram for GLASS-NIRISS..')
    ID_col = 'ID_NIRISS'

    # -------loading the data--------------
    #glass_catalog_filename = args.root_dir / 'glass_data' / 'a2744_spec_cat_niriss_20250401.fits'
    glass_catalog_filename = args.root_dir / 'glass_data' / 'full_internal_em_line_data.fits'
    #glass_catalog_filename = args.root_dir / 'glass_data' / 'hlsp_glass-jwst_jwst_niriss_abell2744_v1.0_cat.fits'
    
    tab = Table(fits.open(glass_catalog_filename)[1].data)
    if 'cdf_z' in tab.columns: tab.remove_column('cdf_z') # because multi dimension column does not agree well with pandas
    df = tab.to_pandas()

    set_arr = []
    label_arr = []
    
    # -----------creating the conditions for Venn diagram--------------
    for line in args.line_list:
        df[f'SNR_{line}'] = df[f'flux_{line}'] / df[f'err_{line}']
        condition = df[f'SNR_{line}'] > args.SNR_thresh
        set_arr, label_arr = make_set(df, condition, f'{line} SNR > {args.SNR_thresh}', set_arr, label_arr, colname=ID_col, silent=True)
 
    # --------adding the secure redshift condition------------
    condition = df['Z_FLAG'] == 4
    set_arr, label_arr = make_set(df, condition, f'Secure redshift', set_arr, label_arr, colname=ID_col, silent=True)

    # ----------plot the venn diagrams----------
    cmap = 'plasma'
    dataset_dict = dict(zip(label_arr, set_arr))

    # ---------manually calling draw_venn() so as to modify petal labels (for 0 counts)----------
    petal_labels = generate_petal_labels(dataset_dict.values(), fmt="{size}")
    petal_labels = {logic: value if int(value) > 0 else '' for logic, value in petal_labels.items()}
    colors = generate_colors(cmap=cmap, n_colors=len(label_arr))
    ax = draw_venn(petal_labels=petal_labels, dataset_labels=dataset_dict.keys(), hint_hidden=False, colors=colors, figsize=(8, 6), fontsize=args.fontsize, legend_loc='upper left', ax=None)

    # ----------annotate and save the diagram----------
    fig = ax.figure
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    fig.text(0.9, 0.9, f'GLASS\nTotal {len(df)} objects', c='k', ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize)

    plot_conditions_text = ''
    if np.array(['snr' in item.lower() for item in label_arr]).any(): plot_conditions_text += ','.join(args.line_list) + f',SNR>{args.SNR_thresh}'
    if np.array(['secure' in item.lower() for item in label_arr]).any(): plot_conditions_text += '_zsecure'
    figname = f'GLASSvenn_diagram_{plot_conditions_text}.png'
    save_fig(fig, figname, args)

    # ----------deriving the intersecting dataframe----------
    intersecting_set = set.intersection(*set_arr)
    intersecting_obj = np.transpose(list(intersecting_set))
    df_int = pd.DataFrame({'ID_NIRISS':intersecting_obj})
    df_int = df.merge(df_int, on='ID_NIRISS', how='inner')

    # --------writing the dataframe------------
    plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
    plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}')
    outfilename = get_output_path('glass', args).parent / 'catalogs' / f'glass-a2744_venn_{plot_conditions_text}_df.csv'
    df_int.to_csv(outfilename, index=None)
    print(f'Saved intersecting dataframe at {outfilename}')

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_passage_venn(fields, args, fontsize=10):
    '''
    Plots and saves the Venn diagram of a list of PASSAGE fields given a set of conditions
    '''
    args.fontsize = fontsize
    print(f'Plotting Venn diagram for {fields}..')
    
    # ---------loading full dataframe for all relevant PASSAGE fields------------
    df_passage = load_full_df(fields, args, cosmos_name=cosmos_name)
    df_int, ax = plot_venn(df_passage, args, silent=True)

    # ----------annotate and save the diagram----------
    fig = ax.figure
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    fig.text(0.9, 0.9, f'PASSAGE\nTotal {len(df_passage)} objects', c='k', ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize)

    has_fields = [str(int(item[3:])) for item in pd.unique(df_passage['field'])]
    has_fields.sort(key=natural_keys)
    field_text = ','.join(has_fields)

    plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
    plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')

    figname = f'Par{field_text}_venn_diagram_{plot_conditions_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def get_sfr_df(objlist, args, Zdiag, Zdiag_branch='low', survey='passage', sum=True):
    '''
    Loads and returns a dataframe that holds all the sfr-related properties for a given list of objects
    Returns dataframe
    '''
    df = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})

    filename = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'logOH_sfr_fits.csv'
    
    df_sfr = pd.read_csv(filename)
    df_sfr = df_sfr.drop_duplicates(subset=['field', 'objid', 'Zdiag', 'Zdiag_branch', 'AGN_diag', 'extent_value', 'extent_unit', 'dered_in_NB', 'snr_text', 'only_seg', 'vorbin_text'], keep='last')
    
    if args.re_limit is not None: df_sfr = df_sfr[(df_sfr['extent_unit'] == 're') & (df_sfr['extent_value'] == args.re_limit)]
    else: df_sfr = df_sfr[(df_sfr['extent_unit'] == 'arcsec') & (df_sfr['extent_value'] == args.arcsec_limit)]
    if args.AGN_diag != 'None': df_sfr = df_sfr[df_sfr['AGN_diag'] == args.AGN_diag]
    df_sfr = df_sfr[(df_sfr['Zdiag'] == Zdiag) & (df_sfr['Zdiag_branch'] == Zdiag_branch) & (df_sfr['dered_in_NB'] == args.dered_in_NB) & (df_sfr['only_seg'] == args.only_seg) & (df_sfr['vorbin_text'] == args.vorbin_text)]
    
    df_sfr = df_sfr.drop(['Zdiag', 'Zdiag_branch', 'AGN_diag', 'dered_in_NB', 'snr_text', 'only_seg', 'vorbin_text'], axis=1)
    df_sfr = pd.merge(df, df_sfr, on=['field', 'objid'], how='left')

    cols_to_extract = ['field', 'objid', 'logZ_logSFR_slope', 'logZ_logSFR_slope_u']
    if sum: cols_to_extract += ['SFR_sum', 'SFR_sum_u']
    else: cols_to_extract += ['SFR_int', 'SFR_int_u']
    df_sfr = df_sfr[cols_to_extract]
    df_sfr = df_sfr.rename(columns = {'SFR_sum': 'SFR', 'SFR_sum_u': 'SFR_u', 'SFR_int': 'SFR', 'SFR_int_u': 'SFR_u'})
 
    return df_sfr

# --------------------------------------------------------------------------------------------------------------------
def get_logOH_df(objlist, args, survey='passage'):
    '''
    Loads and returns a dataframe that holds all the metallicity-related properties for a given list of objects
    Returns dataframe
    '''
    df = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})

    filename = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'logOH_allfits.csv'
    
    df_logOH = pd.read_csv(filename)
    df_logOH = df_logOH.drop_duplicates(subset=['field', 'objid', 'Zdiag_branch', 'Zdiag', 'extent_value', 'extent_unit', 'dered_in_NB', 'snr_text', 'only_seg', 'vorbin_text'], keep='last')
    if args.re_limit is not None: df_logOH = df_logOH[(df_logOH['extent_unit'] == 're') & (df_logOH['extent_value'] == args.re_limit)]
    else: df_logOH = df_logOH[(df_logOH['extent_unit'] == 'arcsec') & (df_logOH['extent_value'] == args.arcsec_limit)]    
    if args.AGN_diag != 'None': df_logOH = df_logOH[df_logOH['AGN_diag'] == args.AGN_diag]
    df_logOH = df_logOH[(df_logOH['dered_in_NB'] == args.dered_in_NB) & (df_logOH['only_seg'] == args.only_seg) & (df_logOH['vorbin_text'] == args.vorbin_text)]

    df_logOH = df_logOH.drop(['extent_unit', 'extent_value', 'AGN_diag', 'dered_in_NB', 'snr_text', 'only_seg', 'vorbin_text'], axis=1)
    df_logOH = pd.merge(df, df_logOH, on=['field', 'objid'], how='inner')

    substrings = ['logOH', 'pa', 'q']
    cols_to_join = [item for item in df_logOH.columns if any([i in item for i in substrings])]
    cols_to_join = [item for item in cols_to_join if not any([item.endswith(i) for i in ['_u', '_ulim', '_llim']])]
    Zdiag_arr = np.hstack([[item] if item in ['NB', 'P25'] else [item + '_C25_low', item + '_low', item + '_high'] for item in args.Zdiag])

    for Zdiag in Zdiag_arr:
        if 'low' in Zdiag: df_sub = df_logOH[(df_logOH['Zdiag'] == Zdiag[:-4]) & (df_logOH['Zdiag_branch'] == 'low')]
        elif 'high' in Zdiag: df_sub = df_logOH[(df_logOH['Zdiag'] == Zdiag[:-5]) & (df_logOH['Zdiag_branch'] == 'high')]
        else: df_sub = df_logOH[(df_logOH['Zdiag'] == Zdiag)].drop_duplicates(subset=['field', 'objid', 'Zdiag'], keep='last')
        df_sub = df_sub.drop(['Zdiag', 'Zdiag_branch'], axis=1)
        
        rename_dict = {}
        for col in cols_to_join:
            rename_dict.update({f'{col}': f'{col}_{Zdiag}'})
            rename_dict.update({f'{col}_u': f'{col}_{Zdiag}_u'})
            rename_dict.update({f'{col}_ulim': f'{col}_{Zdiag}_ulim'})
            rename_dict.update({f'{col}_llim': f'{col}_{Zdiag}_llim'})
        df_sub = df_sub.rename(columns=rename_dict)
        #common_columns = list(set(df.columns) & set(df_sub.columns))
        common_columns = ['field', 'objid']
        df = pd.merge(df, df_sub, on=common_columns, how='left')

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_marker_type(field):
    '''
    Returns Path object with full path to the data based on a given field name
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'

    marker_dict = {'passage': 'o', 'glass': 's'}
    return marker_dict[survey]

# --------------------------------------------------------------------------------------------------------------------
def get_passage_masses_from_cosmos(df, id_col='objid', field_col='field', cosmos_name='web'):
    '''
    Derives stellar masses of PASSAGE galaxies present in the given dataframe from COSMOS-Web catalog
    Returns dataframe
    '''
    passage_fields = [item for item in np.unique(df[field_col]) if 'Par' in item]
    
    df_cosmos = pd.DataFrame()
    for index, thisfield in enumerate(passage_fields):
        if 'web' in cosmos_name:
            cosmosfilename = args.input_dir / 'COSMOS' /  f'cosmoswebb_objects_in_{thisfield}.fits'
            df_cosmos_thisfield = read_COSMOSWebb_catalog(filename=cosmosfilename)
            sed_cols_to_extract = ['passage_id', 'LP_mass_minchi2']
        elif '2020' in cosmos_name:
            cosmosfilename = args.input_dir / 'COSMOS' /  f'cosmos2020_objects_in_{thisfield}.fits'
            df_cosmos_thisfield = read_COSMOS2020_catalog(filename=cosmosfilename)
            sed_cols_to_extract = np.hstack([['passage_id'], np.hstack([[f'lp_{item}_best', f'lp_{item}_med_min68', f'lp_{item}_med_max68'] for item in ['mass', 'SFR']])])

        df_cosmos = pd.concat([df_cosmos, df_cosmos_thisfield[sed_cols_to_extract]])
    
    df['passage_id'] = df[field_col].astype(str) + '-' + df[id_col].astype(str)  # making a unique combination of field and object id
    df = pd.merge(df, df_cosmos, on=['passage_id'], how='left')
    df = df.rename(columns={'LP_mass_minchi2':'lp_mass', 'lp_mass_best':'lp_mass', 'lp_SFR_best':'lp_SFR'})
    for item in ['mass', 'SFR']:
        if f'lp_{item}_med_min68' in df:
            df[f'lp_{item}_u'] = df.apply(lambda row: np.mean([row[f'lp_{item}'] - row[f'lp_{item}_med_min68'], row[f'lp_{item}_med_max68'] - row[f'lp_{item}']]))
            df = df.drop(f'lp_{item}_med_min68', axis=1)
            df = df.drop(f'lp_{item}_med_max68', axis=1)
    df = df.drop('passage_id', axis=1)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_masses_from_mysedfit(df, id_col='objid', field_col='field'):
    '''
    Derives stellar masses of PASSAGE and GLASS galaxies present in the given dataframe from COSMOS-Web catalog
    Returns dataframe
    '''
    df_sed = pd.DataFrame()
    
    for thisfield in np.unique(df[field_col]):
        if 'Par' in thisfield: sed_fit_filename = args.output_dir / 'catalogs' / f'{thisfield}_v0.5_venn_OII,NeIII-3867,Hb,OIII,SNR>2.0,mass_df_withSED_for_paper_only_st.csv'
        elif 'glass' in thisfield: sed_fit_filename = args.root_dir / 'glass_data' / f'GLASS_UNCOVER_photometry_for_paper_only_st_withSED_for_paper_only_st.csv'
        df_sed_thisfield = pd.read_csv(sed_fit_filename)
        df_sed_thisfield = df_sed_thisfield.rename(columns={'ID_NIRISS': id_col})
        df_sed_thisfield[field_col] = thisfield
        sed_cols_to_extract = [field_col, id_col, 'log_mass_bgp', 'log_mass_bgp_u', 'log_sfr_bgp', 'log_sfr_bgp_u']
        df_sed = pd.concat([df_sed, df_sed_thisfield[sed_cols_to_extract]])
    
    df = pd.merge(df, df_sed, on=[field_col, id_col], how='left')
    df = df.rename(columns={'log_mass_bgp':'lp_mass', 'log_mass_bgp_u':'lp_mass_u', 'log_sfr_bgp':'lp_SFR', 'log_sfr_bgp_u':'lp_SFR_u'})

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_glass_masses_from_He2024(df, mass_col='lp_mass', id_col='objid', field_col='field'):
    '''
    Derives stellar masses of GLASS galaxies present in the given dataframe from He+2024 catalog
    Returns dataframe
    '''
    glass_watson_catalog_filename = args.root_dir / 'glass_data' / 'a2744_spec_cat_niriss_20250401.fits'
    watson_id_col = 'ID_NIRISS'
    df_glass_watson = Table(fits.open(glass_watson_catalog_filename)[1].data).to_pandas().rename(columns={'RA':'ra', 'DEC':'dec'})
    df_glass_watson = df_glass_watson[df_glass_watson[watson_id_col].isin(df[df[field_col]=='glass-a2744'][id_col])]

    glass_he_catalog_filename = args.root_dir / 'glass_data' / 'He2024_table_indivi.tex'
    he_id_col = 'ID Grism'
    df_glass_he = Table.read(glass_he_catalog_filename).to_pandas().rename(columns={'RA':'ra', 'DEC':'dec', 'log_mass':mass_col})
    
    df_crossmatch = get_crossmatch(df_glass_watson, df_glass_he, sep_threshold=1., df1_idcol=watson_id_col, df2_idcol=he_id_col)
    df_glass_he_crossmatched = df_glass_he[df_glass_he[he_id_col].isin(df_crossmatch['df2_id'])][[he_id_col, mass_col]].reset_index(drop=True)
    df_glass_watson_crossmatched = df_glass_watson[df_glass_watson[watson_id_col].isin(df_crossmatch['df1_id'])][[watson_id_col]].reset_index(drop=True)
    df_crossmatched = pd.concat([df_glass_watson_crossmatched, df_glass_he_crossmatched], axis=1)
    df_crossmatched[[mass_col, f'{mass_col}_u']] = df_crossmatched[mass_col].apply(lambda x: pd.Series(parse_latex_value(x)))

    df = df.set_index(id_col).combine_first(df_crossmatched[[watson_id_col, mass_col]].set_index(watson_id_col)).reset_index().rename(columns={'index': id_col})
    df = df[df.columns].sort_values(by=field_col)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_sfr_from_oii(OII_flux, redshift):
    '''
    Derives SFR given the OII flux and redshift (to convert flux into luminosity) following Figueira+2022
    Returns SFR
    '''
    distance = cosmo.luminosity_distance(redshift)
    OII_lum = OII_flux * 4 * np.pi * (distance.to('cm').value ** 2) # converting to ergs/s (luminosity)
    
    #SFR_OII = OII_lum * 1.26e-41 # luminosity in ergs/s; SFR in Msun/yr; from Vulcani+2010; Salpeter IMF
    #SFR_OII = OII_lum * ufloat(6.58, 1.65) * 10 ** (-42) # luminosity in ergs/s; SFR in Msun/yr; from Kewley+2004 eq 4; agrees well with Selpeter IMF
    SFR_OII = (10 ** ufloat(-39.69, 0.07)) * (OII_lum ** ufloat(0.96, 0.01)) # luminosity in ergs/s; SFR in Msun/yr; from Figueira+2022 Table 6 second-to-last row; Chabrier IMF

    return SFR_OII

# --------------------------------------------------------------------------------------------------------------------
def make_master_df(objlist, args, sum=True):
    '''
    Creates and returns a dataframe that holds all the relevant global properties for a given list of objects
    Returns dataframe
    '''
    sum_text = '_sum' if sum else '_grizli_int'
    filename = args.root_dir / 'zgrad_paper_plots' / f'master_df_using{sum_text}_with_COSMOS{cosmos_name}.csv'
    
    # -------------making the master df-------------------
    if not os.path.exists(filename) or args.clobber:
        print(f'Could not find{filename}, so making a new one..')
        df = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})
    
        # --------get PASSAGE & GLASS stellar masses----------
        #df = get_passage_masses_from_cosmos(df, id_col='objid', field_col='field', cosmos_name=cosmos_name) # from COSMOS-Web/COSMOS2020 catalog
        #df = get_glass_masses_from_He2024(df) # from He+2024 catalog
        
        df = get_masses_from_mysedfit(df, id_col='objid', field_col='field') # from my SED fitting with Bagpipes

        # -------getting redshift and integrated line flux info: looping through objects--------------
        line_list = ['OIII', 'Hb', 'Ha', 'NeIII-3867', 'OII', 'SII']
        new_cols = np.hstack([['RA', 'Dec', 'redshift', 'EB_V'], np.hstack([[f'{line}', f'{line}_u'] for line in line_list])])
        for col in new_cols: df[col] = np.nan # creating provision for the new columns

        for index, row in df.iterrows():
            field = row['field']
            objid = row['objid']
            print(f'\nDoing object {field}:{objid} which is {index + 1} of {len(df)}..')
            full_hdu = load_full_fits(objid, field, args)
            
            args = load_object_specific_args(full_hdu, args, skip_psf=True, field=field)
            df.loc[index, 'RA'] = full_hdu[0].header['RA']
            df.loc[index, 'Dec'] = full_hdu[0].header['Dec']
            df.loc[index, 'redshift'] = args.z
            df.loc[index, 'EB_V'] = args.EB_V.n
            df.loc[index, 'pix_size_arcsec'] = args.pix_size_arcsec
            df.loc[index, 'npix_in_vorbin'] = len(np.ma.compressed(args.voronoi_bin_IDs))
            df.loc[index, 're_kpc'] = args.re_kpc
            df.loc[index, 're_arcsec'] = args.re_arcsec

            args.only_integrated = True
            for line in line_list:
                try: _, _, line_int, line_sum, _ = get_emission_line_map(line, full_hdu, args, silent=True)
                except UnboundLocalError: line_int, line_sum = ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)
                if not sum: line_sum = line_int # choose the grizli reported integrated values
                df.loc[index, f'{line}'] = line_sum.n
                df.loc[index, f'{line}_u'] = line_sum.s  
            args.only_integrated = False

        # --------get metallicity info----------
        df_logOH_passage = get_logOH_df(passage_objlist, args, survey='passage')
        df_logOH_glass = get_logOH_df(glass_objlist, args, survey='glass')
        df_logOH = pd.concat([df_logOH_passage, df_logOH_glass])
        df = pd.merge(df, df_logOH, on=['field', 'objid'])

        # -------getting SFR info--------------
        df_sfr_passage = get_sfr_df(passage_objlist, args, 'NB', survey='passage', sum=sum)
        df_sfr_glass = get_sfr_df(glass_objlist, args, 'NB', survey='glass', sum=sum)
        df_sfr = pd.concat([df_sfr_passage, df_sfr_glass])
        df = pd.merge(df, df_sfr, on=['field', 'objid'])

        # -------writing out dataframe--------------
        df.to_csv(filename, index=False)
        print(f'Saved df as {filename}')
    
    # -------------reading in existing master df-------------------
    else:
        print(f'Reading in existing {filename}')
        df = pd.read_csv(filename)

    # --------selecting only necessary objects------------
    df_base = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})
    df = df_base.merge(df, on=['field', 'objid'], how='left')
    df['marker'] = df['field'].apply(lambda x: get_marker_type(x))
    
    # --------compute OII-based SFRs-----------
    SFR_OII = get_sfr_from_oii(unp.uarray(df['OII'], df['OII_u']), df['redshift'].values)
    df['SFR_OII'] = unp.nominal_values(SFR_OII)
    df['SFR_OII_u'] = unp.std_devs(SFR_OII)

    # ------computing mixing timescales-----------
    #Zdiag = 'NB'
    #Z_SFR_slope = unp.uarray(df['logZ_logSFR_slope'], df['logZ_logSFR_slope_u']) * (10 ** (unp.uarray(df[f'logOH_sum_{Zdiag}'], df[f'logOH_sum_{Zdiag}_u']) - 8.69)) / unp.uarray(df['SFR'], df['SFR_u']) # computing Z-SFR slope from logZ-logSFR slope
    
    # df['Z_SFR_slope'] = unp.nominal_values(Z_SFR_slope) # now this is in yr/Msun
    # df['Z_SFR_slope_u'] = unp.std_devs(Z_SFR_slope) # now this is in yr/Msun

    # df['log_mmol_T20'] = df.apply(lambda row: get_mmol_from_mstar(row['lp_mass'], row['redshift'], method='T20'), axis=1) 
    # for method in ['C18', 'G20', 'B23', 'C23']:
    #     df['log_mgas_' + method] = df.apply(lambda row: get_mg_from_mstar(row['lp_mass'], log_mmol=row['log_mmol_T20'], method=method), axis=1)
    #     t_mix = Z_SFR_slope * (10 ** df['log_mgas_' + method]) / 1e9 # in Gyr
    #     df['t_mix_' + method] = unp.nominal_values(t_mix)  # in Gyr
    #     df['t_mix_' + method + '_u'] = unp.std_devs(t_mix)  # in Gyr

    method = 'my'
    A, alpha = 2.5e-4, 1.5
    beta = 1 - 1/alpha
    B = A ** (1 / alpha)
    slope = unp.uarray(df['logZ_logSFR_slope'], df['logZ_logSFR_slope_u'])
    df['pix_area'] = df['npix_in_vorbin'] * df['pix_size_arcsec'] / df['redshift'].apply(lambda x: cosmo.arcsec_per_kpc_proper(x).value) # converting arcsec to kpc
    Sigma = unp.uarray(df['SFR'] / df['pix_area'], df['SFR_u'] / df['pix_area'])
    t_mix = slope / (B * (Sigma ** beta) * (beta + slope)) # in yr
    #t_mix = slope / (B * (Sigma ** beta) * beta) # in yr
    df['t_mix_' + method] = unp.nominal_values(t_mix)
    df['t_mix_' + method + '_u'] = unp.std_devs(t_mix)
    # log_t_mix = unp.log10(t_mix)
    # df['log_t_mix_' + method] = unp.nominal_values(log_t_mix)
    # df['log_t_mix_' + method + '_u'] = unp.std_devs(log_t_mix)

    # --------computing dynamical time scale-------
    log_t_dyn = unp.log10(2.1e7 * (df['re_kpc'] ** 1.5) * ((10 ** unp.uarray(df['lp_mass'], df['lp_mass_u']) / 1e9) ** -0.5))
    df['log_t_dyn'] = unp.nominal_values(log_t_dyn)
    df['log_t_dyn_u'] = unp.std_devs(log_t_dyn)

    # log_t_ratio = unp.log10(unp.uarray(df['t_mix_' + method], df['t_mix_' + method + '_u'])) - log_t_dyn
    # df['log_t_ratio'] = unp.nominal_values(log_t_ratio)
    # df['log_t_ratio_u'] = unp.std_devs(log_t_ratio)

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_mmol_from_mstar(log_mstar, redshift, method='T20'):
    '''
    Compute log molecular gas mass, given log stellar mass and a scaling relation method
    Returns log molecular gas mass
    '''
    
    if method == 'T20':
        log_mol_frac = 0.06 - 3.33 * (np.log10(1 + redshift) - 0.65) ** 2 - 0.41 * (log_mstar - 10.7) # eq 9 of Chowdhury+2023, eq 10 of Bera+2023b, this is actually from Tacconi+2020
        log_mmol = log_mstar - log_mol_frac
    
    return log_mmol

# --------------------------------------------------------------------------------------------------------------------
def get_mg_from_mstar(log_mstar, log_mmol=-np.inf, method='C18'):
    '''
    Compute log gas mass, given log stellar mass and a scaling relation method
    Returns log gas mass
    '''
    if method == 'C18':
        gas_frac = get_C18_scaling(log_mstar)
        log_mgas = log_mstar + np.log10(gas_frac)
    elif method == 'G20':
        gas_frac = 0.158 * log_mstar ** 2 - 3.548 * log_mstar + 19.964 # eq 17 of Gullieuszik+2020
        log_mgas = log_mstar + np.log10(gas_frac)
    elif method == 'B23':
        log_mhi = 8.977 + 0.183 * (log_mstar - 9) # eq 10 of Bera+2023a, eq 7 of Bera+2023b
        log_mgas = np.log10(10 ** log_mhi + 10 ** log_mmol)
    elif method == 'C23':
        log_mhi = 0.38 * np.log10(10 ** log_mstar / 1e10) + 10.564 # eq 11 of Chowdhury+2023 (z~1.3)
        log_mgas = np.log10(10 ** log_mhi + 10 ** log_mmol)

    return log_mgas

# --------------------------------------------------------------------------------------------------------------------
def get_C18_scaling(log_mstar):
    '''
    Gas fraction-mstar scaling relation from Catinella+2018
    Values copied from Table 2 of the paper
    Returns interpolated gas fraction at a given log_mstar
    '''
    log_gf_arr = [0.148, 0.040, -0.511, -0.485, -0.785, -0.965, -1.238, -1.496]
    log_mstar_arr = [9.16, 9.44, 9.75, 10.05, 10.34, 10.65, 10.95, 11.21]
    interpolation = interp1d(log_mstar_arr, log_gf_arr, fill_value='extrapolate')

    log_gf = interpolation(log_mstar)
    gas_frac = 10 ** log_gf

    return gas_frac

# --------------------------------------------------------------------------------------------------------------------
def plot_gasmass_comparison(df, args, mass_col='lp_mass', fontsize=10):
    '''
    Makes and saves plot for comparison of gas masses derived using different methods
    '''
    args.fontsize = fontsize
    col_arr = ['g', 'b', 'r', 'salmon']
    methods = [item.split('_')[-1] for item in df.columns if 'log_mgas' in item]

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.95)

    for index, method in enumerate(methods):
     for index2, m in enumerate(pd.unique(df['marker'])):
         df_sub = df[df['marker'] == m]
         ax.scatter(df_sub[mass_col], df_sub['log_mgas_' + method], c=col_arr[index], marker=m, plotnonfinite=True, s =100, lw=1, edgecolor='k', label=method if index2==0 else None)
    
    ax.plot(log_mass_lim, log_mass_lim, ls='dashed', c='k', lw=0.5)

    # ---------annotate axes and save figure-------
    ax.legend(fontsize=args.fontsize)
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log M$_{\rm gas}$/M$_{\odot}$', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    ax.set_ylim(7, 10.5)

    figname = f'log_mgas_comparison.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS(df, args, mass_col='lp_mass', sfr_col='lp_SFR', fontsize=10):
    '''
    Plots and saves the SF-main sequence diagram given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting SFMS...')
    if 'log_' in sfr_col and sfr_col not in df and sfr_col[4:] in df:
        df = break_column_into_uncertainty(df, sfr_col[4:], make_log=True)

    sfr_col2 = 'log_SFR_OII' # None # 
    if sfr_col2 is not None and 'log_' in sfr_col2 and sfr_col2 not in df and sfr_col2[4:] in df:
        df = break_column_into_uncertainty(df, sfr_col2[4:], make_log=True)
        to_replace = df[sfr_col] == df[sfr_col2]
        df.loc[df[sfr_col] == df[sfr_col2], [sfr_col, sfr_col + '_u']] = [np.nan, np.nan]

    df = df.sort_values(by='redshift')
    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.1, top=0.95)

    # ----------plotting----------
    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub[sfr_col], c=df_sub['redshift'], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis', label='GLASS (This work)' if m == 's' else 'PASSAGE (This work)')
        if sfr_col2 is not None: ax.scatter(df_sub[mass_col], df_sub[sfr_col2], marker=m, s=100, lw=1, edgecolor='k', facecolors='none')
    if sfr_col + '_u' in df: ax.errorbar(df[mass_col], df[sfr_col], yerr=df[sfr_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[sfr_col], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if sfr_col2 is not None:
        ax.errorbar(df[mass_col], df[sfr_col2], yerr=df[sfr_col2 + '_u'], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
        #for i, row in df.iterrows(): ax.arrow(row[mass_col], row[sfr_col], 0, row[sfr_col2] - row[sfr_col], lw=0.5, color='grey', length_includes_head=True, head_width=3e-2)
        #ax.vlines(df[mass_col], df[sfr_col], df[sfr_col2], color='k', ls='dashed', lw=0.5, zorder=-10)
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row[mass_col], row[sfr_col], f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')

    # -------fitting-----------------
    if args.fit_correlation:
        df = df.sort_values(by=mass_col)
        linefit_odr = odr_fit(df, quant_x=mass_col, quant_y=sfr_col)
        ax = plot_fitted_line(ax, linefit_odr, df[mass_col], 'salmon', args, quant='', short_label=True, index=0)

    # ----------making colorbar----------
    cbar = plt.colorbar(p, pad=0.01)
    #cbar.set_label(r'$\log$ O III/H$\beta$', fontsize=args.fontsize)
    cbar.set_label(f'Redshift', fontsize=args.fontsize)
    cbar.set_ticklabels([f'{item:.1f}' for item in cbar.get_ticks()], fontsize=args.fontsize)

    # ----------plotting theoretical diagrams----------
    #ax = plot_SFMS_Whitaker14(ax, 2, color='yellowgreen')
    ax = plot_SFMS_Shivaei15(ax, color='salmon')
    ax = plot_SFMS_Popesso23(ax, 2, color='cornflowerblue')
    #ax = plot_SFMS_Popesso23(ax, 3, color='royalblue')

    # ---------annotate axes and save figure-------
    plt.legend(fontsize=args.fontsize, loc='lower right')
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log SFR (M$_{\odot}$/yr)', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    ax.set_ylim(-0.5, 1.2 if sfr_col2 is None else 1.5)

    figname = f'SFMS_colorby_redshift.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_MEx(df, args, mass_col='lp_mass', fontsize=10):
    '''
    Plots and saves the mass-excitation diagram given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting mass-excitation diagram...')
    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.95)

    quant =unp.log10(unp.uarray(df['OIII'], df['OIII_u']) / unp.uarray(df['Hb'], df['Hb_u']))
    df['O3Hb'] = unp.nominal_values(quant)
    df['O3Hb_u'] = unp.std_devs(quant)

    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub['O3Hb'], c=df_sub['redshift'], marker=m, s=100, edgecolor='k', lw=1, cmap='viridis', vmin=1.7, vmax=3.1, label='GLASS (This work)' if m == 's' else 'PASSAGE (This work)')
    ax.errorbar(df[mass_col], df['O3Hb'], yerr=df['O3Hb_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df['O3Hb'], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row[mass_col], row['O3Hb'], f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')

    # ----------making colorbar----------
    cbar = plt.colorbar(p, pad=0.01)
    cbar.set_label('Redshift', fontsize=args.fontsize)
    cbar.set_ticklabels([f'{item:.1f}' for item in cbar.get_ticks()], fontsize=args.fontsize)

    # ---------annotate axes and save figure-------
    plt.legend(fontsize=args.fontsize)
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log O III/H$\beta$', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    ax.set_ylim(-1, 1.7)

    # ---------adding literature lines from Juneau+2014 (https://iopscience.iop.org/article/10.1088/0004-637X/788/1/88/pdf)----------
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    y_up = np.piecewise(x, [x <= 10, x > 10], [lambda x: (0.375 / (x - 10.5)) + 1.14, lambda x: np.poly1d([410.24, -109.333, 9.71731, -0.288244][::-1])(x)]) # J14 eq 1
    y_lo = np.piecewise(x, [x <= 9.6, x > 9.6], [lambda x: (0.375 / (x - 10.5)) + 1.14, lambda x: np.poly1d([352.066, -93.8249, 8.32651, -0.246416][::-1])(x)]) # J14 eq 2
    ax.plot(x, y_up, c='k', ls='dashed', lw=2, label='Juneau+2014')
    ax.plot(x, y_lo, c='brown', ls='dashed', lw=2)
    plt.legend()

    figname = f'mass_excitation_colorby_redshift.png'
    save_fig(fig, figname, args)

    return

# -------------------------------------------------------------------------------------------------------------------
def plot_filled_region(df, xcol, ycol, ax, color=None, noscatter=False, zorder=None, alpha=0.5, label=None):
    '''
    Function to overplot a filled region between multiple y-values for each x-value
    If there is only 1 y-value for each x-value then this will lead to a simple line plot
    This is only efficient/useful when only a few unique x-values exist in the data
    :return: axis handle
    '''    
    xarr = pd.unique(df[xcol])
    if len(xarr) == len(df): # all values of x are unique; then just to do a thick line plot
        p = ax.plot(df[xcol], df[ycol], lw=2, color=color, zorder=zorder, alpha=alpha, label=label)
        if not noscatter: ax.scatter(df['redshift'], df['Zgrad'], c=color, lw=1.0, s=50, ec='k', zorder=zorder + 1 if zorder is not None else None)

    else: # some values of x have multiple y-values
        yup_arr, ylow_arr = [], []
        for thisx in xarr:
            yup_arr.append(df[df[xcol] == thisx][ycol].max())
            ylow_arr.append(df[df[xcol] == thisx][ycol].min())

        p = ax.fill_between(xarr, ylow_arr, yup_arr, color=color, alpha=alpha, edgecolor='k', lw=0.1, zorder=zorder, label=label)
        if not noscatter: ax.scatter(df['redshift'], df['Zgrad'], c=color, lw=0.5, s=10, ec='k', zorder=zorder + 1 if zorder is not None else None)

    return ax, p

# --------------------------------------------------------------------------------------------------------------------
def extract_middle(text, prefix, suffix):
    '''
    Build regex pattern with optional prefix/suffix
    Returns string
    '''
    start = text.find(prefix)
    if start != -1:
        start += len(prefix)
    else:
        start = 0  # no prefix, start from beginning

    end = text.find(suffix, start)
    if end == -1:
        end = len(text)  # no suffix, go till end

    return text[start:end]

# --------------------------------------------------------------------------------------------------------------------
def plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col='logOH_slope_NB', fontsize=10):
    '''
    Plots and saves the mass-metallicity gradient plot given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting MZgrad...')
    zgrad_col = np.atleast_1d(zgrad_col)
    lim_dict = defaultdict(lambda: [None, None], redshift=[1.7, 3.1], lp_mass=[7.0, 11.2])
    literature_dir = args.root_dir / 'zgrad_paper_plots' / 'literature'
    
    # ----------setting up the diagram----------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6)) # one for dex/re and another for dex/kpc
    fig.subplots_adjust(left=0.072, right=0.99, bottom=0.1, top=0.78, wspace=0.2)
    label_dict = smart_dict({'redshift': 'Redshift', 'logOH_sum_NB': r'$\log$ (O/H) + 12 [NB]', 'lp_mass': r'log M$_*$/M$_{\odot}$'})

    # ----------plotting----------
    df = df.sort_values(by=mass_col)
    color_arr = ['darkkhaki', 'teal', 'yellowgreen']
    this_work = []
    for index, this_zgrad_col in enumerate(zgrad_col):
        Zdiag = extract_middle(this_zgrad_col, prefix='logOH_slope_mcmc_', suffix='_low')
        if Zdiag != 'NB':
            if '_' in Zdiag: Zdiag = '; '.join(Zdiag.split('_'))
            else: Zdiag += '; C20'
        
        for m in pd.unique(df['marker']):
            df_sub = df[df['marker'] == m]
            Zdiag_text = f'This work; {Zdiag}'
            if index == 0: this_label = f'{Zdiag_text} (GLASS)' if m == 's' else f'{Zdiag_text} (PASSAGE)'
            elif m == 'o': this_label = Zdiag_text
            else: this_label = None

            axes[0].scatter(df_sub[mass_col], df_sub[this_zgrad_col], c=color_arr[index], marker=m, plotnonfinite=True, s=150 - index * 50, lw=1, edgecolor='k', cmap='viridis', label=this_label, zorder=200)
            h = axes[1].scatter(df_sub[mass_col], df_sub[this_zgrad_col] / df_sub['re_kpc'], c=color_arr[index], marker=m, plotnonfinite=True, s=150 - index * 50, lw=1, edgecolor='k', cmap='viridis', label=this_label, zorder=200)
            this_work.append(h)

        if this_zgrad_col + '_u' in df:
            axes[0].errorbar(df[mass_col], df[this_zgrad_col], yerr=df[this_zgrad_col + '_u'], c=color_arr[index], fmt='none', lw=1, alpha=0.5, zorder=100)
            axes[1].errorbar(df[mass_col], df[this_zgrad_col] / df['re_kpc'], yerr=df[this_zgrad_col + '_u'] / df['re_kpc'], c=color_arr[index], fmt='none', lw=1, alpha=0.5, zorder=100)

        if mass_col + '_u' in df:
            axes[0].errorbar(df[mass_col], df[this_zgrad_col], xerr=df[mass_col + '_u'], c=color_arr[index], fmt='none', lw=1, alpha=0.5, zorder=100)
            axes[1].errorbar(df[mass_col], df[this_zgrad_col] / df['re_kpc'], xerr=df[mass_col + '_u'], c=color_arr[index], fmt='none', lw=1, alpha=0.5, zorder=100)

    for ax in axes: ax.axhline(0, ls='--', c='k', lw=0.5)

    # --------plotting FOGGIE filled region: for dex/kpc----------
    xcol = 'log_mass'
    ycol = 'Zgrad'
    foggie_filename = literature_dir / 'FOGGIE_allhalos.csv'
    df_foggie = pd.read_csv(foggie_filename)
    if xcol != 'redshift': df_foggie = df_foggie[(df_foggie['log_mass'].between(8.5, 11.5)) & ~(df_foggie['halo'] == 8508)]
    df_foggie = df_foggie.sort_values(by=xcol)
    #ax.plot(df_foggie[xcol], df_foggie[ycol], c='sienna', lw=0.1, zorder=-2, alpha=0.5)

    xarr = np.linspace(df_foggie[xcol].min(), df_foggie[xcol].max(), 1000) # uniform grid of redshift values
    xarr = xarr[:-1] + np.diff(xarr) / 2.
    new_df = pd.DataFrame()
    for thishalo in pd.unique(df_foggie['halo']):
        thisdf = df_foggie[df_foggie['halo'] == thishalo]
        thisnewdf = pd.DataFrame()
        thisnewdf[ycol] = np.interp(xarr, thisdf[xcol], thisdf[ycol]) # interpolating data from each halo onto the uniform grid
        thisnewdf[xcol] = xarr
        new_df = pd.concat([new_df, thisnewdf])    
    
    axes[1], foggie = plot_filled_region(new_df, xcol, ycol, axes[1], color='salmon', noscatter=True, label='FOGGIE (Acharyya+25)')

    # ---------for plotting other observed data from literature------------
    legend_dict = {'sami': 'SAMI', 'manga': 'MaNGA', 'califa': 'CALIFA', 'sharda_scaling1': 'S21 scaling 1', 'sharda_scaling2': 'S21 scaling 2', 'mingozzi2020_izi': 'Mingozzi+20 (IZI)', 'wang17': 'Wang+17', 'jones15': 'Jones+15', 'venturi24': 'Venturi+24', 'li25': 'Li+25', 'ju25': 'Ju+25'}
    marker_dict = {'sami': '+', 'manga': 'x', 'califa': '1', 'sharda_scaling1': 'v', 'sharda_scaling2': '^', 'mingozzi2020_izi': '<', 'wang17': 'v', 'jones15': '^', 'venturi24': '>', 'li25': 'D', 'ju25': 'd'}
    ls_dict = {'sami': 'dotted', 'manga': 'dotted', 'califa': 'dotted', 'sharda_scaling1': 'solid', 'sharda_scaling2': 'dashed', 'mingozzi2020_izi': 'dotted', 'wang17': 'dotted', 'jones15': 'dotted', 'venturi24': 'dotted', 'li25': 'dotted', 'ju25': 'dotted'}
    color_dict = {'sami': 'firebrick', 'manga': 'chocolate', 'califa': 'darkgoldenrod', 'sharda_scaling1': 'k', 'sharda_scaling2': 'k', 'mingozzi2020_izi': 'peru', 'wang17': 'brown', 'jones15': 'sandybrown', 'venturi24': 'bisque', 'li25': 'grey', 'ju25': 'lightskyblue'}
 
   # --------plotting Sharda+21 data: for both dex/re and dex/kpc----------
    s21 = []
    for ax, search_text in zip(axes, ['mzgr_*re.csv', 'mzgr_*kpc.csv']):
        literature_files = glob.glob(str(literature_dir / search_text))
        literature_files.sort(key=natural_keys)

        for index, this_file in enumerate(literature_files):
            sample = '_'.join(Path(this_file).stem.split('_')[1:-1])
            df_lit = pd.read_csv(this_file, names=['log_mass', 'Zgrad'], sep=', ')
            if 'scaling' in sample: h = ax.plot(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=1, ls=ls_dict[sample], label=legend_dict[sample])[0]
            else: h = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], ec='k', s=50, lw=2, marker=marker_dict[sample], label=legend_dict[sample])
            if 're' in search_text: s21.append(h)

    # --------plotting Mingozzi+2020 MaNGA data: for dex/re----------
    sample = 'mingozzi2020_izi'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', names=['log_mass', 'Zgrad'], sep=', ')
    m20 = axes[0].scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])

     # --------plotting GLASS-HST data: for dex/re----------
    sample = 'jones15'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv')
    j15 = axes[0].scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    axes[0].errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=0.5, fmt='none')
    
    sample = 'wang17'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv')
    w17 = axes[0].scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    axes[0].errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=0.5, fmt='none')

    # --------plotting Wang+19 data: for dex/kpc----------
    w19 = axes[1].scatter([9.05, 9.41], [0.12, 0.11], color='sienna', marker='*', s=100, lw=0.5, edgecolor='k', label='Wang+19')
    axes[1].errorbar([9.05, 9.41], [0.12, 0.11], xerr=[0.05, 0.01], yerr=[0.01, 0.02], c='gray', fmt='none', lw=1, alpha=0.5)

    # --------plotting Wang+22 data: for dex/kpc----------
    w22 = axes[1].scatter(9.04, 0.165, color='goldenrod', marker='*', s=200, lw=0.5, edgecolor='k', label='Wang+22')
    axes[1].errorbar(9.04, 0.165, xerr=0.05, yerr=0.023, c='gray', fmt='none', lw=1, alpha=0.5)

   # --------NOT plotting Li+25 (NGDEEP) redshift-binnd data: dex/kpc----------
    sample = 'li25'
    redshift_bins = [0.228323913043,0.83149495798,1.27080989010,1.7380614592,2.20353918995,2.72550200322,3.24097725652]
    redshift_bin_err = [0.15068650322,0.11935359080,0.1508593450,0.150354201085,0.135542672091,0.131702882241,0.136858778221]
    zgrad_bin = [-0.03853139,-0.03113908,-0.00728418,0.0219643,0.01319803,-0.03173848,-0.03036552]
    zgrad_bin_err = [0.00374489,0.00333496,0.00166307,0.00397386,0.00409635,0.01259178,0.01507544]

    # --------plotting Li+25 (NGDEEP) data: for dex/kpc vs mass----------
    df_lit = pd.read_csv(literature_dir / 'zgrad_Li25.dat', delim_whitespace=True)

    df_lit_sub = df_lit[df_lit['id'].str.contains('stack')]
    l25s = axes[1].scatter(df_lit_sub['logMass'], df_lit_sub['zgrad'], color=color_dict[sample], lw=1, label=legend_dict[sample] + ' (6 < z < 7; stacked)', marker='h', ec='k', s=100)
    axes[1].errorbar(df_lit_sub['logMass'], df_lit_sub['zgrad'], yerr=df_lit_sub['zgrad_err'], color=color_dict[sample], lw=0.5, fmt='none')

    df_lit_sub = df_lit[~(df_lit['id'].str.contains('stack'))]
    l25 = axes[1].scatter(df_lit_sub['logMass'], df_lit_sub['zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample] + ' (z < 3.5)', marker='h')
    axes[1].errorbar(df_lit_sub['logMass'], df_lit_sub['zgrad'], yerr=df_lit_sub['zgrad_err'], color=color_dict[sample], lw=0.5, fmt='none')

    # --------plotting Li+25 (ASPIRE + FRESCO) mass-binned data: for dex/kpc vs mass----------
    mass_bin = [7.74, 8.44, 9.06, 7.67, 8.14]
    mass_bin_err = [0.21, 0.25, 0.31, 0.19, 0.4]
    zgrad_bin = [-0.42,-0.22,0,-0.62,-0.30]
    zgrad_bin_err = [0.17,0.16,0.18,0.22,0.22]

    l25b = axes[1].scatter(mass_bin, zgrad_bin, color=color_dict[sample], lw=1, label=legend_dict[sample] + ' (z > 5; binned)', marker='D', ec='k', s=100)
    axes[1].errorbar(mass_bin, zgrad_bin, xerr=mass_bin_err, yerr=zgrad_bin_err, color=color_dict[sample], lw=0.5, fmt='none')

    # --------plotting Venturi+24 data: for dex/re and dex/kpc----------
    sample = 'venturi24'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', comment='#', delim_whitespace=True)
    v24 = axes[0].scatter(df_lit['log_mass'], df_lit['Zgrad_re'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], marker=marker_dict[sample], ec='k')
    axes[0].errorbar(df_lit['log_mass'], df_lit['Zgrad_re'], xerr=df_lit['log_mass_u'], yerr=df_lit['Zgrad_re_u'], color=color_dict[sample], lw=0.5, fmt='none')

    axes[1].scatter(df_lit['log_mass'], df_lit['Zgrad_kpc'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], marker=marker_dict[sample], ec='k')
    axes[1].errorbar(df_lit['log_mass'], df_lit['Zgrad_kpc'], xerr=df_lit['log_mass_u'], yerr=df_lit['Zgrad_kpc_u'], color=color_dict[sample], lw=0.5, fmt='none')

    # --------plotting Ju+25 data: for dex/re and dex/kpc----------
    sample = 'ju25'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', comment='#', delim_whitespace=True)
    j25 = axes[0].scatter(df_lit['log_mass'], df_lit['Zgrad_re'], color=color_dict[sample], lw=0.5, label=legend_dict[sample] + ' (z~1)', marker=marker_dict[sample], ec='k')
    axes[0].errorbar(df_lit['log_mass'], df_lit['Zgrad_re'], yerr=df_lit['Zgrad_re_u'], color=color_dict[sample], lw=0.5, fmt='none')

    axes[1].scatter(df_lit['log_mass'], df_lit['Zgrad_kpc'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], marker=marker_dict[sample], ec='k')
    axes[1].errorbar(df_lit['log_mass'], df_lit['Zgrad_kpc'], yerr=df_lit['Zgrad_kpc_u'], color=color_dict[sample], lw=0.5, fmt='none')

    # --------plotting Franchetto+21 data: for dex/kpc----------
    coeff = [-0.199, 0.199 * 10 - 0.432] # Franchetto+21 eq 6
    xarr = log_mass_lim
    f21 = axes[1].plot(xarr, np.poly1d(coeff)(xarr), color='cornflowerblue', ls='dotted', label=f'Franchetto+21')[0]

    # --------plotting Sun+25 FIRE2 data: for dex/kpc----------
    xarr = [8.299, 10.006, 10.743]
    yarr = [-0.03115, -0.01311, -0.03508] # these data grabbed from the red line in Fig 6 left panel of Sun+25 using WebPlotDigititizer
    fire2 = axes[1].plot(xarr, yarr, color='maroon', ls='solid', label=f'FIRE-2 (Sun+25)')[0]

    # ---------annotate axes and save figure-------
    handles = this_work + [j15, w17, w19, m20] + s21 + [w22, v24, l25, l25s, l25b, j25] + [f21, foggie, fire2]
    labels = [h.get_label() for h in handles]
    fig.legend(handles, labels, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.99), fontsize=args.fontsize / args.fontfactor)

    for ax in axes:
        #ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')
        ax.set_xlabel(label_dict[mass_col], fontsize=args.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=args.fontsize)
        ax.set_xlim(lim_dict[mass_col][0], lim_dict[mass_col][1])
    axes[0].set_ylabel(r'log $\nabla$Z$_r$ (dex/R$_e$)', fontsize=args.fontsize)
    axes[1].set_ylabel(r'log $\nabla$Z$_r$ (dex/kpc)', fontsize=args.fontsize)
    #axes[0].set_ylim(-0.5, 1)

    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    xcol_text = 'z' if mass_col == 'redshift' else 'M' if 'mass' in mass_col else mass_col
    figname = f'{xcol_text}Zgrad_Zdiag_{",".join(zgrad_col)}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_MZsfr(df, args, mass_col='lp_mass', zgrad_col='logZ_logSFR_slope', fontsize=10):
    '''
    Plots and saves the mass vs metallicity-SFR slope given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting MZ-SFR...')

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.95)

    # ----------plotting----------
    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub[zgrad_col], c=df_sub['logOH_sum_NB'], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis', vmin=7.5, vmax=9.1, label='GLASS' if m == 's' else 'PASSAGE')
    if zgrad_col + '_u' in df: ax.errorbar(df[mass_col], df[zgrad_col], yerr=df[zgrad_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[zgrad_col], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row[mass_col], row[zgrad_col], f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')
    
    ax.axhline(0, ls='--', c='k' if not args.fortalk else 'w', lw=0.5)

    # -------fitting-----------------
    if args.fit_correlation:
        df = df.sort_values(by=mass_col)
        #linefit = odr_fit(df, quant_x=mass_col, quant_y=zgrad_col)
        #ax = plot_fitted_line(ax, linefit, df[mass_col], 'salmon', args, quant='', short_label=True, index=0)
        df = df.dropna(subset = [mass_col, zgrad_col], axis=0)
        linefit = np.polyfit(df[mass_col], df[zgrad_col], 1)
        xarr = np.linspace(ax.get_xlim()[0] * 1.001, ax.get_xlim()[1] * 0.999, 10)
        y_fitted = np.poly1d(linefit)(xarr)
        ax.plot(xarr, y_fitted, color='salmon', lw=1, ls='dashed')
        ax.text(0.1, 0.05, f' Slope = {linefit[0]: .2f}', c='salmon', fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)

    # ----------making colorbar----------
    cbar = plt.colorbar(p, pad=0.01)
    cbar.set_label(r'$\log$ (O/H) + 12 [NB]', fontsize=args.fontsize)
    cbar.set_ticklabels([f'{item:.1f}' for item in cbar.get_ticks()], fontsize=args.fontsize)

    # ---------annotate axes and save figure-------
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'$\log$ (O/H) + 12 - $\log \Sigma_{\rm SFR}$ slope', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    ax.set_ylim(-3, 3)

    ax.legend(fontsize=args.fontsize, loc='upper left')

    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'MZsfr_colorby_Z_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_MZR(df, args, mass_col='lp_mass', z_col='logOH_sum_NB', colorcol='logOH_slope_mcmc_NB', fontsize=10):
    '''
    Plots and saves the mass vs metallicity given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting MZR...')

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.95)

    # ----------plotting----------
    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub[z_col], c=df_sub[colorcol], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis', label='GLASS (This work)' if m == 's' else 'PASSAGE (This work)')
    if z_col + '_u' in df: ax.errorbar(df[mass_col], df[z_col], yerr=df[z_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[z_col], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row[mass_col], row[z_col], f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')
    
    # -------fitting-----------------
    if args.fit_correlation:
        df = df.sort_values(by=mass_col)
        linefit = odr_fit(df, quant_x=mass_col, quant_y=z_col)
        ax = plot_fitted_line(ax, linefit, df[mass_col], 'salmon', args, quant='', short_label=True, index=0)

    # --------plotting literature relations----------
    #ax = plot_MZR_literature(ax)
    def zahid_func(log_mass, Z0, M0, gamma):
        return Z0 - np.log10(1 + (10**log_mass/10**M0)**(-gamma))
    
    zahid_data_13 = pd.DataFrame({'Sample':['SHELS', 'DEEP2', 'Y12', 'E06'], \
                                    'Redshift':[0.29, 0.78, 1.24, 2.26], \
                                    'Z0':[9.130, 9.161, 9.06, 9.06], \
                                    'Z0_u':[0.007, 0.026, 0.36, 0.27], \
                                    'Z0_fl':[0, 0, 0, 0], \
                                    'M0':[9.304, 9.661, 9.6, 9.7], \
                                    'gamma':[0.77, 0.65, 0.7, 0.6]}) # these are based on KK04 frame and need to be converted to PPN2
    #zahid_data_13['Z0_PPN2'], zahid_data_13['Z0_PPN2_u'], _ = np.vstack(np.array(zahid_data_13.apply(lambda row: convert(row, diagnostic='Z0_KK04', coeff=[-1.3188000, 35.051680, -309.54480, 916.7484], llim=8.2, ulim=9.2), axis=1))).transpose() # K08 table3 col 2 last set of rows
    zahid_data_13['linestyle'] = 'dashed' # because these metallicities have been converted

    zahid_data_14 = pd.DataFrame({'Sample':['SDSS', 'COSMOS'], \
                                    'Redshift':[0.08, 1.55], \
                                    'Z0':[8.710, 8.740], \
                                    'Z0_u':[0.001, 0.042], \
                                    'M0':[8.76, 9.93], \
                                    'gamma':[0.66, 0.88]}) # these are based on PPN2 frame already and need not be coverted
    zahid_data_14['linestyle'] = 'solid' # because these metallicities have NOT been converted

    zahid_data = pd.concat([zahid_data_13, zahid_data_14], join='inner', ignore_index=True).reset_index(drop=True)
    zahid_data = zahid_data.sort_values(by='Redshift').reset_index(drop=True)

    col_ar = ['orange', 'black', 'cyan', 'seagreen', 'blue', 'salmon', 'gray']
    for i in range(len(zahid_data)):
        xarr = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)
        ax.plot(xarr, zahid_func(xarr, zahid_data['Z0'][i], zahid_data['M0'][i], zahid_data['gamma'][i]), color=col_ar[i], lw=2, ls=zahid_data['linestyle'][i], label= 'z = '+str(zahid_data['Redshift'][i]) + '; ' + zahid_data['Sample'][i], zorder=-5)

    # ----------making colorbar----------
    cbar = plt.colorbar(p, pad=0.01)
    cbar.set_label(r'log $\nabla$Z$_r$ (dex/kpc)' if args.re_limit is None else r'log $\nabla$Z$_r$ (dex/R$_e$)', fontsize=args.fontsize)
    cbar.set_ticklabels([f'{item:.1f}' for item in cbar.get_ticks()], fontsize=args.fontsize)

    # ---------annotate axes and save figure-------
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'$\log$ (O/H) + 12 [NB]', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    ax.set_ylim(7.0, 9.1)

    ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'MZR_colorby_{colorcol}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_Mtmix(df, args, mass_col='lp_mass', ycol='t_mix', fontsize=10, mgas_method='G20', colorcol='logZ_logSFR_slope'):
    '''
    Plots and saves the mass vs mixing timescale given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    print(f'Plotting M-t_mix...')
    color_label_dict = {'logZ_logSFR_slope':r'$\log$ Z-$\log \Sigma_{\rm SFR}$ slope', 'SFR': r'SFR (M$_{\odot}$/yr)'}

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.14, right=0.99, bottom=0.1, top=0.95)

    if mgas_method is None:
        mgas_methods = [item.split('_')[-1] for item in df.columns if 'log_mgas' in item]
        col_arr = ['teal', 'cornflowerblue', 'crimson', 'salmon']
    else:
        mgas_methods = [mgas_method]

    # ----------plotting----------
    for index, method in enumerate(mgas_methods):
        if 'log_' in f'{ycol}_{method}' and f'{ycol}_{method}' not in df and f'{ycol}_{method}'[4:] in df:
            df = break_column_into_uncertainty(df, f'{ycol}_{method}'[4:], make_log=True)

        for index2, m in enumerate(pd.unique(df['marker'])):
            df_sub = df[df['marker'] == m]
            p = ax.scatter(df_sub[mass_col], df_sub[f'{ycol}_{method}'], c=df_sub[colorcol] if mgas_method is not None else col_arr[index], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis', label='GLASS' if m == 's' else 'PASSAGE')
            #p = ax.scatter(df_sub[mass_col], df_sub[f'log_t_ratio'], c=df_sub[colorcol] if mgas_method is not None else col_arr[index], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis', label='GLASS' if m == 's' else 'PASSAGE')
            #ax.scatter(df_sub[mass_col], df_sub['log_t_dyn'], c=df_sub[colorcol] if mgas_method is not None else col_arr[index], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', cmap='viridis')
        if f'{ycol}_{method}_u' in df: ax.errorbar(df[mass_col], df[f'{ycol}_{method}'], yerr=df[f'{ycol}_{method}_u'], c='gray', fmt='none', lw=1, alpha=0.5)
        if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[f'{ycol}_{method}'], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
        #if f'log_t_dyn_u' in df: ax.errorbar(df[mass_col], df['log_t_dyn'], yerr=df['log_t_dyn_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    
    ax.axhline(0, c='k', ls='--', lw=0.5)
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row[mass_col], row[f'{ycol}_{method}'], f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')

    # -------fitting-----------------
    if args.fit_correlation:
        df = df.sort_values(by=mass_col)
        linefit_odr = odr_fit(df, quant_x=mass_col, quant_y=f'{ycol}_{method}')
        ax = plot_fitted_line(ax, linefit_odr, df[mass_col], 'salmon', args, quant='', short_label=True, index=0)

    # ----------making colorbar----------
    if mgas_method is not None:
        cbar = plt.colorbar(p, pad=0.01)
        cbar.set_label(color_label_dict[colorcol], fontsize=args.fontsize)
        cbar.set_ticklabels([f'{item:.1f}' for item in cbar.get_ticks()], fontsize=args.fontsize)
    else:
        ax.legend(fontsize=args.fontsize)

    # ---------annotate axes and save figure-------
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log (Metal mixing timescale $t_{mix}$ (yr))' if 'log_' in ycol else r'Metal mixing timescale $t_{mix}$ (yr)', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    ax.set_xlim(log_mass_lim[0], log_mass_lim[1])
    if mgas_method is None or mgas_method == 'C23': ax.set_ylim(-0.7, 4.3)
    elif mgas_method is not None and mgas_method == 'my':
        if 'log_' in ycol: ax.set_ylim(2.2, 3.7) # (None, None) #
        else: ax.set_ylim(-100, 1100) #(None, None) #
    else: ax.set_ylim(-0.5, 0.6)

    ax.legend(fontsize=args.fontsize, loc='upper right')

    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'M_tmix_colorby_{colorcol}_{mgas_method}_upto_{extent_text}.png' if mgas_method is not None else f'M_tmix_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_re_histogram(df, args, fontsize=10):
    '''
    Plots and saves the histogram of R_e (in arcseconds) given a dataframe with list of objects and properties
    '''
    args.fontsize = fontsize
    PSF_FWHM_arcsec_dict = {'F115W': 0.094, 'F150W': 0.092, 'F200W': 0.106} # value of PSF FWHM derived from Table 1 of https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-performance/niriss-point-spread-functions#gsc.tab=0

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)

    # ----------plotting----------
    ax.hist(df['re_arcsec'], histtype='step', lw=2, color='cornflowerblue')
    
    # -----over-plotting F150W PSF radius------------
    color = 'sienna'
    filter = 'F150W'
    ax.axvline(PSF_FWHM_arcsec_dict[filter] / 2, ls='dashed', lw=2, color=color)
    ax.text(1.1 * PSF_FWHM_arcsec_dict[filter] / 2, ax.get_ylim()[0] + 0.5 * np.diff(ax.get_ylim())[0], r'0.5 $\times$ ' + f'FWHM ({filter})', rotation=90, va='center', ha='left', fontsize=args.fontsize, color=color)

    # ---------annotate axes and save figure-------
    ax.set_xlabel(r'R$_e$ (arcseconds)', fontsize=args.fontsize)
    ax.set_ylabel(r'Frequency', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    figname = f'Re_histogram.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def make_latex_table(df, args, Zdiag='NB', sum=True):
    '''
    Makes and saves a latex table for with metallicity and other columns, given a dataframe with list of objects and properties
    Returns the latex table as pandas dataframe
    '''
    column_dict = {'field':'Field', 'objid':'ID', 'redshift':r'$z$', 'lp_mass':r'$\log$ M$_{\star}$/M$_{\odot}$', 'lp_SFR':r'$\log$ SFR (M$_{\odot}$/yr)', 'SFR':r'SFR (M$_{\odot}$/yr)', 'logOH_int_NB':'$\log$ O/H + 12$_{total}$', 'logOH_sum_NB':'$\log$ O/H + 12$_{total}$', 'logOH_slope_NB':r'$\nabla Z$ (dex/kpc)' if args.re_limit is None else r'$\nabla Z$ (dex/R$_e$)', 'logOH_slope_mcmc_NB':r'$\nabla Z$ (dex/kpc)' if args.re_limit is None else r'$\nabla Z$ (dex/R$_e$)'}
    decimal_dict = defaultdict(lambda: 2, redshift=1, lp_mass=1, RA=4, Dec=4)
    base_cols = ['field', 'objid', 'RA', 'Dec', 'redshift', 'lp_mass', 'SFR']
    cols_with_errors = [f'logOH_sum_{Zdiag}' if sum else f'logOH_int_{Zdiag}', f'logOH_slope_mcmc_{Zdiag}', 'lp_mass', 'SFR']

    tex_df = pd.DataFrame()
    for thiscol in base_cols + cols_with_errors:
        if thiscol in cols_with_errors: # special treatment for columns with +/-
            tex_df[thiscol] = df.apply(lambda row: r'$' + f'{row[thiscol]: .{decimal_dict[thiscol]}f} \pm {row[thiscol + "_u"]: .{decimal_dict[thiscol]}f}' + r'$' if np.isfinite(row[thiscol]) else '-', axis=1)
        elif thiscol in ['field', 'objid']:
            tex_df[thiscol] = df[thiscol]
        else:
            tex_df[thiscol] = df[thiscol].map(lambda x:  r'$' + f'{x: .{decimal_dict[thiscol]}f}' + r'$' if np.isfinite(x) else '-')

    tex_df = tex_df.rename(columns=column_dict) # change column names to nice ones
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    tabname = args.root_dir / 'zgrad_paper_plots' / f'paper_table_upto_{extent_text}.tex'
    tex_df.to_latex(tabname, index=False, escape=False)
    print('Saved latex table at', tabname)

    return tex_df

# --------------------------------------------------------------------------------------------------------------------
def get_photoionisation_model_grid(args):
    '''
    Loads and returns a given photoionisation model grid of ratios
    Returns pandas dataframe
    '''
    geom_path_dict = {'s':['spherical', 'sp'], 'p':['plane_par', 'pp']} # to determine directory structures based on geometry and iso parameters
    if args.phot_models.lower() in ['mappings', 'map']:
        grid_filename = args.mappings_dir / 'grids' / f'mappings_grid_{geom_path_dict[args.geometry][1]}_iso_{args.iso}.txt'
        df = pd.read_table(grid_filename, delim_whitespace=True)

    elif args.phot_models.lower() in ['nebulabayes', 'nb']:
        grid_filename = Path(NebulaBayes.__path__[0]) / 'grids' / 'NB_HII_grid.fits.gz'
        df = Table(fits.getdata(grid_filename)).to_pandas()
        df['log q'] = np.round(df['log U'] + np.log10(3e10), 1)

        quant_names_dict = {'Z':'12 + log O/H', 'log(q)':'log q', 'log(P/k)':'log P/k', 'log(U)':'log U'}

        quant_names = [quant_names_dict[quant] for quant in [args.quantity1, args.quantity2, args.quantity3]]
        df = df.sort_values(by=quant_names)
        df = df.rename(columns={v: k for k, v in quant_names_dict.items()})

    print(f'Reading in existing grid from {grid_filename}')

    return df

# --------------------------------------------------------------------------------------------------------------------
def remove_lambda_from_label(label):
    '''
    Remove the wavelength info from a given line label
    Return the rest of the label
    '''
    str_to_find = '$\lambda'
    if str_to_find in label: return label[ : label.find(str_to_find) - 1]
    else: return label

# --------------------------------------------------------------------------------------------------------------------
def get_line_labels(lines, omit_lambda=False):
    '''
    Returns a list of nice labels including latex math mode for a givenlist line name
    '''
    label_dict = smart_dict({'OIII': r'[O III] $\lambda$5007', 'OIII5007': r'[O III] $\lambda$5007', \
                            'NeIII-3867': r'[Ne III]  $\lambda$3869', 'NeIII': r'[Ne III] $\lambda$3869', 'NeIII3869': r'[Ne III] $\lambda$3869', \
                            'OII': r'[O II] $\lambda$3727,29', 'OII3726_29': r'[O II] $\lambda$3727,29', 'OII3727_29': r'[O II] $\lambda$3727,29', \
                            'Hbeta': r'H$\beta$', 'Hb': r'H$\beta$', \
                            'Halpha': r'H$\alpha$', 'Ha': r'H$\alpha$'})
    lines = np.atleast_1d(lines)
    labels = ['+'.join([remove_lambda_from_label(label_dict[item]) if omit_lambda else label_dict[item] for item in line.split(',')]) for line in lines]

    if len(lines) == 1: labels = labels[0]
    return labels

# --------------------------------------------------------------------------------------------------------------------
def get_ratio_labels(ratio_name, omit_lambda=False):
    '''
    Returns a nice label including latex math mode for a given ratio name
    '''
    lines = ratio_name.split('/')
    label = '/'.join(get_line_labels(lines, omit_lambda=omit_lambda))

    return label

# --------------------------------------------------------------------------------------------------------------------
def plot_photoionisation_model_grid(ratio_x, ratio_y, args, fit_y_envelope=False, fontsize=10):
    '''
    Plots and saves the ratio vs ratio parameter space and, optionally, a fit to its envelope, for a given photoionisation model grid
    '''
    args.fontsize = fontsize
    print(f'Plotting photoionisation grid for {ratio_y} vs {ratio_x}..')
    if args.phot_models.lower() in ['mappings', 'map']: line_label_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,9', 'NeIII':'[NeIII]3869', 'NeIII-3867':'[NeIII]3869', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31'}) # to map between user input line labels and line labels used in ratio_list.txt file
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})
    args.xnum_line, args.xden_line = ratio_x.split('/')
    args.ynum_line, args.yden_line = ratio_y.split('/')
    df = get_photoionisation_model_grid(args)

    ratio_x = f'{line_label_dict[args.xnum_line]}/{line_label_dict[args.xden_line]}'
    ratio_y = f'{line_label_dict[args.ynum_line]}/{line_label_dict[args.yden_line]}'
    if ratio_x not in df: df[ratio_x] = df[line_label_dict[args.xnum_line]] / df[line_label_dict[args.xden_line]]
    if ratio_y not in df: df[ratio_y] = df[line_label_dict[args.ynum_line]] / df[line_label_dict[args.yden_line]]

    # ------declare the figure-----------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.1, top=0.95, wspace=0.2)

    # --------plot the model ratios---------------
    ax, xratio_name, yratio_name = plot_ratio_grid(df, ax, args, color1='seagreen', color2='cornflowerblue', color3='slategray')

    # ------annotate figure-----------------------
    #ax.grid(which='both', color='gray', linestyle='solid', linewidth=1, alpha=0.3)
    ax.set_xlabel('Log ' + get_ratio_labels(xratio_name), fontsize=args.fontsize)
    ax.set_ylabel(f'Log {get_ratio_labels(yratio_name)}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    df = df[(df[xratio_name] > 0) & (df[yratio_name] > 0)] # to avoid math errors later while taking log
    xmin = args.xmin if args.xmin is not None else np.log10(np.min(df[xratio_name]) * 0.9)
    xmax = args.xmax if args.xmax is not None else np.log10(np.max(df[xratio_name]) * 1.1)
    ymin = args.ymin if args.ymin is not None else np.log10(np.min(df[yratio_name]) * 0.9)
    ymax = args.ymax if args.ymax is not None else np.log10(np.max(df[yratio_name]) * 1.1)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # ------for fitting to upper y-axis envelope of data------------
    if fit_y_envelope:
        nbins = 25

        def func(x, *popt): return np.poly1d(popt)(x)
        p_init = [1, 1, 1, 1]

        xarr = np.linspace(np.log10(np.min(df[xratio_name])), np.log10(np.max(df[xratio_name])), nbins)
        df['bin'] = pd.cut(np.log10(df[xratio_name]), bins=xarr)
        grouped = df.groupby('bin')
        xbins = np.log10(grouped[xratio_name].mean().values)
        ybins = np.log10(grouped[yratio_name].max().values)
        #ax.plot(xbins, ybins, c='r', lw=0.5)

        color = 'brown'
        popt, pcov = curve_fit(func, xbins, ybins, p0=p_init)
        ax.plot(xarr, func(xarr, *popt), c=color, lw=2)
        ax.text(0.98, 0.01, f'Fit coefficients = [{",".join([f"{item:.2f}" for item in popt])}]', c=color, ha='right', va='bottom', fontsize=args.fontsize, transform=ax.transAxes)

    # ---------saving figure--------------
    if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: phot_model_text = 'NebulaBayes'
    figname = f'{phot_model_text}_grid_{yratio_name.replace("/", "-")}_vs_{xratio_name.replace("/", "-")}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_photoionisation_models(ratio_y, parameter_x, args, fontsize=10):
    '''
    Plots and saves the ratio vs parameter for a given photoionisation model grid and parameter
    '''
    args.fontsize = fontsize
    print(f'Plotting photoionisation models for {ratio_y} vs {parameter_x}..')
    if args.phot_models.lower() in ['mappings', 'map']: line_label_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,9', 'NeIII':'[NeIII]3869', 'NeIII-3867':'[NeIII]3869', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31'}) # to map between user input line labels and line labels used in ratio_list.txt file
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})
    args.ynum_line, args.yden_line = ratio_y.split('/')
    df = get_photoionisation_model_grid(args)

    ratio_y = f'{line_label_dict[args.ynum_line]}/{line_label_dict[args.yden_line]}'
    if ratio_y not in df: df[ratio_y] = df[line_label_dict[args.ynum_line]] / df[line_label_dict[args.yden_line]]

    # ------declare the figure-----------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95, wspace=0.2)

    # --------plot the model ratios---------------
    ax, ratio_name = plot_ratio_model(df, ax, args)

    # ------annotate figure-----------------------
    ax.grid(which='both', color='gray', linestyle='solid', linewidth=1, alpha=0.3)
    ax.legend(fontsize=args.fontsize/2)
    ax.set_xlabel('log(O/H) + 12' if parameter_x == 'Z' else parameter_x, fontsize=args.fontsize)
    ax.set_ylabel(f'Log {get_ratio_labels(ratio_name)}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)


    xmin = args.xmin if args.xmin is not None else np.min(np.unique(df[parameter_x])) * 0.99
    xmax = args.xmax if args.xmax is not None else  np.max(np.unique(df[parameter_x])) * 1.01
    ymin = args.ymin if args.ymin is not None else np.log10(np.min(df[ratio_name]) * 0.9)
    ymax = args.ymax if args.ymax is not None else np.log10(np.max(df[ratio_name]) * 1.1)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # ---------saving figure--------------
    if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: phot_model_text = 'NebulaBayes'
    figname = f'{phot_model_text}_model_{ratio_name.replace("/", "-")}_vs_{parameter_x}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def load_full_fits(objid, field, args):
    '''
    Loads the full.fits file for a given object and returns the hdulist
    '''
    # ------determining directories---------
    product_dir = get_data_path(field, args)
    full_fits_file = product_dir / 'full' / f'{field}_{objid:05d}.full.fits'
    maps_fits_file = product_dir / 'maps' / f'{field}_{objid:05d}.maps.fits'

    if os.path.exists(maps_fits_file): # if the maps.fits files are available
        full_hdu = fits.open(maps_fits_file)

    elif os.path.exists(full_fits_file):  # if the full.fits files are available
        full_hdu = fits.open(full_fits_file)

    return full_hdu

# --------------------------------------------------------------------------------------------------------------------
def load_1d_fits(objid, field, args):
    '''
    Loads the 1D.fits file for a given object and returns the hdulist
    '''
    product_dir = get_data_path(field, args)
    od_filename = product_dir / 'spec1D' / f'{field}_{objid:05d}.1D.fits'
    if not os.path.exists(od_filename): od_filename = Path(str(od_filename).replace('.1D.', '.spec1D.'))
    od_hdu = fits.open(od_filename)

    return od_hdu

# --------------------------------------------------------------------------------------------------------------------
def load_object_specific_args(full_hdu, args, skip_vorbin=False, field=None, sum=True, skip_psf=False):
    '''
    Loads some object specific details into args
    Returns modified args
    '''
    # ----------loading object specific things into args---------------------
    args.z = full_hdu[0].header['REDSHIFT']
    args.id = full_hdu[0].header['ID']
    args.field = field
    line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
    args.pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
    args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
    args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines])  # replace 'OIII-5007' with 'OIII'

    # ----------getting info from the photometric catalog-------------
    if args.use_elliptical_bins:
        if 'glass' in field: survey_name, args.drv = 'glass', 'orig'
        elif 'Par' in field: survey_name, args.drv = 'passage', 'v0.5'
        catalog_file = args.root_dir / f'{survey_name}_data/' / args.drv / args.field / 'Products' / f'{args.field}_photcat.fits'
        catalog = GTable.read(catalog_file)
        
        obj = catalog[catalog['id']==args.id][0]
        args.mag = obj['mag_auto']        
        args.q = obj['b_image'] / obj['a_image']
        args.pa = obj['theta_image'] # radians

    # -------determining true center of object-------------
    args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W')

    # --------determining effective radius of object---------------------
    args.re_kpc, args.re_arcsec = get_re(full_hdu, args, filter='F150W')
    if args.re_limit is not None:
        args.pix_size_re = args.pix_size_arcsec / args.re_arcsec
        offset = args.pix_size_re / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
        args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)
    else:
        offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
        args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)

    # ---------------segmentation map---------------
    segmentation_map = full_hdu['SEG'].data
    segmentation_map = np.roll(segmentation_map, args.ndelta_xpix, axis=0)
    segmentation_map = np.roll(segmentation_map, args.ndelta_ypix, axis=1)
    args.segmentation_map_arcsec = trim_image(segmentation_map, args, skip_re_trim=True)
    args.segmentation_map = trim_image(segmentation_map, args)
    
    # -----------PSF array-----------
    if not skip_psf: args.psf_array = make_psf('F150W', args, supersampling_factor=1)

    ####################################
    if field is not None and 'glass' in field and args.voronoi_snr == 4:
        voronoi_snr = float(2)
        print(f'\nWARNING: Actually voronoi snr={voronoi_snr} for {field}-{args.id}') ##
    else:
        voronoi_snr = float(args.voronoi_snr)
    ####################################

    # ---------------voronoi binning stuff---------------
    if args.vorbin and args.voronoi_line is not None and not skip_vorbin:
        args.voronoi_bin_IDs = get_voronoi_bin_IDs(full_hdu, voronoi_snr, plot=False, quiet=True, args=args)
        args.voronoi_bin_distances = get_voronoi_bin_distances(full_hdu, 'OIII', args)

    # ---------------radially binning stuff---------------
    if args.radbin and args.voronoi_line is not None and not args.only_integrated:
        args.voronoi_bin_IDs = get_radial_bin_IDs(full_hdu, snr_thresh=args.voronoi_snr, plot=False, quiet=True, args=args)
        args.voronoi_bin_distances = get_voronoi_bin_distances(full_hdu, 'OIII', args)
        args.vorbin = True

    # ---------------dust value---------------
    try:
        _, EB_V_int, EB_V_sum = get_EB_V(full_hdu, args, verbose=False, silent=True)
        if sum: args.EB_V = EB_V_sum
        else: args.EB_V = EB_V_int
    except:
        args.EB_V = ufloat(0, 0)

    # ----------global stellar mass----------------------
    if field is not None:
        if 'Par' in field: sed_fit_filename = args.output_dir / 'catalogs' / f'{field}_v0.5_venn_OII,NeIII-3867,Hb,OIII,SNR>2.0,mass_df_withSED_for_paper_only_st.csv'
        elif 'glass' in field: sed_fit_filename = args.root_dir / 'glass_data' / f'GLASS_UNCOVER_photometry_for_paper_only_st_withSED_for_paper_only_st.csv'
        df_sed = pd.read_csv(sed_fit_filename)
        df_sed = df_sed.rename(columns={'ID_NIRISS': 'objid'})
        try: args.log_mass = df_sed[df_sed['objid'] == args.id]['log_mass_bgp'].values[0]
        except IndexError: args.log_mass = -99

    return args

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile_old(image, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, short_label=False, quant='logOH'):
    '''
    Plots and fits the radial profile from a given 2D image in a given axis
    Returns the axis handle and the linefit
    '''
    label_dict = smart_dict({'SFR': r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', 'logOH': r'$\log$ (O/H) + 12', 'Z': r'$\log$ (O/H) + 12'})

    # ----making the dataframe before radial profile plot--------------
    df = pd.DataFrame({'quant': unp.nominal_values(np.ma.compressed(image)), 'quant_u': unp.std_devs(np.ma.compressed(image))})

    # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
    # ----------getting the distance map--------
    distance_map = get_distance_map(np.shape(image), args)
    df['distance'] = np.ma.compressed(np.ma.masked_where(image.mask, distance_map.data))
    
    if args.vorbin:
        df['bin_ID'] = np.ma.compressed(np.ma.masked_where(image.mask, args.voronoi_bin_IDs.data))
        df = df.groupby('bin_ID', as_index=False).agg(np.mean)

    if args.radius_max is not None: df = df[df['distance'] <= args.radius_max]
    df = df.sort_values(by='distance').reset_index(drop=True)

    # -------plotting--------
    ax.scatter(df['distance'], df['quant'], c='grey', s=20, alpha=1)
    ax.errorbar(df['distance'], df['quant'], yerr=df['quant_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)
    ax.set_aspect('auto') 

    # -------radial fitting-------------
    try:
        fit_color = 'salmon'
        linefit, linecov = np.polyfit(df['distance'], df['quant'], 1, cov=True, w=1. / (df['quant_u']) ** 2 if (df['quant_u'] > 0).any() else None)
        y_fitted = np.poly1d(linefit)(df['distance'])
        ax.plot(df['distance'], y_fitted, color=fit_color, lw=1, ls='dashed')
        linefit = np.array([ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))])
        if quant in ['logOH', 'Z']:
            label = r'$\nabla$Z$_r$' + f' = {linefit[0].n: .2f}' if short_label else r'$\nabla$Z$_r$' + f' = {linefit[0]: .2f} dex/kpc'
        else:
            label = r'$\nabla$' + f'{quant}' + r'$_r$' + f' = {linefit[0].s: .2f}' if short_label else r'$\nabla$' + f'{quant}' + r'$_r$' + f' = {linefit[0]: .2f} dex/kpc'
        ax.text(0.1, 0.05, label, c=fit_color, fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)
    except:
        print(f'WARNING: Could not fit radial profile, returning nan fit parameters')
        linefit = np.array([ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)])

    # --------annotating axis--------------
    if xlim is not None: ax.set_xlim(xlim[0], xlim[1]) # kpc
    if ylim is not None: ax.set_ylim(ylim[0], ylim[1])
    if not skip_annotate: ax = annotate_axes(ax, 'Radius (kpc)', label_dict[quant], args, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)
    
    return ax, linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_fitted_line(ax, linefit, xarr, fit_color, args, quant='log_OH', short_label=False, index=0, label=None):
    '''
    Computes and plots the fitted line given the linear fit parameters, on an axisting axis handle
    Returns axis handle
    '''
    xarr = np.linspace(np.min(xarr), np.max(xarr), 10)
    if quant in ['logOH', 'Z', 'log_OH']: main_text = r'$\nabla$Z$_r$'
    else: main_text = r'$\nabla$'
    unit_text = 'dex/kpc' if args.re_limit is None else r'dex/R$_e$'

    if len(np.shape(linefit)) == 1:
        y_fitted = np.poly1d(unp.nominal_values(linefit))(xarr)
        y_low = np.poly1d(unp.nominal_values(linefit) - unp.std_devs(linefit))(xarr)
        y_up = np.poly1d(unp.nominal_values(linefit) + unp.std_devs(linefit))(xarr)
        value_text = f'{linefit[0]: .2f}'
    
    else:
        y_fitted = np.poly1d(linefit[0])(xarr)
        y_low = np.poly1d(linefit[1])(xarr)
        y_up = np.poly1d(linefit[2])(xarr)
        value_text = f'{linefit[0][0]: .2f}' if short_label else r'%.2f$^{%.2f}_{%.2f}$' % (linefit[0][0], linefit[0][0] - linefit[1][0], linefit[2][0] - linefit[0][0])
    
    ax.plot(xarr, y_fitted, color=fit_color, lw=1, ls='dashed')
    ax.fill_between(xarr, y_low, y_up, color=fit_color, lw=0, alpha=0.3)

    if label is None:
        label = main_text + ' = ' + value_text
        if not short_label: label += unit_text
    ax.text(0.1, 0.05 + index * 0.1, label, c=fit_color, fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)
    
    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(df, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, short_label=False, quant='log_OH', Zdiag='NB', do_mcmc=False):
    '''
    Plots and fits the radial profile from a given 2D image in a given axis
    Returns the axis handle and the linefit
    '''
    label_dict = smart_dict({'SFR': r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', 'logOH': r'$\log$ (O/H) + 12', 'Z': r'$\log$ (O/H) + 12', 'log_OH': r'$\log$ (O/H) + 12'})

    if args.re_limit is not None: df['distance'] /= args.re_kpc # converting from kpc to Re
    if args.radius_max is not None: df = df[df['distance'] <= args.radius_max]
    df = df.sort_values(by='distance').reset_index(drop=True)

    # -------fitting by various methods--------------
    quant_x = 'distance'
    
    linefit_original = original_fit(df, quant_x=quant_x, quant_y=quant)
    linefit_odr = odr_fit(df, quant_x=quant_x, quant_y=quant)
    if do_mcmc:
        # run lenstronomy
        linefit_lenstronomy = lenstronomy_fit_wrap(df, args, filter='F150W', supersampling_factor=1, Zdiag=Zdiag, quant_x=quant_x, quant_y=quant, return_intermediate=False)
        #linefit_lenstronomy = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]
        
        # run MCMC
        params_llim, params_median, params_ulim = mcmc_vorbin_fit(df, args, filter='F150W', quant_x=quant_x, quant_y=quant, plot_corner=False, Zdiag=Zdiag)
        
        # convert the unit of best fit slope in order to plot/save
        params_median = convert_slope_unit(params_median, args, quant_x=quant_x) # convert from dex/arcsecond to something else
        params_llim = convert_slope_unit(params_llim, args, quant_x=quant_x) # convert from dex/arcsecond to something else
        params_ulim = convert_slope_unit(params_ulim, args, quant_x=quant_x) # convert from dex/arcsecond to something else
    else:
        linefit_lenstronomy = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]
        nparams = 4
        params_median = np.ones(nparams) * np.nan
        params_llim = np.ones(nparams) * np.nan
        params_ulim = np.ones(nparams) * np.nan
   
    # -------plotting the data and the fits--------
    ax.scatter(df[quant_x], df[quant], c='grey', s=20, alpha=1)
    if quant + '_u' in df: ax.errorbar(df['distance'], df[quant], yerr=df[quant + '_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

    xarr = df[quant_x]
    if do_mcmc:
        #ax = plot_fitted_line(ax, linefit_original, xarr, 'salmon', args, quant=quant, short_label=short_label, index=0)
        #ax = plot_fitted_line(ax, linefit_odr, xarr, 'sienna', args, quant=quant, short_label=short_label, index=1)
        #ax = plot_fitted_line(ax, linefit_lenstronomy, xarr, 'green', args, quant=quant, short_label=short_label, index=2)
        #ax = plot_fitted_line(ax, linefit_mcmc, xarr, 'cornflowerblue', args, quant=quant, short_label=short_label, index=3)
        ax = plot_fitted_line(ax, [params_median[:2], params_llim[:2], params_ulim[:2]], xarr, 'salmon', args, quant=quant, short_label=short_label, index=0)
    else:
        #ax = plot_fitted_line(ax, linefit_original, xarr, 'salmon', args, quant=quant, short_label=short_label, index=0)
        #ax = plot_fitted_line(ax, linefit_odr, xarr, 'sienna', args, quant=quant, short_label=short_label, index=1)
        ax = plot_fitted_line(ax, linefit_odr, xarr, 'salmon', args, quant=quant, short_label=short_label, index=0)

    ax.set_aspect('auto') 

    # --------annotating axis--------------
    if xlim is not None: ax.set_xlim(xlim[0], xlim[1]) # kpc
    if ylim is not None: ax.set_ylim(ylim[0], ylim[1])
    if not skip_annotate: ax = annotate_axes(ax, 'Radius (kpc)' if args.re_limit is None else r'Radius (R$_e$)', label_dict[quant], args, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)
    
    if quant in ['logOH', 'Z', 'log_OH']: return ax, [linefit_original, linefit_odr, linefit_lenstronomy, params_llim, params_median, params_ulim]
    else: return ax, [linefit_original, linefit_odr]

# --------------------------------------------------------------------------------------------------------------------
def plot_linelist(ax, fontsize=10, line_list_file=None, show_log_flux=False):
    '''
    Plots a list of emission line wavelengths on the given axis
    Returns axis handle
    '''
    lines_df = get_linelist(wave_lim=ax.get_xlim(), line_list_file=line_list_file)
    previous_wave = np.min(lines_df['restwave']) * 0.1 # very small value, so that it is far left (blue-ward) of the actual plot
    last_flipped = True # flip switch for determining if last label was placed to the left or right of the vertical line

    for index in range(len(lines_df)):
        this_wave = lines_df.iloc[index]['restwave']
        ax.axvline(this_wave, c='cornflowerblue', lw=1, alpha=0.5)
        #if this_wave - previous_wave < 500 or last_flipped:
        if last_flipped:
            xpos = this_wave - np.diff(ax.get_xlim())[0] * 0.04
        else:
            xpos = this_wave + np.diff(ax.get_xlim())[0] * 0.01
        last_flipped = not last_flipped
        ypos = ax.get_ylim()[1] * 1.0006 if show_log_flux else ax.get_ylim()[1] * 0.98
        ax.text(xpos, ypos, lines_df.iloc[index]['LineID'].strip(), rotation=90, va='top', ha='left', fontsize=fontsize)
        previous_wave = this_wave

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_1D_spectra(od_hdu, ax, args, show_log_flux=False):
    '''
    Plots the 1D spectra for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting 1D spectra..')
    nfilters = sum(['GRISM' in item for item in list(od_hdu[0].header.keys())])
    filters = [od_hdu[0].header[f'GRISM{item + 1:03d}'] for item in range(nfilters)]
    color = 'orangered'
    norm_factor = 1e-19

    # -------plot 1D spectra for each filter-----------
    for index, filter in enumerate(filters):
        print(f'Plotting 1D spectra for filter {filter} which is {index+1} of {nfilters}..')
        table = Table(od_hdu[filter].data).to_pandas()
        table = table[table['wave'].between(table['wave'].min() * (1 + args.trim_filter_by_wavelength_factor), table['wave'].max() * (1 - args.trim_filter_by_wavelength_factor))]

        table['rest_wave'] = table['wave'] / (1 + args.z)
        table['norm_flux'] = table['flux'] / table['flat'] / norm_factor # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        table['norm_flux_u'] = table['err'] / table['flat'] / norm_factor # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        table['norm_cont'] = table['cont'] / table['flat'] / norm_factor # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        #table['norm_flux_u'] *= 10

        if show_log_flux:
            table['norm_cont'] = np.log10(table['norm_cont']) + np.log10(norm_factor)
            log_quant = unp.log10(unp.uarray(table['norm_flux'] * norm_factor, table['norm_flux_u'] * norm_factor))
            table['norm_flux'] = unp.nominal_values(log_quant)
            table['norm_flux_u'] = unp.std_devs(log_quant)

        ax.step(table['rest_wave'], table['norm_flux'], lw=1, c=color, alpha=1, where='pre')
        ax.plot(table['rest_wave'], table['norm_cont'], lw=1, c='grey', alpha=1)
        ax.fill_between(table['rest_wave'], table['norm_cont'] - table['norm_flux_u']/2, table['norm_cont'] + table['norm_flux_u']/2, lw=0, color=color, alpha=0.5, step='pre')#, drawstyle='steps')

    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    if show_log_flux: ylabel = r'$\log$ f$_{\lambda}$ ergs/s/cm$^2$/A)'
    else: ylabel = r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)'
    ax.set_ylabel(ylabel, fontsize=args.fontsize)
    if not show_log_flux: ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune='both'))
    ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    ax = plot_linelist(ax, fontsize=args.fontsize / args.fontfactor, line_list_file=HOME / 'Work/astro/Mappings/labframe.passagelinelist', show_log_flux=show_log_flux)

    # --------decorating axis----------
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def annotate_axes(ax, xlabel, ylabel, args, label='', clabel='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''
    ax.text(0.05, 0.9, label, c='k', fontsize=args.fontsize/args.fontfactor, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9), transform=ax.transAxes)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(xlabel, fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(ylabel, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    if not hide_cbar and p is not None:
        cax = inset_axes(ax, width="5%", height="100%", loc='right', bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.set_label(clabel, fontsize=args.fontsize)

        cbar.locator = ticker.MaxNLocator(integer=cticks_integer, nbins=4)#, prune='both')
        cbar.update_ticks()
        if hide_cbar_ticks:
            cbar.ax.set_yticklabels([])
        else:
            cbar.ax.tick_params(labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def annotate_kpc_scale_bar(kpc, ax, args, label=None, color='k', loc='lower left'):
    '''
    Annotate existing axis with a scale bar corresponding to a given kpc length
    Returns axis handle
    '''
    pix = kpc * cosmo.arcsec_per_kpc_proper(args.z).value  # converting kpc to arcsec
    if args.re_limit is not None: pix /= args.re_arcsec # converting arcsec to Re units
    scalebar = AnchoredSizeBar(ax.transData, pix, label, loc, pad=0.5, color=color, frameon=False, size_vertical=0.01, fontproperties={'size':args.fontsize / args.fontfactor})
    ax.add_artist(scalebar)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_direct_image(full_hdu, filter, args):
    '''
    Loads the direct image for a given filter for a given object
    Returns the image
    '''
    try:
        hdu = full_hdu['DSCI', filter.upper()]
        image = hdu.data
        exptime = 1
    except:
        try:
            hdu = full_hdu['DSCI', f'{filter}-CLEAR']
            image = hdu.data
            exptime = full_hdu[0].header[f'T_{filter.upper()}']
        except:
            try:
                hdu = full_hdu['DSCI', f'{filter.upper()}-{filter.upper()}-CLEAR']
                image = hdu.data
                exptime = full_hdu[0].header[f'T_{filter.upper()}']
            except:
                full_field_filename = get_data_path(args.field, args) / f'{args.field}_{filter.lower()}-clear_drz_sci.fits'
                print(f'{filter.upper()} not found in full_hdu extension. Therefore trying to get cutout from full field image {full_field_filename}')
            
                exptime = fits.open(full_field_filename)[0].header['EXPTIME']
                pos = SkyCoord(full_hdu[0].header['RA'], full_hdu[0].header['DEC'], unit = 'deg')
                size = 2 * args.arcsec_limit * u.arcsec
                target_header = full_hdu['DSCI', 'F140W'].header
                
                temp1, temp2 = args.only_seg, args.vorbin
                args.only_seg, args.vorbin = False, False
                image = get_cutout(full_field_filename, pos, size, target_header, args, plot_test_axes=None, skip_re_trim=True)
                args.only_seg, args.vorbin = temp1, temp2

    image = np.roll(image, args.ndelta_xpix, axis=0)
    image = np.roll(image, args.ndelta_ypix, axis=1)

    image = trim_image(image, args, skip_re_trim=True)
 
    return image, exptime

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_image(full_hdu, filter, ax, args, cmap='Greens'):
    '''
    Plots the direct image for a given filter for a given object, in a given axis
    Returns the axis handle
    '''
    image, _ = get_direct_image(full_hdu, filter, args)
    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, alpha=1)
    ax.set_aspect('auto')    
    
    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_rgb_image(full_hdu, filters, ax, args, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, hide_filter_names=False, hide_pa=False):
    '''
    Plots the direct image as an RGB image, combining the given list of filters, for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting the RGB images with filters {filters}..')
    PSF_FWHM_arcsec_dict = {'F115W': 0.094, 'F150W': 0.092, 'F200W': 0.106} # value of PSF FWHM derived from Table 1 of https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-performance/niriss-point-spread-functions#gsc.tab=0

    # -------get image for each filter---------
    image_arr, exptime_arr = [], []
    col_arr = ['r', 'lightgreen', 'cornflowerblue']
    for index, filter in enumerate(filters):
        image, exptime = get_direct_image(full_hdu, filter, args)
        image_arr.append(image)
        exptime_arr.append(exptime)
        if not hide_filter_names:
            ax.text(0.05, 0.95 - index * 0.1, f'{filter}', c=col_arr[index], fontsize=args.fontsize / args.fontfactor, ha='left', va='top', transform=ax.transAxes, zorder=10)
            if filter == 'F150W': ax.add_patch(plt.Circle((0.35, 0.93 - index * 0.1), PSF_FWHM_arcsec_dict[filter]/2, color=col_arr[index], fill=True, lw=0.5, transform=ax.transAxes)) # for over-plotting PSF size

    # -------normalising each image to exptime---------
    for index in range(len(exptime_arr)):
        image_arr[index] = image_arr[index] * exptime_arr[0] / exptime_arr[index]

    # -------create RGB image---------
    pctl, maximum = 99.9, 0.
    for img in image_arr:
        val = np.percentile(img, pctl)
        if val > maximum: maximum = val

    rgb_image = make_rgb(image_arr[0], image_arr[1], image_arr[2], interval=ManualInterval(vmin=0, vmax=maximum), stretch=SqrtStretch())

    # -------plot RGB image---------
    offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
    extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)
    p = ax.imshow(rgb_image, origin='lower', extent=extent, alpha=1)
    ax.set_aspect('auto') 

    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.text(0.05, 0.05, f'z={args.z:.2f}', c='w', fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)

    # ----------annotate axis---------------
    if args.re_limit is not None:
        limit = args.re_limit * args.re_arcsec
        ax.add_patch(plt.Rectangle((-limit, -limit), 2 * limit, 2 * limit, lw=0.5, color='r', fill=False))
        #ax.add_patch(plt.Circle((0, 0), args.re_arcsec, color='k', fill=False, lw=0.5))

    if not hide_pa:
        pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
        ax = annotate_PAs(pa_arr, ax, fontsize=args.fontsize / args.fontfactor, color='w')

    ax.set_xlim(-args.arcsec_limit, args.arcsec_limit)  # arcsec
    ax.set_ylim(-args.arcsec_limit, args.arcsec_limit)  # arcsec

    if 'segmentation_map_arcsec' in args: ax.contour(args.segmentation_map_arcsec != args.id, levels=0, colors='w', extent=extent, linewidths=0.5) # for over-plotting segmentation map
    if not skip_annotate: ax = annotate_axes(ax, 'arcsec', 'arcsec', args, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax, rgb_image

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''
    image = np.ma.masked_where(image.mask, unp.nominal_values(image.data))

    if takelog:
        new_mask = image <= 0
        image[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        image = np.ma.masked_where(image.mask | new_mask, np.log10(image.data))

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.set_aspect('auto') 
    
    if 'segmentation_map' in args: ax.contour(args.segmentation_map != args.id, levels=0, colors='w' if args.fortalk else 'k', extent=args.extent, linewidths=0.5) # demarcating the segmentation map zone
    if not skip_annotate:
        xylabel = 'arcsec' if args.re_limit is None else r'Offset (R$_e$)'
        ax = annotate_axes(ax, xylabel, xylabel, args, label=label, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_maps(full_hdu, ratio_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, takelog=False):
    '''
    Plots the 2D line ratio maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting the line ratio maps {ratio_labels}..')
    for index, ax in enumerate(axes):
        ratio = ratio_labels[index]
        num_line, den_line = ratio.split('/')
        line_map_num, _, _, _, _ = get_emission_line_map(num_line, full_hdu, args, dered=True, silent=True)
        line_map_den, _, _, _, _ = get_emission_line_map(den_line, full_hdu, args, dered=True, silent=True)

        bad_mask = unp.nominal_values(line_map_den.data) == 0
        line_map_den[bad_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        line_map_den = np.ma.masked_where(line_map_den.mask | bad_mask, line_map_den.data)

        ratio_map = np.ma.masked_where(line_map_num.mask | line_map_den.mask, line_map_num.data / line_map_den.data)
        log_text = r'log ' if takelog else ''
        ax = plot_2D_map(ratio_map, ax, get_ratio_labels(ratio), args, takelog=takelog, clabel=log_text + 'flux ratio', cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_line_maps(full_hdu, line_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, takelog=False):
    '''
    Plots the 2D line flux maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting the line maps {line_labels}..')
    for index, ax in enumerate(axes):
        line = line_labels[index]
        line_map, _, _, _, _ = get_emission_line_map(line, full_hdu, args, dered=True, silent=True)
        log_text = r'$\log$ ' if takelog else ''
        ax = plot_2D_map(line_map, ax, get_line_labels(line), args, takelog=takelog, clabel=log_text + r'SB (ergs/s/cm$^2$/kpc$^2$)', cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_galaxy_example_fig(objid, field, args, fontsize=10, show_log_spectra=False, show_log_map=False):
    '''
    Plots and saves a single figure with the direct image, 1D spectra, emission line maps and emission line ratio maps for a given object
    '''
    args.fontsize = fontsize
    print(f'\nPlotting example galaxy {field}:{objid}..')
    ncol = len(args.line_list) # one each for OII, OIII, Hb and NeIII line

    # -------setting up the figure layout------------
    if len(args.line_list) > 3: fig = plt.figure(figsize=(10, 7) if args.plot_ratio_maps else (10, 5))
    else: fig = plt.figure(figsize=(7.5, 7) if args.plot_ratio_maps else (8, 5))
    fig.subplots_adjust(left=0.09, right=0.9, top=0.94, bottom=0.07)
    outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2 if args.plot_ratio_maps else 1], figure=fig, hspace=0.2) # Outer GridSpec: 2 rows – top (loose), bottom (tight)

    # ---------setting up direct image (1x1 square) and 1D spectra (1x3 wide) subplots----------------
    top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_gs[0], width_ratios=[1, ncol - 1], wspace=0.07)
    ax_dirimg = fig.add_subplot(top_gs[0, 0])
    ax_1dspec = fig.add_subplot(top_gs[0, 1])

    # ---------setting up emission line map axes: ncol x 2 (or ncol x 1) tight subplots with shared axes----------------
    bottom_gs = gridspec.GridSpecFromSubplotSpec(2 if args.plot_ratio_maps else 1, ncol, subplot_spec=outer_gs[1], wspace=0., hspace=0.) # Use sub-GridSpec with tight spacing
    axes_line_maps = [fig.add_subplot(bottom_gs[0, item]) for item in np.arange(ncol)]
    if args.plot_ratio_maps:
        axes_ratio_maps = [fig.add_subplot(bottom_gs[1, item]) for item in np.arange(ncol)]
        #ratios_to_plot = ['OII/NeIII-3867', 'OIII/OII', 'OII/Hb', 'OIII/Hb']
        if len(args.line_list) > 3: ratios_to_plot = ['OIII/OII', 'OII/Hb', 'OIII/Hb', 'Ha/Hb']
        else: ratios_to_plot = ['OIII/OII', 'OII/Hb', 'OIII/Hb']

    # ----------loading the full.fits and 1D.fits files--------------
    full_hdu = load_full_fits(objid, field, args)
    od_hdu = load_1d_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)

    # ----------plotting direct image--------------
    ax_dirimg, _ = plot_rgb_image(full_hdu, ['F200W', 'F150W', 'F115W'], ax_dirimg, args)

    # ----annotating with scalebar----------------
    kpc = 2
    pix = kpc * cosmo.arcsec_per_kpc_proper(args.z).value  # converting kpc to arcsec
    scalebar = AnchoredSizeBar(ax_dirimg.transData, pix, f'{kpc} kpc', 'lower right', pad=0.5, color='w', frameon=False, size_vertical=0.01, fontproperties={'size':args.fontsize / args.fontfactor})
    ax_dirimg.add_artist(scalebar)

    # ----------plotting 1D spectra--------------
    ax_1dspec = plot_1D_spectra(od_hdu, ax_1dspec, args, show_log_flux=show_log_spectra)

    # ----------plotting line flux maps--------------
    flux_map_lim = [10**(-18-0.2), 10**(-16.8+0.2)]
    if show_log_map: flux_map_lim = np.log10(np.array(flux_map_lim))
    #flux_map_lim = [None, None] ##
    axes_line_maps = plot_line_maps(full_hdu, args.line_list, axes_line_maps, args, cmap='pink', vmin=flux_map_lim[0], vmax=flux_map_lim[1], hide_xaxis=args.plot_ratio_maps, takelog=show_log_map)

    # -------zoom-in annotation------------------
    if args.re_limit is not None:
        limit = args.re_limit * args.re_arcsec
        connect1 = ConnectionPatch(xyA=(-limit, -limit), xyB=(axes_line_maps[0].get_xlim()[0], axes_line_maps[0].get_ylim()[1]), coordsA='data', coordsB='data', axesA=ax_dirimg, axesB=axes_line_maps[0], color='r', lw=0.5)
        connect2 = ConnectionPatch(xyA=(limit, -limit), xyB=(axes_line_maps[0].get_xlim()[1], axes_line_maps[0].get_ylim()[1]), coordsA='data', coordsB='data', axesA=ax_dirimg, axesB=axes_line_maps[0], color='r', lw=0.5)
        axes_line_maps[0].add_artist(connect1)
        axes_line_maps[0].add_artist(connect2)
        
    # ----------plotting line ratio image--------------
    if args.plot_ratio_maps:
        ratio_map_lim = [10**-1, 10**1]
        if show_log_map: ratio_map_lim = np.log10(np.array(ratio_map_lim))
        #ratio_map_lim = [None, None] ##
        axes_ratio_maps = plot_line_ratio_maps(full_hdu, ratios_to_plot, axes_ratio_maps, args, cmap='bone', vmin=ratio_map_lim[0], vmax=ratio_map_lim[1], hide_xaxis=False, takelog=show_log_map)

    # ----------annotating and saving figure--------------
    fig.text(0.05, 0.98, f'{field}: ID {objid}', fontsize=args.fontsize, c='k', ha='left', va='top')
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'{field}_{objid}_example_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def get_distance_map_from_AGN_line(x_num, x_den, y_num, y_den, args):
    '''
    Computes and returns the distance of each elements in (x_num, x_den, y_num, y_den) form the AGN demarcation line
    Returned array has same shape as x_num, etc.
    '''
    theoretical_lines, line_labels = get_AGN_func_methods(args)
    dist_method = theoretical_lines[0]

    y_ratio = take_safe_log_ratio(y_num, y_den)
    x_ratio = take_safe_log_ratio(x_num, x_den)

    net_mask = y_ratio.mask | x_ratio.mask

    # ------distance of pixels from AGN line--------------
    sign_map = (unp.nominal_values(y_ratio.data) > AGN_func(unp.nominal_values(x_ratio.data), dist_method)).astype(int)
    sign_map[sign_map == 0] = -1
    distance_from_AGN_line_map = sign_map * get_distance_from_line(unp.nominal_values(x_ratio.data), unp.nominal_values(y_ratio.data), AGN_func, dist_method)
    distance_from_AGN_line_map = np.ma.masked_where(net_mask, distance_from_AGN_line_map)

    return distance_from_AGN_line_map

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_ax(x_num, x_den, y_num, y_den, ax, args, color=None, marker='o', size=20, lw=0.5):
    '''
    Plots line ratio vs line ratio on a given axis
    Returns axis handle and the scatter plot handle
    '''
    if len(np.atleast_1d(x_num).flatten()) > 1:
        distance_map = get_distance_map(np.shape(y_num), args)
    else:
        distance_map = 0

    df = pd.DataFrame({'xnum': unp.nominal_values(np.atleast_1d(x_num)).flatten(), \
                       'xnum_u': unp.std_devs(np.atleast_1d(x_num)).flatten(), \
                       'xden': unp.nominal_values(np.atleast_1d(x_den)).flatten(), \
                       'xden_u': unp.std_devs(np.atleast_1d(x_den)).flatten(), \
                       'ynum': unp.nominal_values(np.atleast_1d(y_num)).flatten(), \
                       'ynum_u': unp.std_devs(np.atleast_1d(y_num)).flatten(), \
                       'yden': unp.nominal_values(np.atleast_1d(y_den)).flatten(), \
                       'yden_u': unp.std_devs(np.atleast_1d(y_den)).flatten(), \
                       'distance': np.atleast_1d(distance_map).flatten()
                       })
    df = df.drop_duplicates().reset_index(drop=True)
    df = df[(df['xnum'] > 0) & (df['xden'] > 0) & (df['ynum'] > 0) & (df['yden'] > 0)]

    y_ratio = unp.log10(unp.uarray(df['ynum'], df['ynum_u']) / unp.uarray(df['yden'], df['yden_u']))
    x_ratio = unp.log10(unp.uarray(df['xnum'], df['xnum_u']) / unp.uarray(df['xden'], df['xden_u']))

    if color is None:
        color = get_distance_map_from_AGN_line(df['xnum'], df['xden'], df['ynum'], df['yden'], args).data
        cmap, vmin, vmax = args.diverging_cmap, -1, 1
    elif color == 'distance':
        color = df['distance']
        cmap, vmin, vmax = 'viridis', None, None

    p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, cmap=cmap, vmin=vmin, vmax=vmax, marker=marker, s=size, lw=lw, edgecolor='w' if args.fortalk else 'k', zorder=10)
    #p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=df['distance'], cmap='viridis', vmin=0, vmax=5, marker=marker, s=size, lw=lw, edgecolor='w' if args.fortalk else 'k', zorder=10)
    ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c='gray', fmt='none', lw=lw, alpha=0.5)

    return ax, p

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_object(full_hdu, args, ax, marker='o', fontsize=10, size=50):
    '''
    Plots the spatially resolved AGN demarcation for a given object on a given axis
    Returns axis handle
    '''    
    # -----------getting the fluxes------------------
    try:
        ynum_map, _, ynum_int, ynum_sum, _ = get_emission_line_map(args.ynum_line, full_hdu, args, silent=True)
        yden_map, _, yden_int, yden_sum, _ = get_emission_line_map(args.yden_line, full_hdu, args, silent=True)

        xnum_map, _, xnum_int, xnum_sum, _ = get_emission_line_map(args.xnum_line, full_hdu, args, silent=True)
        xden_map, _, xden_int, xden_sum, _ = get_emission_line_map(args.xden_line, full_hdu, args, silent=True)
    except:
        print(f'Required emission lines not available for {args.id} with {args.AGN_diag} AGN diagnostic. So skipping this object')
        return ax, None, None
    
    if not args.do_not_correct_flux and args.AGN_diag in ['H21', 'B22'] and args.xden_line == 'Ha': # special treatment for H-alpha line, in order to add the NII 6584 component back
        factor = 0.823  # from grizli source code
        print(f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
        xden_map = np.ma.masked_where(xden_map.mask, xden_map.data / factor)
        xden_int = xden_int / factor
        xden_sum = xden_sum / factor

   # -----------spatially_resolved-----------------------
    ax, scatter_plot_handle = plot_AGN_demarcation_ax(xnum_map, xden_map, ynum_map, yden_map, ax, args, marker=marker, size=size, lw=0.5, color=args.AGN_colorby if 'AGN_colorby' in args else None)

    # -----------integrated-----------------------
    #ax, _ = plot_AGN_demarcation_ax(xnum_int, xden_int, ynum_int, yden_int, ax, args, marker=marker, size=4*size, lw=2)
    ax, _ = plot_AGN_demarcation_ax(xnum_sum, xden_sum, ynum_sum, yden_sum, ax, args, marker=marker, size=4*size, lw=2)

     # -----------2D map inset-----------------------
    ax_inset = ax.inset_axes([0.68, 0.68, 0.3, 0.3])
    distance_from_AGN_line_map = get_distance_map_from_AGN_line(xnum_map, xden_map, ynum_map, yden_map, args)
    ax_inset = plot_2D_map(distance_from_AGN_line_map, ax_inset, '', args, takelog=False, cmap=args.diverging_cmap, vmin=-1, vmax=1, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)
    
    return ax, scatter_plot_handle, ax_inset

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_figure_single(objid, field, args, fontsize=10):
    '''
    Plots and saves the spatially resolved AGN demarcation for a single object
    '''
    args.fontsize = fontsize
    print(f'Plotting AGN demarcation diagram for single object {objid}..')

    # -------setting up the figure--------------------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.1, top=0.98)

    # --------loading the data-------------------------
    full_hdu = load_full_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)
 
    # --------plotting the data-------------------------
    ax, scatter_plot_handle, ax_inset = plot_AGN_demarcation_object(full_hdu, args, ax, marker='o' if 'Par' in field else 's', size=50, fontsize=fontsize)
    ax_inset = annotate_kpc_scale_bar(2, ax_inset, args, label='2 kpc', loc='lower left')

    # -----------annotating axes-----------------------
    theoretical_lines, line_labels = get_AGN_func_methods(args)
    cbar = plt.colorbar(scatter_plot_handle)
    unit_text = r'kpc' if args.re_limit is None else r'R_e'
    cbar.set_label(f'Distance from center ({unit_text})' if 'AGN_colorby' in args and args.AGN_colorby == 'distance' else f'Distance from {theoretical_lines[0]} line', fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax.set_xlim(-1.6, 0.5)
    ax.set_ylim(-1, 2)
    ax.set_xlabel(f'Log {get_ratio_labels("NeIII-3867/OII")}' if args.AGN_diag == 'Ne3O2' else f'Log {get_ratio_labels("SII/NII,Ha")}' if args.AGN_diag == 'H21' else f'Log {get_ratio_labels(f"{args.xnum_line}/{args.xden_line}")}', fontsize=args.fontsize)
    ax.set_ylabel(f'Log {get_ratio_labels("OIII/Hb")}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---------adding literature AGN demarcation lines----------
    color_arr = ['brown', 'darkgreen', 'dodgerblue', 'cyan', 'sienna']
    for index, (theoretical_line, line_label) in enumerate(zip(theoretical_lines, line_labels)):
        overplot_AGN_line_on_BPT(ax, theoretical_line=theoretical_line, label=line_label, color=color_arr[index], fontsize=args.fontsize, lw=0.5 if index else 1, ls='solid')

   # -----------saving figure------------
    ax.text(0.05, 0.95, f'ID #{full_hdu[0].header["ID"]}', fontsize=args.fontsize, c='k', ha='left', va='top', transform=ax.transAxes)
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'BPT_{full_hdu[0].header["ID"]:05d}_upto_{extent_text}_AGNdiag_{args.AGN_diag}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_figure_integrated(df_input, args, fontsize=10):
    '''
    Plots and saves the integrated AGN demarcation for a single object
    '''
    args.fontsize = fontsize
    print(f'Plotting integrated AGN demarcation diagram..')
    theoretical_lines, line_labels = get_AGN_func_methods(args)

    # -------setting up the figure--------------------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.1, top=0.98)

    # -------------correcting Ha back to NII + Ha-----------
    if not args.do_not_correct_flux and args.AGN_diag in ['H21', 'B22'] and args.xden_line == 'Ha': # special treatment for H-alpha line, in order to add the NII 6584 component back
        factor = 0.823  # from grizli source code
        print(f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
        df['Ha'] /= factor
        df['Ha_u'] /= factor

    # --------computing the log ratios-------------------------
    df = df_input.copy()
    df = df[(df[args.xnum_line] > 0) & (df[args.xden_line] > 0) & (df[args.ynum_line] > 0) & (df[args.yden_line] > 0)]
    df['y_ratio'] = unp.log10(unp.uarray(df[args.ynum_line], df[args.ynum_line + '_u']) / unp.uarray(df[args.yden_line], df[args.yden_line + '_u']))
    df['x_ratio'] = unp.log10(unp.uarray(df[args.xnum_line], df[args.xnum_line + '_u']) / unp.uarray(df[args.xden_line], df[args.xden_line + '_u']))
    df['color'] = get_distance_map_from_AGN_line(df[args.xnum_line], df[args.xden_line], df[args.ynum_line], df[args.yden_line], args).data

    # --------plotting the data-------------------------
    for marker in pd.unique(df['marker']):
        df_sub = df[df['marker'] == marker]
        scatter_plot_handle = ax.scatter(unp.nominal_values(df_sub['x_ratio']), unp.nominal_values(df_sub['y_ratio']), c=df_sub['color'], cmap=args.diverging_cmap, vmin=-1, vmax=1, marker=marker, s=100, lw=2, edgecolor='w' if args.fortalk else 'k', label='GLASS (This work)' if marker == 's' else 'PASSAGE (This work)')
        ax.errorbar(unp.nominal_values(df_sub['x_ratio']), unp.nominal_values(df_sub['y_ratio']), xerr=unp.std_devs(df_sub['x_ratio']), yerr=unp.std_devs(df_sub['y_ratio']), c='gray', fmt='none', lw=2, alpha=0.5, zorder=-1)
    
    # --------annotating every point-------------------------
    if args.annotate:
        for index, row in df.iterrows(): ax.text(row['x_ratio'].n, row['y_ratio'].n, f'{row["objid"]}', fontsize=args.fontsize/1.5, c='r', ha='left', va='top')

    # -----------annotating axes-----------------------
    cbar = plt.colorbar(scatter_plot_handle)
    cbar.set_label(f'Distance from maximum upper envelope', fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax.set_xlim(-2.5, 0.5)
    ax.set_ylim(-1, 1)
    ax.set_xlabel(f'Log {get_ratio_labels("NeIII-3867/OII")}' if args.AGN_diag == 'Ne3O2' else f'Log {get_ratio_labels("SII/NII,Ha")}' if args.AGN_diag == 'H21' else f'Log {get_ratio_labels(f"{args.xnum_line}/{args.xden_line}")}', fontsize=args.fontsize)
    ax.set_ylabel(f'Log {get_ratio_labels("OIII/Hb")}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---------adding literature AGN demarcation lines----------
    color_arr = ['brown', 'darkgreen', 'dodgerblue', 'cyan', 'sienna']
    for index, (theoretical_line, line_label) in enumerate(zip(theoretical_lines, line_labels)):
        overplot_AGN_line_on_BPT(ax, theoretical_line=theoretical_line, label=line_label, color=color_arr[index], fontsize=args.fontsize, lw=0.5 if index else 1, ls='solid')

   # -----------saving figure------------
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'BPT_integrated_upto_{extent_text}_AGNdiag_{args.AGN_diag}.png'
    save_fig(fig, figname, args)

    return df

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_figure_multiple(objlist, args, fontsize=10, exclude_ids=[]):
    '''
    Plots and saves the spatially resolved AGN demarcation for a single object
    '''
    args.fontsize = fontsize
    objlist = [item for item in objlist if not item[1] in exclude_ids]
    print(f'Plotting multi-panel AGN demarcation diagram for {len(objlist)} objects..')
    theoretical_lines, line_labels = get_AGN_func_methods(args)        
    color_arr = ['brown', 'darkgreen', 'dodgerblue', 'cyan', 'sienna']

    # -------setting up the figure--------------------
    nrow, ncol = 2, 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 5), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.92, bottom=0.1, top=0.98, wspace=0., hspace=0.)

    # --------looping over all objects-------------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'Doing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')
        thisrow = int(index / ncol)
        thiscol = index % ncol
        ax = axes[thisrow][thiscol]
    
        # --------loading the data-------------------------
        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args, field=field)
    
        # --------plotting the data-------------------------
        ax, scatter_plot_handle, ax_inset = plot_AGN_demarcation_object(full_hdu, args, ax, marker='o' if 'Par' in field else 's', size=20, fontsize=fontsize)
        if ax_inset is not None: ax_inset = annotate_kpc_scale_bar(2, ax_inset, args, label=None, color='brown', loc='upper left')

        # -----------annotating axes-----------------------
        ax.set_xlim(-1.6, 0.5)
        ax.set_ylim(-1, 2)

        if thisrow == nrow - 1: ax.set_xlabel(f'Log {get_ratio_labels("NeIII-3867/OII")}' if args.AGN_diag == 'Ne3O2' else f'Log {get_ratio_labels("SII/NII,Ha")}' if args.AGN_diag == 'H21' else f'Log {get_ratio_labels(f"{args.xnum_line}/{args.xden_line}")}', fontsize=args.fontsize)
        if thiscol == 0: ax.set_ylabel(f'Log {get_ratio_labels("OIII/Hb")}', fontsize=args.fontsize)
        
        ax.tick_params(axis='both', which='major', labelsize=args.fontsize)
        ax.text(0.05, 0.95, f'ID #{objid}', fontsize=args.fontsize, c='k', ha='left', va='top', transform=ax.transAxes)

       # ---------adding literature AGN demarcation lines----------
        for index2, (theoretical_line, line_label) in enumerate(zip(theoretical_lines, line_labels)):
            overplot_AGN_line_on_BPT(ax, theoretical_line=theoretical_line, label=line_label if index == len(objlist) - 1 else None, color=color_arr[index2], fontsize=args.fontsize, lw=0.5 if index2 else 1, ls='solid')

    # -------adding colorbar for entire figure-----------
    cax = fig.add_axes([0.92 + 0.01, 0.1, 0.01, 0.98 - 0.1])
    if scatter_plot_handle is not None:
        cbar = plt.colorbar(scatter_plot_handle, cax=cax, orientation='vertical')
        unit_text = r'kpc' if args.re_limit is None else r'R_e'
        cbar.set_label(f'Distance from center ({unit_text})' if 'AGN_colorby' in args and args.AGN_colorby == 'distance' else f'Distance from {theoretical_lines[0]} line', fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)
    
    # -----------saving figure------------
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'BPT_all_objects_upto_{extent_text}_AGNdiag_{args.AGN_diag}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_map(field, objid, Zdiag, args):
    '''
    Loads the pre-computed 2D metallicity map for a given object
    Returns 2D image
    '''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{",".join(args.exclude_lines)}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    AGN_diag_text = f'_AGNdiag_{args.AGN_diag}'if args.AGN_diag != 'None' else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in Zdiag else ''
    dered_text = '_dered' if args.dered_in_NB and 'NB' in Zdiag else ''
    output_fitsname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'{field}_{objid:05d}_logOH_map_upto_{extent_text}{args.snr_text}{args.only_seg_text}{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}{AGN_diag_text}{exclude_text}{C25_text}{dered_text}.fits'
    ####################################
    if 'glass' in field and 'SNR_4.0' in str(output_fitsname):
        output_fitsname = Path(str(output_fitsname).replace('SNR_4.0', 'SNR_2.0'))
        print(f'\nWARNING: Actually choosing metallicity map corresponding to vorbin SNR=2 for {field}-{objid}') ##
    ####################################

    # ------reading in existing outputfile--------------
    print(f'\nReading metallicity map from existing fits file {output_fitsname}')
    hdul = fits.open(output_fitsname)
    hdu = hdul['log_OH']
    logOH_map = np.ma.masked_where(np.isnan(hdu.data), unp.uarray(hdu.data, hdul['log_OH_u'].data))
    if hdu.header['log_oh_int'] is None: logOH_int = ufloat(np.nan, np.nan)
    else: logOH_int = ufloat(hdu.header['log_oh_int'], hdu.header['log_oh_int_err'])
    if hdu.header['log_oh_sum'] is None: logOH_sum = ufloat(np.nan, np.nan)
    else: logOH_sum = ufloat(hdu.header['log_oh_sum'], hdu.header['log_oh_sum_err'])

    return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_df(field, objid, Zdiag, args, ensure_sii=False):
    '''
    Loads the pre-computed list of spatially resolved metallicity values for a given object as a dataframe
    Returns dataframe
    '''

    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    Zbranch_text = '' if np.array([item in Zdiag for item in ['NB', 'P25']]).any() else f'-{args.Zbranch}'
    exclude_text = f'_without_{",".join(args.exclude_lines)}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''
    AGN_diag_text = f'_AGNdiag_{args.AGN_diag}'if args.AGN_diag != 'None' else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in Zdiag else ''
    dered_text = '_dered' if args.dered_in_NB and 'NB' in Zdiag else ''
    output_fitsname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'{field}_{objid:05d}_logOH_map_upto_{extent_text}{args.snr_text}{args.only_seg_text}{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}{AGN_diag_text}{exclude_text}{C25_text}{dered_text}.fits'
    ####################################
    if 'glass' in field and 'SNR_4.0' in str(output_fitsname):
        output_fitsname = Path(str(output_fitsname).replace('SNR_4.0', 'SNR_2.0'))
        print(f'\nWARNING: Actually choosing metallicity map corresponding to vorbin SNR=2 for {field}-{objid}') ##
    ####################################

    # ------reading in existing outputfile--------------
    print(f'\nReading metallicity map from existing fits file {output_fitsname}')
    hdul = fits.open(output_fitsname)
    hdu = hdul['log_OH']
    if ensure_sii and 'SII' not in hdu.header['LINES']:
        print(f'Required the presence of SII lines but no SII line in dits file for {field}:{objid}, therefore returning empty df')
        return None, ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)

    logOH_df = Table(hdul['tab'].data).to_pandas()

    if hdu.header['log_oh_int'] is None: logOH_int = ufloat(np.nan, np.nan)
    else: logOH_int = ufloat(hdu.header['log_oh_int'], hdu.header['log_oh_int_err'])
    if hdu.header['log_oh_sum'] is None: logOH_sum = ufloat(np.nan, np.nan)
    else: logOH_sum = ufloat(hdu.header['log_oh_sum'], hdu.header['log_oh_sum_err'])

    return logOH_df, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def append_to_logOH_fit_df(field, objid, logOH_int, logOH_sum, linefit_original, linefit_odr, linefit_lenstronomy, linefit_mcmc, params_llim, params_ulim, params_median, Zdiag, args):
    '''
    Looks for existing csv file with radial fitting info for metallicity and appends new fitted parameters to the csv file
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    output_dfname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'logOH_allfits.csv'
    
    df_logOH_fit = pd.DataFrame({'field': field, 'objid': objid, \
                                    'logOH_sum': logOH_sum.n, 'logOH_sum_u': logOH_sum.s, 'logOH_int': logOH_int.n, 'logOH_int_u': logOH_int.s, \
                                    'logOH_cen': linefit_original[1].n, 'logOH_cen_u': linefit_original[1].s, 'logOH_slope': linefit_original[0].n, 'logOH_slope_u': linefit_original[0].s, \
                                    'logOH_cen_wls': linefit_odr[1].n, 'logOH_cen_wls_u': linefit_odr[1].s, 'logOH_slope_wls': linefit_odr[0].n, 'logOH_slope_wls_u': linefit_odr[0].s, \
                                    'logOH_cen_lenstronomy': linefit_lenstronomy[1].n, 'logOH_cen_lenstronomy_u': linefit_lenstronomy[1].s, 'logOH_slope_lenstronomy': linefit_lenstronomy[0].n, 'logOH_slope_lenstronomy_u': linefit_lenstronomy[0].s, \
                                    'logOH_cen_mcmc': linefit_mcmc[1].n, 'logOH_cen_mcmc_u': linefit_mcmc[1].s, 'logOH_slope_mcmc': linefit_mcmc[0].n, 'logOH_slope_mcmc_u': linefit_mcmc[0].s, \
                                    'logOH_cen_mcmc_llim': params_llim[1], 'logOH_cen_mcmc_ulim': params_ulim[1], 'logOH_slope_mcmc_llim': params_llim[0], 'logOH_slope_mcmc_ulim': params_ulim[0], \
                                    'q_mcmc': params_median[2], 'q_mcmc_llim': params_llim[2], 'q_mcmc_ulim': params_ulim[2], \
                                    'pa_mcmc': params_median[3], 'pa_mcmc_llim': params_llim[3], 'pa_mcmc_ulim': params_ulim[3], \
                                    'Zdiag': f'{Zdiag}_C25' if args.use_C25 and 'NB' not in Zdiag else f'{Zdiag}', 'Zdiag_branch': args.Zbranch, 'AGN_diag': args.AGN_diag, 'dered_in_NB': args.dered_in_NB, \
                                    'extent_unit': 'arcsec' if args.re_limit is None else 're', 'extent_value': args.arcsec_limit if args.re_limit is None else args.re_limit, \
                                    'snr_text': args.snr_text, 'only_seg': args.only_seg, 'vorbin_text': args.vorbin_text}, \
                                    index=[0])
    df_logOH_fit.to_csv(output_dfname, index=None, mode='a', header=not os.path.exists(output_dfname))
    print(f'Appended metallicity fit to catalog file {output_dfname}')

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fig_single(objid, field, Zdiag, args, fontsize=10, do_mcmc=False):
    '''
    Plots and saves a single figure with the direct image, 2D metallicity map and metallicity radial profile for a given object
    '''
    args.fontsize = fontsize
    print(f'Plotting metallicity ({Zdiag}) figure for {objid}..')

    # -----------loading the data---------------
    full_hdu = load_full_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)
    logOH_map, _, _ = load_metallicity_map(field, objid, Zdiag, args)
    logOH_df, _, _ =  load_metallicity_df(field, objid, Zdiag, args)

    # --------setting up the figure------------
    fig = plt.figure(figsize=(9, 3))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.96, bottom=0.14, wspace=0.3)
    outer_gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1], figure=fig)

    left_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_gs[0], wspace=0.1)
    right_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_gs[1])
    axes = [fig.add_subplot(left_gs[0, 0]), fig.add_subplot(left_gs[0, 1]), fig.add_subplot(right_gs[0, 0])]

    # -----plotting direct image-----------------
    axes[0], _ = plot_rgb_image(full_hdu, ['F200W', 'F150W', 'F115W'], axes[0], args)

    # -----plotting 2D metallicity map-----------
    Zlim = [7.0, 9.0]
    axes[1] = plot_2D_map(logOH_map, axes[1], f'log O/H + 12 ({Zdiag})', args, clabel='', takelog=False, cmap='cividis', vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=False, hide_yaxis=True, hide_cbar=False)
    axes[1] = annotate_kpc_scale_bar(2, axes[1], args, label='2 kpc', loc='lower right')

    # ------plotting metallicity radial profile-----------
    axes[2], _ = plot_radial_profile(logOH_df, axes[2], args, ylim=Zlim, xlim=[0, 5 if args.re_limit is None else args.re_limit], Zdiag=Zdiag, do_mcmc=do_mcmc)

    # -------zoom-in annotation------------------
    if args.re_limit is not None:
        limit = args.re_limit * args.re_arcsec
        connect1 = ConnectionPatch(xyA=(limit, limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[1]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
        connect2 = ConnectionPatch(xyA=(limit, -limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[0]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
        axes[1].add_artist(connect1)
        axes[1].add_artist(connect2)
        
    # -----------saving figure------------
    axes[2].text(0.05, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[2].transAxes)
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in Zdiag else ''
    figname = f'metallicity_{objid:05d}_upto_{extent_text}{C25_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fig_multiple(objlist, Zdiag, args, fontsize=10, do_mcmc=False):
    '''
    Plots and saves a single multi-panel figure with the direct image, 2D metallicity map and metallicity radial profile for a given list of objects
    Also saves the fitted metallicity gradient and related info in a csv file
    '''
    args.fontsize = fontsize
    print(f'Plotting metallicity ({Zdiag}) figure for {len(objlist)} objects..')

   # -------setting up the figure--------------------
    Zlim, cmap = [7.0, 9.1], 'cividis'
    ncol = 2
    nrow = int(np.ceil(len(objlist) / ncol))
    fig = plt.figure(figsize=(10, 6))
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.07, top=0.93, wspace=0., hspace=0.)
    outer_gs = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.33, hspace=0.)

    # --------looping over all objects-------------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'\nDoing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')
        thisrow = int(index / ncol)
        thiscol = index % ncol
        show_xaxis = (thisrow == nrow - 1) or (index == len(objlist) - 1) or ((len(objlist) - index - 1 ) < ncol)
        hide_xaxis = not show_xaxis

        # --------setting up the sub-figure------------
        object_gs = outer_gs[index].subgridspec(1, 3, wspace=0., hspace=0.)
        axes = [fig.add_subplot(object_gs[0, item]) for item in range(3)]

        # -----------loading the data---------------
        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args, field=field)
        logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)
        logOH_df, logOH_int, logOH_sum =  load_metallicity_df(field, objid, Zdiag, args)

        # -----plotting direct image-----------------
        axes[0], _ = plot_rgb_image(full_hdu, ['F200W', 'F150W', 'F115W'], axes[0], args, hide_xaxis=hide_xaxis, hide_yaxis=False, hide_cbar=True, hide_filter_names=index, hide_pa=True)
        axes[0].text(0.95, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='w', ha='right', va='top', transform=axes[0].transAxes)
        
        # -----plotting 2D metallicity map-----------
        axes[1] = plot_2D_map(logOH_map, axes[1], None, args, clabel='', takelog=False, cmap=cmap, vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=hide_xaxis, hide_yaxis=True, hide_cbar=True)
        axes[1] = annotate_kpc_scale_bar(2, axes[1], args, label='2 kpc', color='brown', loc='lower right')
    
        # -------zoom-in annotation------------------
        if args.re_limit is not None:
            limit = args.re_limit * args.re_arcsec
            connect1 = ConnectionPatch(xyA=(limit, limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[1]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
            connect2 = ConnectionPatch(xyA=(limit, -limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[0]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
            axes[1].add_artist(connect1)
            axes[1].add_artist(connect2)
        
        # ------plotting metallicity radial profile-----------
        axes[2], linefit = plot_radial_profile(logOH_df, axes[2], args, ylim=Zlim, xlim=[0, 5 if args.re_limit is None else args.re_limit], hide_xaxis=hide_xaxis, hide_yaxis=False, hide_cbar=True, short_label=False, Zdiag=Zdiag, do_mcmc=do_mcmc)
        linefit_original, linefit_odr, linefit_lenstronomy, params_llim, params_median, params_ulim = linefit
        linefit_mcmc = unp.uarray(params_median[:2], np.mean([np.array(params_median[:2]) - np.array(params_llim[:2]), np.array(params_ulim[:2]) - np.array(params_median[:2])], axis=0))
        axes[2].yaxis.set_label_position('right')
        axes[2].yaxis.tick_right()
   
        # ----------append fit to dataframe and save-----------
        append_to_logOH_fit_df(field, objid, logOH_int, logOH_sum, linefit_original, linefit_odr, linefit_lenstronomy, linefit_mcmc, params_llim, params_ulim, params_median, Zdiag, args)

    # -------making colorbars for the entire fig----------
    cax = fig.add_axes([0.07, 0.96, 0.93 - 0.07, 0.01])    
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=mplcolors.Normalize(vmin=Zlim[0], vmax=Zlim[1]), orientation='horizontal')
    cbar.set_label(r'$\log$' + f' (O/H) + 12 [{Zdiag}]', fontsize=args.fontsize)
    cax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=args.fontsize)
    
    # -----------saving figure------------
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in Zdiag else ''
    dered_text = '_dered' if args.dered_in_NB and 'NB' in Zdiag else ''
    figname = f'metallicity_multi_panel{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}_upto_{extent_text}{C25_text}{dered_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def correct_ha_C20(N2_plus_Ha_map, logOH_map):
    '''
    Extract the H-alpha flux from the N2+Halpha compund, following C20 metallicity calibration
    '''
    #log_N2Ha_logOH_poly_coeff = [1, -10] # from approx fit to MAPPINGS models
    log_N2Ha_logOH_poly_coeff = [-0.489, 1.513, -2.554, -5.293, -2.867][::-1] # from Table 2 N2 row of Curti+2019

    if hasattr(N2_plus_Ha_map, "__len__"): # if it is an array
        logOH_arr = logOH_map.data.flatten()
        log_N2Ha_arr = []
        for logOH in logOH_arr:
            try:
                #log_N2Ha = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH)
                log_N2Ha = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH - 8.69) # because using C19 N2 calibration
                log_N2Ha_arr.append(log_N2Ha)
            except:
                log_N2Ha_arr.append(np.nan)
        log_N2Ha_map = np.ma.masked_where(logOH_map.mask, np.reshape(log_N2Ha_arr, np.shape(logOH_map)))
        N2Ha_map = np.ma.masked_where(log_N2Ha_map.mask, 10 ** log_N2Ha_map.data)

        Ha_data = N2_plus_Ha_map.data / (1 + N2Ha_map.data)
        Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask | N2Ha_map.mask, Ha_data)
    else: # for single values
        log_N2Ha_int = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH_map - 8.69)
        N2Ha_int = 10 ** log_N2Ha_int
        Ha_map = N2_plus_Ha_map / (1 + N2Ha_int)

    return Ha_map

# --------------------------------------------------------------------------------------------------------------------
def correct_ha_F18(N2_plus_Ha_map, args):
    '''
    Extract the H-alpha flux from the N2+Halpha compund, following F18 N2/Ha calibration
    '''
    def F18_func(x, redshift, log_mass): # Eq 3 of Faisst+2018
        psi = x + 0.138 - 0.042 * (1 + redshift) **2 
        solve = 3.696 * psi + 3.236 * psi**(-1) + 0.729 * psi ** (-2) + 14.928 + 0.156 * (1 + redshift)**2 - log_mass
        return solve

    log_N2Ha = brentq(F18_func, -2, 0, args=(args.z, args.log_mass))
    N2Ha = 10 ** log_N2Ha

    if hasattr(N2_plus_Ha_map, "__len__"): # if it is an array
        Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask, N2_plus_Ha_map.data / (1 + 1.333 * N2Ha)) # 1.333 factor is because F18 consider BOTH NII components, assuming NII 6548 = N II 6584 / 3
        Ha_map.data[Ha_map.data==0] = np.nan
    else: # for single values
        Ha_map = N2_plus_Ha_map / (1 + 1.333 * N2Ha)

    return Ha_map

# --------------------------------------------------------------------------------------------------------------------
def get_corrected_ha(N2_plus_Ha_map, logOH_map, args):
    '''
    Computes the corrected Halpha map given the metallicity map and N2+Halpha map
    Returns corrected and uncorrcted Ha and SFR maps
    '''
    Ha_map = correct_ha_C20(N2_plus_Ha_map, logOH_map)
    #Ha_map = correct_ha_F18(N2_plus_Ha_map, args)

    return Ha_map

# --------------------------------------------------------------------------------------------------------------------
def get_corrected_sfr(full_hdu, logOH_map, args, logOH_int=None, logOH_sum=None):
    '''
    Computes the corrected SFR map given the metallicity map and full_hdu (from whic it obtains the uncorrected Halpha map)
    Returns corrected and uncorrcted Ha and SFR maps
    '''
    distance = cosmo.luminosity_distance(args.z)

    # --------deriving the Ha-based SFR where Halpha available---------
    if 'Ha' in args.available_lines:
        # -------deriving H-alpha map-----------
        N2_plus_Ha_map, _, N2_plus_Ha_int, N2_plus_Ha_sum, _ = get_emission_line_map('Ha', full_hdu, args, silent=True)
        sfr_map = compute_SFR(N2_plus_Ha_map, distance) # N2_plus_Ha_map here is really Ha_map, because the correction has not been undone yet
        sfr_map.data[sfr_map.data==0] = np.nan

        # ----------correcting Ha map------------
        if not args.do_not_correct_flux:
            factor = 0.823 # from James et al. 2023?
            N2_plus_Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask, N2_plus_Ha_map.data / factor)

        Ha_map = get_corrected_ha(N2_plus_Ha_map, logOH_map, args)

        # -------deriving SFR map-----------
        sfr_map_corrected = compute_SFR(Ha_map, distance)

        # -------deriving the integrated SFR------------
        if logOH_int is not None:
            Ha_int = get_corrected_ha(N2_plus_Ha_int, logOH_int, args)
            sfr_int_corrected = compute_SFR(Ha_int, distance)
        else:
            sfr_int_corrected = ufloat(np.nan, np.nan)

        if logOH_sum is not None:
            Ha_sum = get_corrected_ha(N2_plus_Ha_sum, logOH_sum, args)
            sfr_sum_corrected = compute_SFR(Ha_sum, distance)
        else:
            sfr_sum_corrected = ufloat(np.nan, np.nan)

        return N2_plus_Ha_map, Ha_map, sfr_map, sfr_map_corrected, sfr_int_corrected, sfr_sum_corrected

    # --------deriving the OII-based SFR in case Halpha unavailable---------
    elif 'OII' in args.available_lines:
        print(f'\nH-alpha unavailable for ID {args.id}, so deriving SFR from OII..')
        OII_map, _, OII_int, OII_sum, _ = get_emission_line_map('OII', full_hdu, args, silent=True)
        sfr_map = get_sfr_from_oii(OII_map.data, args.z)
        sfr_map = np.ma.masked_where(OII_map.mask, sfr_map)

        sfr_int = get_sfr_from_oii(OII_int, args.z)
        sfr_sum = get_sfr_from_oii(OII_sum, args.z)

        return None, OII_map, None, sfr_map, sfr_int, sfr_sum

# --------------------------------------------------------------------------------------------------------------------
def append_to_sfr_fit_df(field, objid, sfr_int, sfr_sum, linefit, Zdiag, args):
    '''
    Looks for existing csv file with radial fitting info for SFR-metallicity correlations and appends new fitted parameters to the csv file
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    output_dfname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'logOH_sfr_fits.csv'

    df_Zsfr_fit = pd.DataFrame({'field': field, 'objid': objid, 'SFR_int': sfr_int.n, 'SFR_int_u': sfr_int.s, 'SFR_sum': sfr_sum.n, 'SFR_sum_u': sfr_sum.s, 'logZ_logSFR_cen': linefit[1].n, \
                                'logZ_logSFR_cen_u': linefit[1].s, 'logZ_logSFR_slope': linefit[0].n, 'logZ_logSFR_slope_u': linefit[0].s, 'Zdiag': Zdiag, 'Zdiag_branch': args.Zbranch, 'AGN_diag': args.AGN_diag, \
                                'extent_unit': 'arcsec' if args.re_limit is None else 're', 'extent_value': args.arcsec_limit if args.re_limit is None else args.re_limit, 'dered_in_NB': args.dered_in_NB, \
                                'snr_text': args.snr_text, 'only_seg': args.only_seg, 'vorbin_text': args.vorbin_text}, \
                                index=[0])
    df_Zsfr_fit.to_csv(output_dfname, index=None, mode='a', header=not os.path.exists(output_dfname))
    print(f'Appended metallicity-sfr fit to catalog file {output_dfname}')

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_fig_single(objid, field, Zdiag, args, fontsize=10):
    '''
    Plots and saves a single figure with the SFR map, metallicity vs SFR plot for a given object
    '''
    args.fontsize = fontsize
    print(f'Plotting metallicity-SFR figure for {objid}..')
        
    # -----------loading the data---------------
    full_hdu = load_full_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)
    logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)
    Zlim = [7.0, 9.0]
    log_sfr_lim = [-2, 0.7]

    # --------setting up the figure------------
    fig, axes = plt.subplots(3 if args.debug_Zsfr else 1, 2, figsize=(5.7, 7) if args.debug_Zsfr else (7, 3))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.15, wspace=0.5, hspace=0.3)

    # ----------getting the SFR maps------------------
    N2_plus_Ha_map, Ha_map, sfr_map, sfr_map_corrected, sfr_int, sfr_sum = get_corrected_sfr(full_hdu, logOH_map, args, logOH_int=logOH_int, logOH_sum=logOH_sum)
    log_sfr_map = np.ma.masked_where(sfr_map_corrected.mask | logOH_map.mask, unp.log10(sfr_map_corrected.data))
   
    # ------plotting native Ha and SFR maps------------
    if args.debug_Zsfr:
        axes[0][0] = plot_2D_map(N2_plus_Ha_map, axes[0][0], r'H$\alpha$ (orig)', args, takelog=True, cmap='RdPu_r', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=False, hide_cbar=False)
        axes[0][1] = plot_2D_map(sfr_map, axes[0][1], f'log SFR (orig)', args, takelog=True, cmap='winter', vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=True, hide_yaxis=True, hide_cbar=False)

    # -----plotting 2D metallicity map and corrected Halpha map-----------
    if args.debug_Zsfr:
        axes[1][0] = plot_2D_map(logOH_map, axes[1][0], f'log O/H + 12 ({Zdiag})', args, takelog=False, cmap='cividis', vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=True, hide_yaxis=False, hide_cbar=False)
        axes[1][1] = plot_2D_map(Ha_map, axes[1][1], r'H$\alpha$ (corrected)', args, takelog=True, cmap='RdPu_r', vmin=-20, vmax=-18, hide_xaxis=False, hide_yaxis=True, hide_cbar=False)
        axes = axes[2]

    # -----plotting 2D SFR map-----------
    axes[0] = plot_2D_map(log_sfr_map, axes[0], r'$\log(\Sigma_{\rm SFR})$ (corrected)' if args.debug_Zsfr else r'$\log(\Sigma_{\rm SFR})$', args, clabel=r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', takelog=False, cmap='winter', vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=False, hide_yaxis=False, hide_cbar=False)

    # ------plotting metallicity vs SFR-----------
    df = pd.DataFrame({'logOH': unp.nominal_values(np.ma.compressed(logOH_map)), 'logOH_u': unp.std_devs(np.ma.compressed(logOH_map)), 'log_sfr':unp.nominal_values(np.ma.compressed(log_sfr_map)), 'log_sfr_u': unp.std_devs(np.ma.compressed(log_sfr_map))})
    if args.vorbin:
        df['bin_ID'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, args.voronoi_bin_IDs.data))
        df = df.groupby('bin_ID', as_index=False).agg(np.mean)

    axes[1].scatter(df['log_sfr'], df['logOH'], c='grey', s=20, alpha=1)
    axes[1].errorbar(df['log_sfr'], df['logOH'], xerr=df['log_sfr_u'], yerr=df['logOH_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

    # -------Z vs SFR fitting-------------
    fit_color = 'cornflowerblue'
    linefit = odr_fit(df, quant_x='log_sfr', quant_y='logOH')
    axes[1] = plot_fitted_line(axes[1], linefit, df['log_sfr'], fit_color, args, short_label=False, index=0, label=f'Slope = {linefit[0]: .2f}')

    # --------annotating axis--------------
    axes[1].set_xlim(log_sfr_lim[0], log_sfr_lim[1]) # kpc
    axes[1].set_ylim(Zlim[0], Zlim[1])
    axes[1].set_box_aspect(1)
    axes[1] = annotate_axes(axes[1], r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', r'$\log$ (O/H) + 12', args)

    # ----------append fit to dataframe and save-----------
    append_to_sfr_fit_df(field, objid, sfr_int, sfr_sum, linefit, Zdiag, args)

    # -----------saving figure------------
    axes[1].text(0.05, 0.85, r'$\log$(M/M$_{\odot}$) =' + f'{args.log_mass: .2f}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[1].transAxes)
    axes[1].text(0.05, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[1].transAxes)
    debug_text = '_debug' if args.debug_Zsfr else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'metallicity-sfr_{objid:05d}{debug_text}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_fig_multiple(objlist, Zdiag, args, fontsize=10, exclude_ids=[]):
    '''
    Plots and saves a single figure with the SFR map, metallicity vs SFR plot for a given list of objects
    '''
    args.fontsize = fontsize
    objlist = [item for item in objlist if not item[1] in exclude_ids]
    print(f'Plotting metallicity-SFR figure for {len(objlist)} objects..')

    Zlim = [7.0, 9.0]
    log_sfr_lim = [-2, 0.7]

   # -------setting up the figure--------------------
    nrow, ncol = 4, 2
    fig = plt.figure(figsize=(9, 6) if nrow == 3 else (9, 7.5))
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.08, top=0.9, wspace=0., hspace=0.)
    outer_gs = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.4, hspace=0.)

    # --------looping over all objects-------------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'Doing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')

        # --------setting up the sub-figure------------
        object_gs = outer_gs[index].subgridspec(1, 2, wspace=0., hspace=0.)
        axes = [fig.add_subplot(object_gs[0, item]) for item in range(2)]

        # -----------loading the data---------------
        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args, field=field)
        logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)

        # ----------getting the SFR maps------------------
        _, _, _, sfr_map, sfr_int, sfr_sum = get_corrected_sfr(full_hdu, logOH_map, args, logOH_int=logOH_int, logOH_sum=logOH_sum)
        sfr_map.data[sfr_map.data < 0] = np.nan
        log_sfr_map = np.ma.masked_where(sfr_map.mask | logOH_map.mask, unp.log10(sfr_map.data))
    
        # -----plotting 2D SFR map-----------
        cmap = 'winter'
        axes[0] = plot_2D_map(log_sfr_map, axes[0], r'$\log$(M/M$_{\odot}$) =' + f'{args.log_mass: .2f}', args, clabel=r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', takelog=False, cmap=cmap, vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=index < nrow - 1, hide_yaxis=False, hide_cbar=True)
        axes[0] = annotate_kpc_scale_bar(2, axes[0], args, label='2 kpc', color='brown', loc='lower right')

        # ------plotting metallicity vs SFR-----------
        df = pd.DataFrame({'logOH': unp.nominal_values(np.ma.compressed(logOH_map)), 'logOH_u': unp.std_devs(np.ma.compressed(logOH_map)), 'log_sfr':unp.nominal_values(np.ma.compressed(log_sfr_map)), 'log_sfr_u': unp.std_devs(np.ma.compressed(log_sfr_map))})
        if args.vorbin:
            df['bin_ID'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, args.voronoi_bin_IDs.data))
            df = df.groupby('bin_ID', as_index=False).agg(np.mean)
        df = df.dropna(axis=0)

        axes[1].scatter(df['log_sfr'], df['logOH'], c='grey', s=20, alpha=1)
        axes[1].errorbar(df['log_sfr'], df['logOH'], xerr=df['log_sfr_u'], yerr=df['logOH_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

        # -------Z vs SFR fitting-------------
        fit_color = 'cornflowerblue'
        linefit = odr_fit(df, quant_x='log_sfr', quant_y='logOH')
        axes[1] = plot_fitted_line(axes[1], linefit, df['log_sfr'], fit_color, args, short_label=False, index=0, label=f'Slope = {linefit[0]: .2f}')

        # --------annotating axis--------------
        axes[1].text(0.05, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[1].transAxes)
        axes[1].set_xlim(log_sfr_lim[0], log_sfr_lim[1]) # kpc
        axes[1].set_ylim(Zlim[0], Zlim[1])
        axes[1] = annotate_axes(axes[1], r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', r'$\log$ (O/H) + 12', args, hide_xaxis=index < nrow - 1, hide_yaxis=False, hide_cbar=True)
        axes[1].yaxis.set_label_position('right')
        axes[1].yaxis.tick_right()

        # ----------append fit to dataframe and save-----------
        append_to_sfr_fit_df(field, objid, sfr_int, sfr_sum, linefit, Zdiag, args)

   # -------making colorbars for the entire fig----------
    cax = fig.add_axes([0.07, 0.95, 0.93 - 0.07, 0.01])    
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=mplcolors.Normalize(vmin=log_sfr_lim[0], vmax=log_sfr_lim[1]), orientation='horizontal')
    cbar.set_label(r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', fontsize=args.fontsize)
    cax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=args.fontsize)
    
     # -----------saving figure------------
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'metallicity-sfr_multi_panel_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_radial_profile_fig_single(objid, field, Zdiag, args, fontsize=10):
    '''
    Plots and saves a single figure with the 2D metallicity map, metallicity radial profile, SFR map, SFR radial profile and metallicity vs SFR plot for a given object
    '''
    args.fontsize = fontsize
    args.fontfactor = 1.3
    print(f'Plotting metallicity-SFR radial profile figure for {field}-{objid}..')

    # -------setting up the figure--------------------
    fig = plt.figure(figsize=(10, 5))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.1, top=0.93, wspace=0., hspace=0.)
    outer_gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], figure=fig, wspace=0.41, hspace=0.)

    # --------setting up the sub-figure------------
    map_radprof_gs = outer_gs[0].subgridspec(2, 2, wspace=0., hspace=0.)
    axes = [fig.add_subplot(map_radprof_gs[0, 0]), fig.add_subplot(map_radprof_gs[0, 1]), fig.add_subplot(map_radprof_gs[1, 0]), fig.add_subplot(map_radprof_gs[1, 1])]
    Zsfr_gs = outer_gs[1].subgridspec(1, 1, wspace=0., hspace=0.)
    axes += [fig.add_subplot(Zsfr_gs[0, 0])]
    
    # -----------loading the data---------------
    full_hdu = load_full_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)
    logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)
    logOH_df, _, _ = load_metallicity_df(field, objid, Zdiag, args)
    Zlim = [7., 8.9]
    log_sfr_lim = [-1., 0.5]

    # ----------getting the SFR maps------------------
    _, _, _, sfr_map, sfr_int, sfr_sum = get_corrected_sfr(full_hdu, logOH_map, args, logOH_int=logOH_int, logOH_sum=logOH_sum)
    log_sfr_map = np.ma.masked_where(sfr_map.mask | logOH_map.mask, unp.log10(sfr_map.data))

    # -----plotting 2D metallicity map-----------
    axes[0] = plot_2D_map(logOH_map, axes[0], 'log(O/H) + 12 [NB]', args, clabel='', takelog=False, cmap='cividis', vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=True, hide_yaxis=False, hide_cbar=False, hide_cbar_ticks=True, cticks_integer=False)
    axes[0] = annotate_kpc_scale_bar(2, axes[0], args, label='2 kpc', color='brown', loc='lower right')

    # # ------plotting metallicity radial profile-----------
    axes[1], _ = plot_radial_profile(logOH_df, axes[1], args, ylim=Zlim, xlim=[0, 5 if args.re_limit is None else args.re_limit], hide_xaxis=True, hide_yaxis=False, hide_cbar=True, short_label=False, do_mcmc=True)
    axes[1].yaxis.set_label_position('right')
    axes[1].yaxis.tick_right()
    
    # -----plotting 2D SFR map-----------
    axes[2] = plot_2D_map(log_sfr_map, axes[2], r'$\log$ $\Sigma_{\rm SFR}$', args, clabel='', takelog=False, cmap='winter', vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=False, hide_yaxis=False, hide_cbar=False, hide_cbar_ticks=True, cticks_integer=False)

    # ----making SFR dataframe before radial profile plot--------------
    log_sfr_df = pd.DataFrame({'SFR': unp.nominal_values(np.ma.compressed(log_sfr_map)), 'SFR_u': unp.std_devs(np.ma.compressed(log_sfr_map))})

    # ----------getting the distance map--------
    distance_map = get_distance_map(np.shape(log_sfr_map), args)
    log_sfr_df['distance'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, distance_map.data))
    if args.vorbin:
        log_sfr_df['bin_ID'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, args.voronoi_bin_IDs.data))
        log_sfr_df = log_sfr_df.groupby('bin_ID', as_index=False).agg(np.mean)

    # -----plotting SFR radial profile-----------
    axes[3], _ = plot_radial_profile(log_sfr_df, axes[3], args, quant='SFR', ylim=log_sfr_lim, xlim=[0, 5 if args.re_limit is None else args.re_limit], hide_xaxis=False, hide_yaxis=False, hide_cbar=True, short_label=False)
    axes[3].yaxis.set_label_position('right')
    axes[3].yaxis.tick_right()
  
    # ------plotting metallicity vs SFR-----------
    df = pd.DataFrame({'logOH': unp.nominal_values(np.ma.compressed(logOH_map)), 'logOH_u': unp.std_devs(np.ma.compressed(logOH_map)), 'log_sfr':unp.nominal_values(np.ma.compressed(log_sfr_map)), 'log_sfr_u': unp.std_devs(np.ma.compressed(log_sfr_map))})
    if args.vorbin:
        df['bin_ID'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, args.voronoi_bin_IDs.data))
        df = df.groupby('bin_ID', as_index=False).agg(np.mean)
    
    # ----excluding the "turnover" metallicity values from the fit (comment out next 2 lines to avoid this)--------
    #axes[4].scatter(df['log_sfr'], df['logOH'], c='w', lw=1, ec='k', s=20, alpha=1)
    #df = df[~df['logOH'].between(8.3, 8.4)]
    
    axes[4].scatter(df['log_sfr'], df['logOH'], c='grey', s=20, alpha=1)
    axes[4].errorbar(df['log_sfr'], df['logOH'], xerr=df['log_sfr_u'], yerr=df['logOH_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

    # -------Z vs SFR fitting-------------
    fit_color = 'cornflowerblue'
    linefit = odr_fit(df, quant_x='log_sfr', quant_y='logOH')
    axes[4] = plot_fitted_line(axes[4], linefit, df['log_sfr'], fit_color, args, short_label=False, index=0, label=f'Slope = {linefit[0]: .2f}')

    # --------annotating axis--------------
    axes[4].text(0.05, 0.95, r'$\log$(M/M$_{\odot}$) =' + f'{args.log_mass: .2f}', fontsize=args.fontsize, c='k', ha='left', va='top', transform=axes[4].transAxes)
    axes[4].text(0.95, 0.05, f'ID #{objid}', fontsize=args.fontsize, c='k', ha='right', va='bottom', transform=axes[4].transAxes)
    axes[4].set_xlim(log_sfr_lim[0], log_sfr_lim[1]) # kpc
    axes[4].set_ylim(Zlim[0], Zlim[1])
    axes[4] = annotate_axes(axes[4], r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', r'$\log$ (O/H) + 12', args, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)
    axes[4].set_aspect('equal')

    # ----------append fit to dataframe and save-----------
    append_to_sfr_fit_df(field, objid, sfr_int, sfr_sum, linefit, Zdiag, args)

    # -----------saving figure------------
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'metallicity-sfr_radial_profile_{objid:05d}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_radial_profile_fig_multiple(objlist, Zdiag, args, fontsize=10, exclude_ids=[]):
    '''
    Plots and saves a single figure with the 2D metallicity map, metallicity radial profile, SFR map, SFR radial profile and metallicity vs SFR plot for a given list of objects
    '''
    objlist = [item for item in objlist if not item[1] in exclude_ids]
    print(f'Plotting metallicity-SFR figure for {len(objlist)} objects..')

    # --------looping over all objects-------------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        plot_metallicity_sfr_radial_profile_fig_single(objid, field, Zdiag, args, fontsize=fontsize)

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_comparison_fig(objlist, Zdiag_arr, args, Zbranch='low', fontsize=10):
    '''
    Plots and saves a single figure with corner plots for comparisons across a given list of different metallicity diagnostics for a given list of objects
    Repeats that for both high and low metallicity branch solutions
    '''
    args.fontsize = fontsize
    print(f'Plotting metallicity diagnostics comparison for {Zdiag_arr}..')
    args.Zbranch = Zbranch

    # -------setting limits and colors-----------
    Z_limits = [7.0, 9.1]
    color_lim_dict = {'color':[None, None, '', ''], 'bin_ID':[None, None, 'Voronoi bin ID', 'rainbow'], 'distance':[0, 5, 'Galactocentric\ndistance (kpc)', 'cividis'], 'distance_re':[0, 1.5, 'Galactocentric\ndistance ' + r'(R$_e$)', 'cividis'], 'agn_dist':[-1, 1, f'Distance from {args.AGN_diag} SF line', args.diverging_cmap]}
    color = 'brown'

    # --------setting up full figure----------------------------------
    Zbranch_text = '' if np.array([item in ['NB', 'P25'] for item in Zdiag_arr]).all() else f'-{args.Zbranch}'
    nrow, ncol = len(Zdiag_arr) - 1, len(Zdiag_arr) - 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(8, 7), layout='tight')  # layout = 'tight' or 'constrained'
    axes = np.atleast_2d(axes)

    # --------looping over objects---------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args, field=field)
        markersize = 80
        marker='o' if 'Par' in field else 's'

        # ------now plotting for every diag combination--------------
        for col_index in range(ncol):
            for row_index in range(nrow):
                ax = axes[row_index][col_index]
                if row_index < col_index:
                    try: fig.delaxes(ax)
                    except: pass
                else:
                    Z1_index, Z2_index = col_index, row_index + 1
                    Zdiag1 = Zdiag_arr[Z1_index]
                    Zdiag2 = Zdiag_arr[Z2_index]

                    # --------reading in the two dataframes-------
                    df_Zdiag1, logOH_int_Zdiag1, logOH_sum_Zdiag1 = load_metallicity_df(field, objid, Zdiag1, args)
                    df_Zdiag1 = df_Zdiag1.rename(columns={'log_OH': f'log_OH_{Zdiag1}', 'log_OH_u': f'log_OH_{Zdiag1}_err'})
                    df_Zdiag2, logOH_int_Zdiag2, logOH_sum_Zdiag2 = load_metallicity_df(field, objid, Zdiag2, args)
                    df_Zdiag2 = df_Zdiag2.rename(columns={'log_OH': f'log_OH_{Zdiag2}', 'log_OH_u': f'log_OH_{Zdiag2}_err'})
                    df = pd.merge(df_Zdiag1, df_Zdiag2, on = ['distance', 'bin_ID'])
                    df['distance_re'] = df['distance'] / args.re_kpc # converting from kpc to Re
                    df['color'] = color
                
                    if len(df) > 0:
                        # ---plotting integrated--
                        ax.scatter(logOH_sum_Zdiag1.n, logOH_sum_Zdiag2.n, s=markersize, c=color, lw=1, edgecolor='k', marker=marker)
                        ax.errorbar(logOH_sum_Zdiag1.n, logOH_sum_Zdiag2.n, xerr=logOH_sum_Zdiag1.s, yerr=logOH_sum_Zdiag2.s, c='grey', fmt='none', lw=0.5, alpha=0.5)

                        # ---plotting spatially resolved--
                        p = ax.scatter(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], s=markersize/4, c=df[args.colorcol], lw=0, marker=marker, cmap=color_lim_dict[args.colorcol][3], vmin=color_lim_dict[args.colorcol][0], vmax=color_lim_dict[args.colorcol][1])
                        ax.errorbar(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], xerr=df[f'log_OH_{Zdiag1}_err'], yerr=df[f'log_OH_{Zdiag2}_err'], c='grey', fmt='none', lw=0.1, alpha=0.5)

                        ax.plot(Z_limits, Z_limits, ls='dotted', lw=0.1, c='k')

                        if 'NB' in Zdiag2: [x.set_linewidth(2) for x in ax.spines.values()] # making thicker borders for NB axes

                        # ----annotate axis----------
                        ax.set_xlim(Z_limits)
                        ax.set_ylim(Z_limits)
                        ax = annotate_axes(ax, Zdiag1, Zdiag2, args, hide_xaxis=Z2_index < nrow, hide_yaxis=Z1_index)
                    else:
                        print(f'Cannot plot {Zdiag1} vs {Zdiag2} comparison for object {objid}')

    # ------making colorbar----------------------
    if 'p' in locals() and args.colorcol != 'color' and index == len(objlist) - 1:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.45, 0.95, 0.3, 0.02]) # left, bottom, width, height
        cbar = fig.colorbar(p, cax=cbar_ax, orientation='horizontal')
        cbar.set_label(color_lim_dict[args.colorcol][2], fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)

    # -----low/high branch text-------------
    fig.text(0.95, 0.75, f'{Zbranch.title()} metallicity\nbranch chosen', fontsize=args.fontsize, c='k', ha='right', va='top')

    # ------------saving the full figure--------------------------
    colorby_text = f'_colorby_{args.colorcol}' if args.colorcol != 'color' else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'Zdiag{Zbranch_text}_comparison{colorby_text}{args.vorbin_text}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_nb_comparison_sii(objlist, args, fontsize=10):
    '''
    Plots and saves a single figure with comparing NebulaBayes metallicity with and without including SII lines, for a given list of objects
    '''
    args.fontsize = fontsize
    print(f'Plotting NB comparison wth and without SII..')

    # -------setting limits and colors-----------
    Z_limits = [7.0, 9.1]
    color_lim_dict = {'color':[None, None, '', ''], 'bin_ID':[None, None, 'Voronoi bin ID', 'rainbow'], 'distance':[0, 5, 'Galactocentric distance (kpc)', 'cividis'], 'distance_re':[0, 1.5, r'Galactocentric distance (R$_e$)', 'cividis'], 'agn_dist':[-1, 1, f'Distance from {args.AGN_diag} SF line', args.diverging_cmap]}
    color = 'brown'
    counter = 0

    # --------setting up full figure----------------------------------
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.12, top=0.98)
   
    # --------looping over objects---------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args, field=field)
        markersize = 80
        marker='o' if 'Par' in field else 's'

        # --------looping over diagnostics---------------------
        args.exclude_lines = []
        df_wsii, this_logOH_int_wsii, this_logOH_sum_wsii = load_metallicity_df(field, objid, 'NB', args, ensure_sii=True)
        if df_wsii is None: continue
        else: counter += 1
        args.exclude_lines = ['SII']
        df_wosii, this_logOH_int_wosii, this_logOH_sum_wosii = load_metallicity_df(field, objid, 'NB', args)

        df = pd.merge(df_wsii, df_wosii, on = ['distance', 'bin_ID'], suffixes = ('_with_SII', '_without_SII'))
        df['distance_re'] = df['distance'] / args.re_kpc # converting from kpc to Re

        # ---plotting spatially resolved--
        p = ax.scatter(df[f'log_OH_with_SII'], df[f'log_OH_without_SII'], s=markersize/4, c=df[args.colorcol], lw=0, marker=marker, cmap=color_lim_dict[args.colorcol][3], vmin=color_lim_dict[args.colorcol][0], vmax=color_lim_dict[args.colorcol][1])
        ax.errorbar(df[f'log_OH_with_SII'], df[f'log_OH_without_SII'], xerr=df[f'log_OH_u_with_SII'], yerr=df[f'log_OH_u_without_SII'], c='grey', fmt='none', lw=0.1, alpha=0.5)

        ax.plot(Z_limits, Z_limits, ls='dotted', lw=0.1, c='k')

        # ---plotting integrated--
        ax.scatter(this_logOH_sum_wsii.n, this_logOH_sum_wosii.n, s=markersize, c=color, lw=1, edgecolor='k', marker=marker)
        ax.errorbar(this_logOH_sum_wsii.n, this_logOH_sum_wosii.n, xerr=this_logOH_sum_wsii.s, yerr=this_logOH_sum_wosii.s, c='grey', fmt='none', lw=0.5, alpha=0.5)

        # ----annotate axis----------
        ax.set_xlim(Z_limits)
        ax.set_ylim(Z_limits)
        ax = annotate_axes(ax, r'$\log$ (O/H) + 12 (NB: including S II)', r'$\log$ (O/H) + 12 (NB: excluding S II)', args, hide_xaxis=False, hide_yaxis=False)

    print(f'Eventually only plotted only {counter} out of {len(objlist)} objects, that have both with and without SII data')
    
    # ------making colorbar----------------------
    if 'p' in locals() and args.colorcol != 'color':
        cax = inset_axes(ax, width="3%", height="100%", loc='right', bbox_to_anchor=(0.03, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        cbar = fig.colorbar(p, cax=cax, orientation='vertical')
        cbar.set_label(color_lim_dict[args.colorcol][2], fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)

    # ------------saving the full figure--------------------------
    colorby_text = f'_colorby_{args.colorcol}' if args.colorcol != 'color' else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'Zdiag_NB_SII_comparison{colorby_text}{args.vorbin_text}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def get_data_path(field, args):
    '''
    Returns Path object with full path to the data based on a given field name
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'

    return args.root_dir / f'{survey}_data/' / f'{args.version_dict[survey]}' / f'{field}' / 'Products'

# --------------------------------------------------------------------------------------------------------------------
def get_output_path(field, args):
    '''
    Returns Path object with full path to the relevant output based on a given field name
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'

    return args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / f'{field}'

# --------------------------------------------------------------------------------------------------------------------
def get_line_ratio_df(objlist, ratios, args):
    '''
    Creates (or loads from existing) a dataframe for all objects with line fluxes and ratios for a given list of ratios for all spaxels and saves it as a fits file
    Returns dataframe
    '''
    outfitsname = args.root_dir / 'zgrad_paper_plots' / f'line_ratios_df{args.snr_text}{args.only_seg_text}{args.vorbin_text}.fits'
    
    if not os.path.exists(outfitsname) or args.clobber:
        print(f'Making df of {ratios} for {len(objlist)} objects..')

        all_obj_df = pd.DataFrame()
        all_obj_df_int = pd.DataFrame()
        
        # --------looping over all objects-------------------------
        for index2, obj in enumerate(objlist):
            field = obj[0]
            objid = obj[1]
            print(f'\nDoing object {field}-{objid} which is {index2 + 1} of {len(objlist)} objects..')
    
            # ---------loading the data--------------
            full_hdu = load_full_fits(objid, field, args)
            args = load_object_specific_args(full_hdu, args, field=field)

            # -------loop over all ratios---------------
            ratios = np.atleast_1d(ratios)
            all_ratio_df = pd.DataFrame()
            all_ratio_int_arr, all_ratio_sum_arr = [], []

            for index, ratio in enumerate(ratios):
                # ---------getting the ratio maps--------------
                num_lines, den_lines = ratio.split('/')
                num_lines_arr = num_lines.split(',')
                den_lines_arr = den_lines.split(',')
                dummy_map, _, _, _, _ = get_emission_line_map('OII', full_hdu, args, silent=True)
                num_map = np.ma.masked_where(False, np.zeros(np.shape(dummy_map)))
                den_map = np.ma.masked_where(False, np.zeros(np.shape(dummy_map)))
                num_int, num_sum, den_int, den_sum = 0, 0, 0, 0

                # ------getting numerator fluxes--------
                for num_line in num_lines_arr:
                    this_map, _, this_int, this_sum, _ = get_emission_line_map(num_line, full_hdu, args, silent=True)
                    
                    # --------special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back----------
                    if ('OII' in num_lines_arr and 'Hb' in den_lines_arr) or ('OII' in den_lines_arr and 'OIII' in den_lines_arr):
                        if not args.do_not_correct_flux and num_line == 'OIII':
                            ratio_5007_to_4959 = 2.98  # from grizli source code
                            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                            print(f'For ratio {ratio}: un-correcting OIII numerator to include the 4959 component, by factor of {factor:.3f}')
                            this_map = np.ma.masked_where(this_map.mask, this_map.data / factor)
                            this_int /= factor
                            this_sum /= factor
                    
                    num_map = np.ma.masked_where(num_map.mask | this_map.mask, np.sum([num_map.data, this_map.data], axis=0))
                    num_int += this_int
                    num_sum += this_sum

                # ------getting denominator fluxes--------
                for den_line in den_lines_arr:
                    this_map, _, this_int, this_sum, _ = get_emission_line_map(den_line, full_hdu, args, silent=True)

                    # --------special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back----------
                    if ('OII' in den_lines_arr and 'OIII' in den_lines_arr):
                        if not args.do_not_correct_flux and den_line == 'OIII':
                            ratio_5007_to_4959 = 2.98  # from grizli source code
                            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                            print(f'For ratio {ratio}: un-correcting OIII denominator to include the 4959 component, by factor of {factor:.3f}')
                            this_map = np.ma.masked_where(this_map.mask, this_map.data / factor)
                            this_int /= factor
                            this_sum /= factor

                    den_map = np.ma.masked_where(den_map.mask | this_map.mask, np.sum([den_map.data, this_map.data], axis=0))
                    den_int += this_int
                    den_sum += this_sum

                # --------getting the other relevant maps--------------
                net_mask = num_map.mask | den_map.mask
                num_map = np.ma.masked_where(net_mask, num_map.data)
                den_map = np.ma.masked_where(net_mask, den_map.data)
                
                if args.vorbin: vorbin_id_map = np.ma.masked_where(net_mask, args.voronoi_bin_IDs.data)
                else: vorbin_id_map = np.ma.masked_where(net_mask, np.arange(np.shape(distance_map)[0] * np.shape(distance_map)[1]).reshape(np.shape(distance_map)))

                # --------making the dataframe--------------
                df = pd.DataFrame({'bin_ID': np.ma.compressed(vorbin_id_map), \
                                f'{num_lines}': unp.nominal_values(np.ma.compressed(num_map)), f'{num_lines}_u': unp.std_devs(np.ma.compressed(num_map)), \
                                f'{den_lines}': unp.nominal_values(np.ma.compressed(den_map)), f'{den_lines}_u': unp.std_devs(np.ma.compressed(den_map)), \
                                })
                if index == 0: # doing the distance only the first time, as it should be the same every time really
                    distance_map = get_distance_map(np.shape(dummy_map), args)
                    distance_map = np.ma.masked_where(net_mask, distance_map)
                    df['distance'] = np.ma.compressed(distance_map)
                
                df = df.groupby('bin_ID', as_index=False).agg(np.mean)
                df['field'] = field
                df['objid'] = objid

                df = df[(df[f'{num_lines}'] > 0) & (df[f'{den_lines}'] > 0)]
                log_ratio = unp.log10(unp.uarray(df[f'{num_lines}'], df[f'{num_lines}_u']) / unp.uarray(df[f'{den_lines}'], df[f'{den_lines}_u']))
                df[f'log_{ratio}'] = unp.nominal_values(log_ratio)
                df[f'log_{ratio}_u'] = unp.std_devs(log_ratio)

                # -----computing integrated ratios------------
                try: log_ratio_int = unp.log10(num_int / den_int).tolist()
                except: log_ratio_int = ufloat(np.nan, np.nan)
                try: log_ratio_sum = unp.log10(num_sum / den_sum).tolist()
                except: log_ratio_sum = ufloat(np.nan, np.nan)

                # ------appending df and integrated arrays-------
                all_ratio_int_arr.append(log_ratio_int)
                all_ratio_sum_arr.append(log_ratio_sum)
                if index: all_ratio_df = pd.merge(all_ratio_df, df, on=['field', 'objid', 'bin_ID'], how='outer')
                else: all_ratio_df = df
        
            # --------appending df and integrated arrays---------
            all_obj_df = pd.concat([all_obj_df, all_ratio_df])
            
            columns_dict = {'field': field, 'objid':objid}
            columns_dict.update(dict(zip([f'{item}_sum' for item in ratios], unp.nominal_values(all_ratio_sum_arr))))
            columns_dict.update(dict(zip([f'{item}_sum_u' for item in ratios], unp.std_devs(all_ratio_sum_arr))))
            columns_dict.update(dict(zip([f'{item}_int' for item in ratios], unp.nominal_values(all_ratio_int_arr))))
            columns_dict.update(dict(zip([f'{item}_int_u' for item in ratios], unp.std_devs(all_ratio_int_arr))))
            df_int = pd.DataFrame(columns_dict, index=[0])
            all_obj_df_int = pd.concat([all_obj_df_int, df_int])

        # -------saving the fits file-----------
        hdr1 = fits.Header()
        hdr1['field'] = args.field
        hdr1['object'] = args.id
        hdr1['redshift'] = args.z
        primary_hdu = fits.PrimaryHDU(header=hdr1)

        spaxels_df_hdu = fits.BinTableHDU(Table.from_pandas(all_obj_df), name='spaxels')
        int_df_hdu = fits.BinTableHDU(Table.from_pandas(all_obj_df_int), name='int')
        hdul = fits.HDUList([primary_hdu, spaxels_df_hdu, int_df_hdu])
        hdul.writeto(outfitsname, overwrite=True)
        print(f'Saved line ratio dataframes in {outfitsname}')
    
    else:
        print(f'Reading from existing {outfitsname}')
    
    data = fits.open(outfitsname)
    spaxels_df = Table(data['spaxels'].data).to_pandas()
    int_df = Table(data['int'].data).to_pandas()

    return spaxels_df, int_df

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_histogram(full_df_spaxels, objlist, Zdiag_arr, args, fontsize=10, full_df_int=None):
    '''
    Plots and saves the histogram for a given list of ratios for a given list of objects and appropriate dataframes
    Returns a datafrmae
    '''
    args.fontsize = fontsize
    print(f'Plotting histogram of {Zdiag_arr} for {len(objlist)} objects..')
    turnover_ratio_dict = {'OII/Hb':np.nan, 'OIII/Hb':0.87, 'OII,OIII/Hb':1.02, 'OIII/OII':np.nan} # these are the limiting line ratios obtained using Cataldi+25 calibrations
    monosolution_uplim_dict = {'OII/Hb':0.38, 'OIII/Hb':0.67, 'OII,OIII/Hb':0.94, 'OIII/OII':1.31} # these are the limiting line ratios obtained using Cataldi+25 calibrations
    monosolution_lolim_dict = {'OII/Hb':-0.77, 'OIII/Hb':0.29, 'OII,OIII/Hb':0.50, 'OIII/OII':0.38} # these are the limiting line ratios obtained using Cataldi+25 calibrations
    Zdiag_ratios_dict = {'R2':'OII/Hb', 'R3':'OIII/Hb', 'R23':'OII,OIII/Hb', 'O3O2':'OIII/OII'}
    
    if np.array(['/' in item for item in np.atleast_1d(Zdiag_arr)]).any(): ratios = np.atleast_1d(Zdiag_arr)
    else: ratios = [Zdiag_ratios_dict[item] for item in np.atleast_1d(Zdiag_arr)]

    # ----------filtering the dataframes-----------      
    df_base = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})
    if full_df_int is not None: df_int = pd.merge(df_base, full_df_int.dropna(axis=0), on=['field', 'objid'], how='left')
    df_spaxels = pd.merge(df_base, full_df_spaxels.dropna(axis=0), on=['field', 'objid'], how='left')
    
    # ---------setting up the fig--------------
    fig, axes = plt.subplots(1, len(ratios), figsize=(14, 5), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.13, wspace=0.)
    label = 'Voronoi bins' if args.vorbin else 'Pixels'
    axes = np.atleast_1d(axes)

    # ---------looping over ratios----------
    for index, ratio in enumerate(ratios):  
        ax = axes[index]
        
        # ---------plotting the histogram--------------
        if args.histbycol is None:
            average, bins = np.histogram(df_spaxels[f'log_{ratio}'], bins=30)
            histcol_int = 'dummy'
        else:
            if args.histbycol.lower() == 'snr':
                df_spaxels[f'{ratio}_snr'] = 10 ** df_spaxels[f'log_{ratio}'] / 10 ** df_spaxels[f'log_{ratio}_u']
                if full_df_int is not None:
                    df_int[f'{ratio}_snr_sum'] = 10 ** df_int[f'{ratio}_sum'] / 10 ** df_int[f'{ratio}_sum_u']
                    df_int[f'{ratio}_snr_int'] = 10 ** df_int[f'{ratio}_int'] / 10 ** df_int[f'{ratio}_int_u']
                
                histcol, histcol_int = f'{ratio}_snr', f'{ratio}_snr'
            else:
                histcol, histcol_int = args.histbycol, args.histbycol
            
            bin_lolim, bin_uplim = np.min(df_spaxels[f'log_{ratio}']), np.max(df_spaxels[f'log_{ratio}'])
            
            if full_df_int is not None:
                bin_lolim, bin_uplim = min(bin_lolim, np.min(df_int[f'{ratio}_int'] - df_int[f'{ratio}_int_u'])), max(bin_uplim, np.max(df_int[f'{ratio}_int'] + df_int[f'{ratio}_int_u']))
                bin_lolim, bin_uplim = min(bin_lolim, np.min(df_int[f'{ratio}_sum'] - df_int[f'{ratio}_sum_u'])), max(bin_uplim, np.max(df_int[f'{ratio}_sum'] + df_int[f'{ratio}_sum_u']))
            
            bins = np.linspace(bin_lolim, bin_uplim, 30)
            sum_values, _ = np.histogram(df_spaxels[f'log_{ratio}'], bins=bins, weights=df_spaxels[histcol])
            count_entries, _ = np.histogram(df_spaxels[f'log_{ratio}'], bins=bins)
            average = np.zeros_like(sum_values)
            nonzero = count_entries > 0
            average[nonzero] = sum_values[nonzero] / count_entries[nonzero]
        
        bin_centers = 0.5 * (bins[1:] + bins[:-1])
        ax.plot(bin_centers, average, drawstyle='steps-mid', color='cornflowerblue', label=label if index == 0 else None, lw=2, zorder=200) 
        ylim = [-0.2, 10] if args.histbycol is not None and args.histbycol.lower() == 'snr' else [0, ax.get_ylim()[1]]

        # ---------plotting the integrated measurements--------------
        if full_df_int is not None:
            df_int['marker'] = df_int['field'].apply(lambda x: get_marker_type(x))

            if args.histbycol is None or f'{histcol_int}_int' not in df_int:
                yextent = np.diff(ax.get_ylim())[0]
                df_int[f'{histcol_int}_sum'] = ax.get_ylim()[0] + 0.3 * yextent
                df_int[f'{histcol_int}_int'] = df_int[f'{histcol_int}_sum'] + 0.1 * yextent
                
            for index2, m in enumerate(pd.unique(df_int['marker'])):
                df_sub = df_int[df_int['marker'] == m]                
                ax.scatter(df_sub[f'{ratio}_int'], df_sub[f'{histcol_int}_int'], c='navy', marker=m, s=50, ec='k', lw=0.5, label='Grizli-integrated' if index == 0 and index2 == 0 else None)
                ax.scatter(df_sub[f'{ratio}_sum'], df_sub[f'{histcol_int}_sum'], c='brown', marker=m, s=50, ec='k', lw=0.5, label='2D map summed' if index == 0 and index2 == 0 else None)
            
            ax.errorbar(df_int[f'{ratio}_int'], df_int[f'{histcol_int}_int'], xerr = df_int[f'{ratio}_int_u'], color='grey', alpha=0.5, lw=0.5, fmt='none')     
            ax.errorbar(df_int[f'{ratio}_sum'], df_int[f'{histcol_int}_sum'], xerr = df_int[f'{ratio}_sum_u'], color='grey', alpha=0.5, lw=0.5, fmt='none')

        # ---------plotting the different zones--------------
        shade_alpha = 0.2
        patches = []
        if ratio in turnover_ratio_dict:
            if np.isfinite(turnover_ratio_dict[ratio]):
                ax.fill_betweenx([-50, 50], turnover_ratio_dict[ratio], ax.get_xlim()[1], color='salmon', alpha=shade_alpha, lw=0, label='_no_legend', zorder=-10)
                patches.append(matplotlib.patches.Patch(facecolor='salmon', edgecolor='black', linewidth=1, alpha=shade_alpha, label='Outside model\n(no solution)' if index == len(ratios) - 1 else None))
            elif ratio in monosolution_uplim_dict:
                ax.fill_betweenx([-50, 50], monosolution_uplim_dict[ratio], ax.get_xlim()[1], color='salmon', alpha=shade_alpha, lw=0, label='_no_legend', zorder=-10)
                patches.append(matplotlib.patches.Patch(facecolor='salmon', edgecolor='black', linewidth=1, alpha=shade_alpha, label='Outside model\n(no solution)' if index == len(ratios) - 1 else None))

        if ratio in monosolution_uplim_dict:
            if ax.get_xlim()[0] < monosolution_lolim_dict[ratio]: ax.fill_betweenx([-50, 50], ax.get_xlim()[0], monosolution_lolim_dict[ratio], color='salmon', alpha=shade_alpha, lw=0, label='_no_legend_', zorder=-10)
            ax.fill_betweenx([-50, 50],  monosolution_lolim_dict[ratio], monosolution_uplim_dict[ratio], color='limegreen', alpha=shade_alpha, lw=0, label='_no_legend_', zorder=-10)
            patches.append(matplotlib.patches.Patch(facecolor='limegreen', edgecolor='black', linewidth=1, alpha=shade_alpha, label='Only high-Z branch\n(one solution)' if index == len(ratios) - 1 else None))

        ax.set_ylim(ylim)
        ax.set_xlabel(f'Log {ratio.replace(",", "+")}', fontsize=args.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---------annotating the plot--------------
    if len(objlist) == 1: label = f'{len(df_spaxels)} spaxels\nfrom ID #{objlist[0][1]}'
    else: label = f'{len(df_spaxels)} spaxels\nfrom {len(objlist)} objects'
    axes[0].text(0.05 if args.histbycol is None else 0.95, 0.95, label, fontsize=args.fontsize / args.fontfactor, c='k', ha='left' if args.histbycol is None else 'right', va='top', transform=axes[0].transAxes, zorder=250)
    axes[0].set_ylabel('Counts' if args.histbycol is None else f'Average {args.histbycol}', fontsize=args.fontsize)

    axes[0].legend(fontsize=args.fontsize / args.fontfactor, loc='upper right' if args.histbycol is None else 'upper left', framealpha=1)
    axes[-1].legend(handles=patches, fontsize=args.fontsize / args.fontfactor, loc='upper right' if args.histbycol is None else 'upper left', framealpha=1)

    # ---------saving the fig--------------
    histbycol_text = '' if args.histbycol is None else f'_histby_{args.histbycol.lower()}'
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'histogram_ratios_{"_".join(ratios).replace("/", "-")}{histbycol_text}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def bin_2D(value_array, id_array):
    '''
    Bin value_array based on unique bin IDs from id_array
    Return binned array
    '''
    valid_ids = id_array.data[~id_array.mask].astype(int)
    valid_values = value_array[~id_array.mask]

    unique_ids = np.unique(valid_ids)
    means = ndimage.mean(valid_values, labels=valid_ids, index=unique_ids)

    lut = np.full(int(np.nanmax(id_array.data)) + 1, np.nan)
    lut[unique_ids] = means

    result_array = np.full(id_array.shape, np.nan)
    result_array[~id_array.mask] = lut[valid_ids]

    return result_array

# --------------------------------------------------------------------------------------------------------------------
def compute_model(params, psf_array, bin_IDs_map, pixel_scale, quant_y='log_OH'):
    '''
    Compute the model voronoi binned metallicity data from convolving a 2D metallicity profile (built with input gradient and morphology parameters)
    with a PSF array
    '''
    # -------creating an elliptical metallicity map--------------
    center_pix = np.array(np.shape(bin_IDs_map)) / 2.
    y, x = np.indices(np.shape(bin_IDs_map))
    x_rel = x - center_pix[0]
    y_rel = (y - center_pix[1]) / params[2]
    theta = np.radians(params[3])
    x_rel = x_rel * np.cos(theta) + y_rel * np.sin(theta)
    y_rel = -x_rel * np.sin(theta) + y_rel * np.cos(theta)

    distance_map = np.sqrt(x_rel ** 2 + y_rel ** 2) * pixel_scale # arcsec
    logOH_map_model = np.poly1d(params[:2])(distance_map)
    
    # ------smoothing metallicity model with PSF------------------
    logOH_map_model_smoothed = convolve(logOH_map_model, psf_array, boundary='extend')

    # -------voronoi binning the smoothed map as per given vorbin segmentation----------
    logOH_map_model_binned = bin_2D(logOH_map_model_smoothed, bin_IDs_map)
    logOH_df_model = pd.DataFrame({'distance_arcsec': distance_map.flatten(), quant_y + '_model': logOH_map_model_binned.flatten(), 'bin_ID': bin_IDs_map.data.flatten()})
    logOH_df_model = logOH_df_model.dropna().reset_index(drop=True)
    logOH_df_model['bin_ID'] = logOH_df_model['bin_ID'].astype(int)
    logOH_df_model = logOH_df_model.groupby(['bin_ID'], as_index=False).agg('mean')

    return logOH_df_model, distance_map, logOH_map_model, logOH_map_model_smoothed, logOH_map_model_binned

# --------------------------------------------------------------------------------------------------------------------
def log_prior(params, bounds):
    '''
    Uniform priors for running MCMC on elliptical metallicity models
    '''
    for index, this_param in enumerate(params):
        if  this_param < bounds[index][0] or this_param > bounds[index][1]:
            return -np.inf
    
    return 0

# --------------------------------------------------------------------------------------------------------------------
def log_likelihood(params, logOH_df, psf_array, bin_IDs_map, bounds, pixel_scale, quant_y='log_OH'):
    '''
    Likelihood function for running MCMC on elliptical metallicity models
    '''
    lnpr = log_prior(params, bounds)
    if not np.isfinite(lnpr):
        return -np.inf
    
    y = logOH_df[quant_y]
    yerr = logOH_df[quant_y + '_u']

    logOH_df_model, _, _, _, _= compute_model(params, psf_array, bin_IDs_map, pixel_scale, quant_y=quant_y)
    logOH_df_merged = logOH_df.merge(logOH_df_model, on=['bin_ID'])
    model = logOH_df_merged[quant_y + '_model']

    log_like = -0.5 * np.sum((y - model) ** 2 / (yerr ** 2) + np.log(2 * np.pi * yerr ** 2))

    return log_like

# --------------------------------------------------------------------------------------------------------------------
def make_psf(filter, args, supersampling_factor=1):
    '''
    Computes NIRISS PSF in a given filter
    Returns 2D array
    '''
    niriss = webbpsf.NIRISS()
    niriss.filter = filter
    niriss.pixelscale = args.pix_size_arcsec
    if args.re_limit is None: fov_arcsec = args.arcsec_limit * 2
    else: fov_arcsec = args.re_limit * args.re_arcsec * 2
    psf = niriss.calc_psf(fov_arcsec=fov_arcsec + niriss.pixelscale, oversample=supersampling_factor)
    psf_array  = psf[0].data
    
    if np.shape(psf_array)[0] % 2 == 0:
        psf = niriss.calc_psf(fov_arcsec=fov_arcsec, oversample=supersampling_factor)
        psf_array  = psf[0].data
    
    return psf_array

# --------------------------------------------------------------------------------------------------------------------
def convert_slope_unit(params, args, quant_x='distance_arcsec'):
    '''
    Convert the slope (first element of the given parameter list) from dex/arcsecond to dex/kpc or dex/Re
    Return the full parameter list
    '''
    if quant_x != 'distance_arcsec':
        if args.re_limit is None: # distance is in kpc
            params[0] *= cosmo.arcsec_per_kpc_proper(args.z).value # converting from dex/arcsec to dex/kpc
        else: # distance is in Re
            params[0] *= args.re_arcsec # converting from dex/arcsec to dex/Re
    
    return params

# --------------------------------------------------------------------------------------------------------------------
def mcmc_vorbin_fit(logOH_df, args, filter='F150W', quant_x='distance_arcsec', quant_y='log_OH', plot_corner=True, Zdiag=''):
    '''
    Fits the given x and y quantities by taking into account PSF smearing (using MCMC) and the given Voronoi bin segmentation map
    Returns fitted parameters
    '''
    # ----------converting distance column of observed dataframe to arcseconds, if not already in arcseconds----------
    if quant_x != 'distance_arcsec':
        if args.re_limit is None: # distance is in kpc
            logOH_df['distance_arcsec'] = logOH_df[quant_x] * cosmo.arcsec_per_kpc_proper(args.z).value  # converting kpc to arcsec
        else: # distance is in Re
            logOH_df['distance_arcsec'] = logOH_df[quant_x] * args.re_arcsec # converting from Re to arcsec

    # -----------initialising the parameters for MCMC---------------
    bin_IDs_map = args.voronoi_bin_IDs

    labels = ['Zslope', 'Zcen', 'q', 'pa']
    params_init = [-0.5, 8, 0.7, 45] # [slope (in dex/arcsec), Z_cen, q=1-ellipticity, position angle in deg]
    bounds = [[-10, 10], [6, 10], [0, 1], [0, 180]]
    ndim = len(params_init)

    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in Zdiag else ''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'_{args.Zbranch}'
    exclude_text = f'_without_{",".join(args.exclude_lines)}' if len(args.exclude_lines) > 0 and 'NB' in Zdiag else ''
    dered_text = '_dered' if args.dered_in_NB and 'NB' in Zdiag else ''
    AGN_diag_text = f'_AGNdiag_{args.AGN_diag}'if args.AGN_diag != 'None' else ''
 
    mcmc_filename = args.root_dir / 'zgrad_paper_plots' / f'{args.field}_{args.id:05d}_logOH{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}_MCMC_ndim_{ndim}_nwalkers_{args.nwalkers}_niter_{args.niter}_upto_{extent_text}{C25_text}{AGN_diag_text}{exclude_text}{dered_text}.h5py'

    # ----------------checking if saved file exists--------------
    if not os.path.exists(mcmc_filename) or args.clobber_mcmc:
        # -------setup backend-----------
        if os.path.exists(mcmc_filename): os.remove(mcmc_filename)
        backend = emcee.backends.HDFBackend(mcmc_filename)
        backend.reset(args.nwalkers, ndim)
        
        # -------running MCMC-----------
        pos = params_init + 1e-4 * np.random.randn(args.nwalkers, ndim)
        sampler = emcee.EnsembleSampler(args.nwalkers, ndim, log_likelihood, args=(logOH_df, args.psf_array, bin_IDs_map, bounds, args.pix_size_arcsec, quant_y), backend=backend)
        dummy = sampler.run_mcmc(pos, args.niter, progress=True)
    else:
        print(f'\nReading saved MCMC results from {mcmc_filename}')
        sampler = emcee.backends.HDFBackend(mcmc_filename)

    # --------extracting results-------------
    flat_samples = sampler.get_chain(discard=args.ndiscard, thin=15, flat=True)

    # --------getting the final results--------------
    params_llim = [np.percentile(flat_samples[:, item], 16) for item in range(np.shape(flat_samples)[1])]
    params_median = [np.percentile(flat_samples[:, item], 50) for item in range(np.shape(flat_samples)[1])]
    params_ulim = [np.percentile(flat_samples[:, item], 84) for item in range(np.shape(flat_samples)[1])]
    if plot_corner: fig_corner = corner.corner(flat_samples, labels=labels, truths=params_median)
    
    return [params_llim, params_median, params_ulim]

# --------------------------------------------------------------------------------------------------------------------
def lenstronomy_fit_wrap(logOH_df, args, filter='F150W', supersampling_factor=1, Zdiag='NB', quant_x='distance_arcsec', quant_y='log_OH', return_intermediate=False):
    '''
    Wrapper function to read in the metallicity map and direct image map, given a set of args, and then call lenstronomy_fit() to fit the metallicity map accounting for PSF smearing
    Returns fitted parameters
    '''

    # -----------loading the data---------------
    full_hdu = load_full_fits(args.id, args.field, args)
    args = load_object_specific_args(full_hdu, args, field=args.field)
    logOH_map, _, _ = load_metallicity_map(args.field, args.id, Zdiag, args)

    # -----reading in direct image-----------------
    filter_image, exptime = get_direct_image(full_hdu, filter, args)
    center_pix = int(np.shape(filter_image)[0] / 2)
    farthest_pix = int(np.shape(logOH_map)[0] / 2)
    filter_image = filter_image[center_pix - farthest_pix : center_pix + farthest_pix, center_pix - farthest_pix : center_pix + farthest_pix]

    # ------interpolating voronoi binned radial profile------------------
    # distance_array = np.linspace(np.min(logOH_df[quant_x]), np.max(logOH_df[quant_x]), 10)
    # interp_func = interp1d(logOH_df[quant_x], logOH_df[quant_y], assume_sorted=False, fill_value='extrapolate', bounds_error=False)
    # interp_func_u = interp1d(logOH_df[quant_x], logOH_df[quant_y + '_u'], assume_sorted=False, fill_value=0, bounds_error=False)
    # distance_map = get_distance_map(np.shape(filter_image), args, for_distmap=True)
    # distance_map = distance_map * cosmo.arcsec_per_kpc_proper(args.z).value # converting from kpc to arcsec
    # logOH_map_lenstronomy = np.ma.masked_where(False, unp.uarray(interp_func(distance_map), interp_func_u(distance_map)))

    # ------replacing nans in logOH map with zeros------------------
    logOH_map_lenstronomy = logOH_map.copy()
    logOH_map_lenstronomy.data[np.isnan(unp.nominal_values(logOH_map_lenstronomy.data))] = ufloat(0, 1e20)

    # -----calling lenstronomy fitter----------------------
    linefit_lenstronomy = lenstronomy_fit(filter_image, logOH_map_lenstronomy, args.psf_array, exptime=exptime, pixel_scale=args.pix_size_arcsec, supersampling_factor=supersampling_factor) # the output slope is dex/arcsec

    if quant_x == 'distance':
        linefit_lenstronomy[0] *= cosmo.arcsec_per_kpc_proper(args.z).value # converting from dex/arcsec to dex/kpc
        if args.re_limit is not None: linefit_lenstronomy[0] *= args.re_kpc # converting from dex/kpc to dex/Re
    linefit_lenstronomy = unp.uarray(linefit_lenstronomy, np.zeros(len(linefit_lenstronomy)))

    if return_intermediate:
        return linefit_lenstronomy, filter_image, logOH_map, logOH_map_lenstronomy
    else:
        return linefit_lenstronomy

# --------------------------------------------------------------------------------------------------------------------
def lenstronomy_fit(light_map, logOH_map, psf_array, exptime=500, pixel_scale=0.066, supersampling_factor=1):
    '''
    Fits the x and y columns using Lenstronomy
    Returns fitted parameters
    This code is adapted from https://github.com/astrobenji/lenstronomy-metals-notebooks/blob/main/tracer_module_tutorial.ipynb
    '''

    # -----setting up the PSF-----------------
    psf_type   = 'PIXEL'
    kwargs_psf = {'psf_type': psf_type, 'kernel_point_source':psf_array, 'point_source_supersampling_factor':supersampling_factor}
    
    # ---------defining numeric kwargs-----------------------------
    kwargs_numerics = {'supersampling_factor':supersampling_factor, 
                   'point_source_supersampling_factor':supersampling_factor, 
                   'supersampling_convolution': True}

    # ---------defining lens model-----------------------------
    lens_params = [[], [], [], [], []] # Since our model is no lens, use a bunch of empty lists

    # ---------defining source model params-----------------------------
    # Source params -- fixed to true input values; We know that lenstronomy can fit galaxy properties, we are not testing this here.
    fixed_source        = [{'amp': 16, 'R_sersic': 1, 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0}]
    kwargs_source_init  = [{'amp': 16, 'R_sersic': 1, 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0}]
    kwargs_source_sigma = [{'amp': 1,  'R_sersic': 1, 'n_sersic': 1, 'e1': 1, 'e2': 1, 'center_x': 1, 'center_y': 1}]
    kwargs_lower_source = [{'amp': 16, 'R_sersic': 1, 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0}]
    kwargs_upper_source = [{'amp': 16, 'R_sersic': 1, 'n_sersic': 1, 'e1': 0, 'e2': 0, 'center_x': 0, 'center_y': 0}]

    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

    # ---------defining tracer model params-----------------------------
    fixed_tracer_source = [{'center_x': 0., 'center_y': 0}]
    kwargs_tracer_source_init  = [{'center_x': 0., 'center_y': 0, 'amp': 8, 'k': 0}] 
    kwargs_tracer_source_sigma = [{'center_x': 1., 'center_y': 1, 'amp': 0.5, 'k': 0.3}]
    kwargs_lower_tracer_source = [{'center_x': -2, 'center_y': -2, 'amp': 6,   'k': -1}]
    kwargs_upper_tracer_source = [{'center_x': 2, 'center_y': 2, 'amp': 9,   'k': +1}]

    tracer_source_params = [kwargs_tracer_source_init, kwargs_tracer_source_sigma, fixed_tracer_source, kwargs_lower_tracer_source, kwargs_upper_tracer_source]

    # ---------combing all models-----------------------------
    kwargs_params = {'lens_model': lens_params, 'source_model': source_params, 'tracer_source_model': tracer_source_params}           
    kwargs_likelihood = {'source_marg': False, 'tracer_likelihood': True}
    kwargs_model = {'lens_model_list': [], 'source_light_model_list': ['SERSIC_ELLIPSE'],  'tracer_source_model_list': ['LINEAR'], 'tracer_source_band': 0, 'tracer_partition':[[0,0]]}

    # --------generate the coordinate grid and image properties (we only read out the relevant lines we need)-----------
    _, _, ra_at_xy_0, dec_at_xy_0, _, _, Mpix2coord, _ = util.make_grid_with_coordtransform(numPix=np.shape(logOH_map)[0], deltapix=pixel_scale, center_ra=0, center_dec=0, subgrid_res=1, inverse=False)

    # ---------reading in observed light map-----------------------------
    kwargs_image_data = {'background_rms':.005 ,  # rms of background noise; background noise per pixel
                     'exposure_time': exptime,  # exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
                     'ra_at_xy_0': ra_at_xy_0,  # RA at (0,0) pixel
                     'dec_at_xy_0': dec_at_xy_0,  # DEC at (0,0) pixel 
                     'transform_pix2angle': Mpix2coord,  # matrix to translate shift in pixel in shift in relative RA/DEC (2x2 matrix). Make sure it's units are arcseconds or the angular units you want to model.
                     'image_data': light_map   # 2d data vector
                     }

    # ---------reading in observed metallicity map-----------------------------
    kwargs_tracer_data = {'noise_map': unp.std_devs(logOH_map.data),
                      'ra_at_xy_0': ra_at_xy_0,  # RA at (0,0) pixel
                      'dec_at_xy_0': dec_at_xy_0,  # DEC at (0,0) pixel 
                      'transform_pix2angle': Mpix2coord,  # matrix to translate shift in pixel in shift in relative RA/DEC (2x2 matrix). Make sure it's units are arcseconds or the angular units you want to model.
                      'image_data': unp.nominal_values(logOH_map.data)}  # 2d data vector

    # ---------fitting process-----------------------------
    multi_band_list = [[kwargs_image_data, kwargs_psf, kwargs_numerics]]
    kwargs_data_joint = {'multi_band_list': multi_band_list, 
                        'multi_band_type': 'single-band',  # 'multi-linear': every imaging band has independent solutions of the surface brightness, 'joint-linear': there is one joint solution of the linear coefficients demanded across the bands.
                        'tracer_data': [kwargs_tracer_data, kwargs_psf, kwargs_numerics]}
    kwargs_constraints = {'linear_solver': True}

    fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    fitting_kwargs_list = [['PSO',  {'sigma_scale': 1, 'n_particles': 100, 'n_iterations': 500}], ['MCMC', {'n_burn': 100, 'n_run': 1000, 'n_walkers': 100, 'sigma_scale': .1}]]

    chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
    kwargs_result = fitting_seq.best_fit()
    linefit = [kwargs_result['kwargs_tracer_source'][0]['k'], kwargs_result['kwargs_tracer_source'][0]['amp']]
    
    return linefit

# --------------------------------------------------------------------------------------------------------------------
def wls_fit(df_input, quant_x='distance_arcsec', quant_y='log_OH'):
    '''
    Fits the x and y columns using WLS
    Returns fitted parameters
    '''
    df = df_input[[quant_x, quant_y, quant_y + '_u']]
    df = df.dropna(axis=0)
    try:
        ones = np.ones(len(df))
        covariates = np.vstack((ones, df[quant_x])).T
        wls_model_fit = WLS(df[quant_y], covariates, weights=1/(df[quant_y + '_u']**2) if (df[quant_y + '_u'] != 0).any() else 1).fit()
        linefit = [ufloat(wls_model_fit.params['x1'], wls_model_fit.bse['x1']), ufloat(wls_model_fit.params['const'], wls_model_fit.bse['const'])]
    except:
        print(f'WARNING: Could not fit radial profile via WLS, returning nan fit parameters')
        linefit = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]

    return linefit

# --------------------------------------------------------------------------------------------------------------------
def original_fit(df, quant_x='distance_arcsec', quant_y='log_OH'):
    '''
    Fits the x and y columns using np.polyfit
    Returns fitted parameters
    '''
    try:
        linefit, linecov = np.polyfit(df[quant_x], df[quant_y], 1, cov=True, w=1. / (df[quant_y + '_u']) ** 2 if (df[quant_y + '_u'] > 0).any() else None)
        linefit = [ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))]
    except:
        print(f'WARNING: Could not fit radial profile via polyfit, returning nan fit parameters')
        linefit = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]
    
    return linefit

# --------------------------------------------------------------------------------------------------------------------
def odr_fit(df_input, quant_x='distance_arcsec', quant_y='log_OH'):
    '''
    Fits the x and y columns using WLS
    Returns fitted parameters
    '''
    df = df_input[[quant_x, quant_y, quant_y + '_u']]
    recenter_y = False
    # ---------rescaling y-axis for better fit----------
    if recenter_y and 'OH' in quant_y:
        recenter_to = 8.0
        df[quant_y] = df[quant_y] - recenter_to
        print(f'\nDeb3541: in odfit(), recentering {quant_y} to {recenter_to}..')

    df = df.dropna(axis=0)
    try:
        model = Model(lambda param, x: np.poly1d(param)(x))
        data = RealData(df[quant_x], df[quant_y], sy=df[quant_y + '_u'] if quant_y + '_u' in df and np.sum(df[quant_y + '_u']) > 0 else None, sx=df[quant_x + '_u'] if quant_x + '_u' in df and np.sum(df[quant_x + '_u']) > 0 else None)
        odr = ODR(data, model, beta0=[0., 0.])  # Initial guess for slope and intercept
        output = odr.run()
        linefit = [ufloat(output.beta[0], output.sd_beta[0]), ufloat(output.beta[1], output.sd_beta[1])]
        if recenter_y and 'OH' in quant_y: linefit[1] = linefit[1] + recenter_to
    except:
        print(f'WARNING: Could not fit radial profile via ODR, returning nan fit parameters')
        linefit = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]

    return linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fit_tests(objid, field, Zdiag, args, filter='F150W', fontsize=10):
    '''
    Test various radial fitting methods (polyfit, WLS, lenstronomy)
    Saves the figure
    '''
    quant_x = 'distance_arcsec'
    quant_y = 'log_OH'
    args.fontsize = fontsize
    print(f'Plotting radial fit comparisons for ({Zdiag}) for {objid}..')

    # -----------loading the data---------------    
    full_hdu = load_full_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args, field=field)
    logOH_df, _, _ =  load_metallicity_df(field, objid, Zdiag, args)
    logOH_df['distance_arcsec'] = logOH_df['distance'] * cosmo.arcsec_per_kpc_proper(args.z).value # converting from kpc to arcsec
    if args.re_limit is None: arcsec_limit = args.arcsec_limit
    else: arcsec_limit = args.re_limit * args.re_arcsec

    # -------doing the radial fits--------------
    linefit_original = original_fit(logOH_df, quant_x=quant_x, quant_y=quant_y)
    linefit_odr = odr_fit(logOH_df, quant_x=quant_x, quant_y=quant_y)
    # run lenstronomy
    linefit_lenstronomy, filter_image, logOH_map, logOH_map_lenstronomy = lenstronomy_fit_wrap(logOH_df, args, filter=filter, supersampling_factor=1, Zdiag=Zdiag, quant_x=quant_x, quant_y=quant_y, return_intermediate=True)
    
    # run MCMMC
    params_llim, params_median, params_ulim = mcmc_vorbin_fit(logOH_df, args, filter=filter, quant_x=quant_x, quant_y=quant_y, plot_corner=False)

    # convert the unit of best fit slope in order to plot/save
    params_median = convert_slope_unit(params_median, args, quant_x=quant_x) # convert from dex/arcsecond to something else
    linefit_mcmc = unp.uarray(params_median[:2], np.mean([np.array(params_median[:2]) - np.array(params_llim[:2]), np.array(params_ulim[:2]) - np.array(params_median[:2])], axis=0))

    # ---------setting up the figure----------------
    fig, axes = plt.subplots(1, 5, figsize=(16, 2.8))
    fig.subplots_adjust(left=0.05, right=0.99, top=0.96, bottom=0.14, wspace=0.4)
    cmap = 'cividis'

    # -----plotting direct image-----------------
    axes[0], _ = plot_rgb_image(full_hdu, ['F200W', 'F150W', 'F115W'], axes[0], args)

    # -----plotting trimmed direct image---------------
    p = axes[1].imshow(filter_image, extent=(-arcsec_limit, arcsec_limit, -arcsec_limit, arcsec_limit), origin='lower')
    axes[1].scatter(0, 0, marker='x', s=10, c='grey')
    axes[1].set_aspect('auto') 
    axes[1] = annotate_axes(axes[1], 'arcsec', 'arcsec', args, label=f'{filter} image', clabel='', hide_cbar=False, p=p)

    # -------zoom-in annotation------------------
    if args.re_limit is not None:
        connect1 = ConnectionPatch(xyA=(arcsec_limit, arcsec_limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[1]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
        connect2 = ConnectionPatch(xyA=(arcsec_limit, -arcsec_limit), xyB=(axes[1].get_xlim()[0], axes[1].get_ylim()[0]), coordsA='data', coordsB='data', axesA=axes[0], axesB=axes[1], color='r', lw=0.5)
        axes[1].add_artist(connect1)
        axes[1].add_artist(connect2)
        
    # -----plotting 2D metallicity map-----------
    logOH_map_to_plot = np.ma.masked_where(logOH_map.mask, unp.nominal_values(logOH_map.data))
    clim = [np.min(np.ma.compressed(logOH_map_to_plot)), np.max(np.ma.compressed(logOH_map_to_plot))]
    p = axes[2].imshow(logOH_map_to_plot, extent=(-arcsec_limit, arcsec_limit, -arcsec_limit, arcsec_limit), cmap=cmap, origin='lower')
    axes[2].scatter(0, 0, marker='x', s=10, c='grey')
    axes[2].set_aspect('auto') 
    axes[2] = annotate_axes(axes[2], 'arcsec', 'arcsec', args, label='Metallicity', clabel='', hide_cbar=False, p=p, hide_yaxis=True)

    # -----plotting interpolated 2D metallicity map-----------
    logOH_map_to_plot = np.ma.masked_where(logOH_map_lenstronomy.mask, unp.nominal_values(logOH_map_lenstronomy))
    p = axes[3].imshow(logOH_map_to_plot, extent=(-arcsec_limit, arcsec_limit, -arcsec_limit, arcsec_limit), cmap=cmap, origin='lower', vmin=clim[0], vmax=clim[1])
    axes[3].scatter(0, 0, marker='x', s=10, c='grey')
    axes[3].set_aspect('auto') 
    axes[3] = annotate_axes(axes[3], 'arcsec', 'arcsec', args, label='Lenstronomy input metallicity', clabel='', hide_cbar=False, p=p, hide_yaxis=True)

    # ------plotting metallicity radial profile-----------
    radprof_ax = axes[4]
    radprof_ax.scatter(logOH_df[quant_x], logOH_df[quant_y], s=20, c='grey', lw=0, alpha=1)
    radprof_ax.errorbar(logOH_df[quant_x], logOH_df[quant_y], yerr=logOH_df[quant_y + '_u'], fmt='none', c='grey', lw=0.5, alpha=0.2)
    radprof_ax.set_aspect('auto') 
    radprof_ax = annotate_axes(radprof_ax, 'Radius (arcsec)', r'$\log$ O/H + 12', args)

    # -------ploting the radial fits--------------
    xarr = logOH_df[quant_x]
    
    radprof_ax.plot(xarr, np.poly1d(unp.nominal_values(linefit_original))(xarr), color='salmon', lw=1, ls='dashed', label=r'$\nabla$Z$_r$ = ' + f'{linefit_original[0].n:.2f}')
    radprof_ax.plot(xarr, np.poly1d(unp.nominal_values(linefit_odr))(xarr), color='sienna', lw=1, ls='dashed', label=r'$\nabla$Z$_r$ = ' + f'{linefit_odr[0].n:.2f}')
    radprof_ax.plot(xarr, np.poly1d(unp.nominal_values(linefit_lenstronomy))(xarr), color='green', lw=1, ls='dashed', label=r'$\nabla$Z$_r$ = ' + f'{linefit_lenstronomy[0].n:.2f}')
    radprof_ax.plot(xarr, np.poly1d(unp.nominal_values(linefit_mcmc))(xarr), color='cornflowerblue', lw=1, ls='dashed', label=r'$\nabla$Z$_r$ = ' + f'{linefit_mcmc[0].n:.2f}')
    radprof_ax.legend(fontsize=args.fontsize / args.fontfactor, loc='upper left')

    # -----------saving figure------------
    radprof_ax.text(0.95, 0.95, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='right', va='top', transform=radprof_ax.transAxes)
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    figname = f'metallicity_fit_tests_{objid:05d}_upto_{extent_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.colorcol == 'ez_z_phot': args.colorcol = 'color'

    # -----------setting up hard-coded values-------------------
    args.fontsize = 10
    args.fontfactor = 1.3
    args.only_seg = True
    args.plot_ratio_maps = True
    args.exclude_lines = [] #['SII']
    
    args.radbin = True
    args.nbins = 5
    args.use_elliptical_bins = True
    
    #args.vorbin = True
    #args.voronoi_line = 'OII'
    #args.voronoi_snr = 3.
    
    args.re_limit = 2.5
    args.use_C25 = True
    args.dered_in_NB = True

    args.nwalkers = 100 # 100
    args.niter = 5000 # 5000
    args.ndiscard = 100

    primary_Zdiag = 'NB'
    cosmos_name = 'web' # choose between '2020' (i.e. COSMOS2020 catalog) or 'web' (i.e. COSMOSWeb catalog))
    args.version_dict = {'passage': 'v0.5', 'glass': 'orig'}

    args.plot_conditions = 'SNR,photometry'.split(',')
    #args.line_list = 'OII,NeIII-3867,Hb,OIII'.split(',')
    #args.line_list = 'OII,Hb,OIII,Ha'.split(',')
    args.line_list = 'OII,Hb,OIII'.split(',')
    args.SNR_thresh = 2
    
    args.Zdiag = 'R2,R3,R23,NB'.split(',')
    #args.Zdiag += ['O3O2']
    args.colorcol = 'distance' if args.re_limit is None else 'distance_re'
    args.phot_models = 'nb'
    log_mass_lim = [7.5, 10]

    # -------setting up objects to plot--------------
    Par28_objects = [300, 1303, 1849, 2867]
    #Par28_objects = [2727] + Par28_objects
    glass_objects = [1721, 1983, 1991, 1333]
    #glass_objects = [2128] + glass_objects

    passage_objlist = [['Par028', item] for item in Par28_objects]
    glass_objlist = [['glass-a2744', item] for item in glass_objects]
    objlist = passage_objlist + glass_objlist

    # -----------setting up global properties-------------------
    args.snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    args.only_seg_text = '_onlyseg' if args.only_seg else ''
    args.vorbin_text = f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}' if args.vorbin else f'_ellbin_{args.nbins}bins' if args.radbin and args.use_elliptical_bins else f'_radbin_{args.nbins}bins' if args.radbin else ''

    # ------additional variables, for clarity-------------
    all_ratios = ['OII/Hb', 'OIII/Hb', 'OII,OIII/Hb', 'OIII/OII', 'OIII/OII,OIII'] #'NeIII-3867/OII', 
    SEL_Zdiags = [item for item in args.Zdiag if 'NB' not in item]

    # ---------venn diagram plot----------------------
    #plot_passage_venn(['Par028'], args, fontsize=10)
    #plot_glass_venn(args, fontsize=10)

    # ---------photoionisation model plots----------------------
    #plot_photoionisation_model_grid('NeIII/OII', 'OIII/Hb', args, fit_y_envelope=True, fontsize=15)
    #plot_photoionisation_models('OIII/Hb', 'Z', args, fontsize=15)

    # ---------single galaxy plot: example galaxy----------------------
    #plot_galaxy_example_fig(1303, 'Par028', args, fontsize=12, show_log_spectra=False, show_log_map=True)
    
    # ---------single galaxy plot: AGN demarcation----------------------
    #args.AGN_colorby = 'distance'
    #plot_AGN_demarcation_figure_single(1303, 'Par028', args, fontsize=15) # AGN demarcation
    # temp_objlistg = [['glass-a2744', item] for item in [1333, 1983, 2128, 1721, 2273, 1991, 1504, 562, 3015, 34, 1622, 2928, 2492]]
    # temp_objlistp = [['Par028', item] for item  in [235,300,502,668,1302,1303,1634,1675,1798,1849,2038,2171,2180,2706,2727,2797,2867]]
    # objlist = temp_objlistg + temp_objlistp
    
    # --------multi-panel AGN demarcation plots------------------
    #plot_AGN_demarcation_figure_multiple(objlist, args, fontsize=10)#, exclude_ids=[1303])

    # ---------single galaxy plot: Z map and gradient----------------------
    #plot_metallicity_fig_single(1303, 'Par028', primary_Zdiag, args, fontsize=10, do_mcmc=False) # zgrad plot
    #plot_metallicity_fig_single(2128, 'glass-a2744', primary_Zdiag, args, fontsize=10, do_mcmc=False) # zgrad plot
    #plot_metallicity_fit_tests(2867, 'Par028', primary_Zdiag, args, fontsize=10)
    
    # --------to create dataframe of all metallicity quantities including radial fits, etc------------------
    # for Zdiag in args.Zdiag:
    #     if 'NB' in Zdiag: continue
    #     for args.Zbranch in ['low', 'high']:
    #         plot_metallicity_fig_multiple(objlist, Zdiag, args, fontsize=10)

    # --------multi-panel Z map plots------------------
    #plot_metallicity_fig_multiple(objlist, primary_Zdiag, args, fontsize=10, do_mcmc=True)
    #plot_metallicity_fig_multiple(objlist, 'R23', args, fontsize=10, do_mcmc=True)

    # ---------metallicity comparison plots----------------------
    #plot_metallicity_comparison_fig(objlist, args.Zdiag, args, Zbranch='low', fontsize=20)
    #plot_metallicity_comparison_fig(objlist, args.Zdiag, args, Zbranch='high', fontsize=20)

    # ---------single galaxy plot: SFR map and correlation----------------------
    #plot_metallicity_sfr_fig_single(300, 'Par028', primary_Zdiag, args, fontsize=10) # z-sfr plot
    
    # ---------single galaxy plot: SFR-Z radial profile----------------------
    #plot_metallicity_sfr_radial_profile_fig_single(1303, 'Par028', primary_Zdiag, args, fontsize=13)

    # ---------loading master dataframe with only objects in objlist------------
    df = make_master_df(objlist, args, sum=True)
    
    # --------multi-panel SFR-Z plots------------------
    #plot_metallicity_sfr_fig_multiple(objlist, primary_Zdiag, args, fontsize=10)#, exclude_ids=[1303])
    # objlist_ha = [[row['field'], row['objid']] for index, row in df[df['redshift'] < 2.8].iterrows()]
    # #plot_nb_comparison_sii(objlist_ha, args, fontsize=15)

    # ---------metallicity latex table for paper----------------------
    #df_latex = make_latex_table(df, args, sum=True)

    # ---------full population plots----------------------
    #plot_gasmass_comparison(df, args, mass_col='lp_mass', fontsize=15)
    #plot_MEx(df, args, mass_col='lp_mass', fontsize=15)
    #df_agn = plot_AGN_demarcation_figure_integrated(df, args, fontsize=15)
    #plot_SFMS(df, args, mass_col='lp_mass', sfr_col='log_SFR', fontsize=15)
    #plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col='logOH_slope_mcmc_NB', fontsize=15)
    #plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col=['logOH_slope_mcmc_NB', 'logOH_slope_mcmc_R23_low', 'logOH_slope_mcmc_R23_C25_low'], fontsize=15)
    #plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col=['logOH_slope_mcmc_NB', 'logOH_slope_lenstronomy_NB'], fontsize=15)
    #plot_MZsfr(df, args, mass_col='lp_mass', zgrad_col='logZ_logSFR_slope', fontsize=15)
    #plot_MZR(df, args, mass_col='lp_mass', z_col='logOH_sum_NB', colorcol='logOH_slope_mcmc_NB', fontsize=15)
    #plot_Mtmix(df, args, mass_col='lp_mass', ycol='log_t_mix', fontsize=15, colorcol='SFR', mgas_method='my')
    #plot_re_histogram(df, args, fontsize=15)

    # -----------line ratio histograms--------------
    #full_df_spaxels, full_df_int = get_line_ratio_df(objlist, all_ratios, args)
    #plot_line_ratio_histogram(full_df_spaxels, objlist, SEL_Zdiags, args, fontsize=15, full_df_int=full_df_int)
    #plot_line_ratio_histogram(full_df_spaxels, objlist, ['OIII/OII,OIII', 'OIII/OII'], args, fontsize=15, full_df_int=full_df_int)

    # -----------other diagnostics--------------
    #plot_metallicity_sfr_radial_profile_fig_multiple(objlist_ha, primary_Zdiag, args, fontsize=13, exclude_ids=[1303])

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
