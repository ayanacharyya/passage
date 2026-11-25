'''
    Filename: explore_ngdeep.py
    Notes: Explore the official NGDEEP catalog and data products to find bright extended galaxies, for testing, for GRACE
    Author : Ayan
    Created: 23-11-26
    Example: run explore_ngdeep.py --snr_cut_o3 3 --sn_cut_ha 10 --sn_cut_hb 2
             run explore_ngdeep.py --snr_cut_all 3 --sn_cut_ha 50 --log_mass_cut 10
'''

from header import *
from util import *
setup_plot_style()

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_hist(df, col_arr, xlabel, outfilename, args, range=None, bins=20):
    '''
    Plot histograms of a series of columns in the same figure
    Saves figure
    Returns figure handle
    '''
    col_arr = np.atleast_1d(col_arr)
    fig, ax = plt.subplots(1)
    for col in col_arr:
        ax.hist(df[f'{col}'], histtype='step', lw=2, label=f'{col}', range=range, bins=bins)

    ax = annotate_axes(ax, xlabel, 'Frequency', args=args)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.legend(loc='best', fontsize=args.fontsize / args.fontfactor)
    fig = save_fig(fig, args.output_dir / 'plots', outfilename, args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def get_observable_redshift_range(required_lines, lam_min=10000.0, lam_max=19100.0):
    '''
    Compute the redshift window where all required_lines simultaneously
    fall within [lam_min, lam_max] (in Angstroms).
    
    Defaults: Roman WFI G150 grism effective bandpass (1.00 - 1.91 microns).
    Returns redshift range as a tuple
    '''
    # Look up rest wavelengths for the requested lines
    rest_lams = {line: rest_wave_dict[line] * 10 for line in required_lines} # converting from nm to Angstroms 
    
    # The bluest line sets the minimum redshift; the reddest sets the maximum redshift
    bluest_line = min(rest_lams, key=rest_lams.get)
    reddest_line = max(rest_lams, key=rest_lams.get)
    
    z_min = (lam_min / rest_lams[bluest_line]) - 1.0
    z_max = (lam_max / rest_lams[reddest_line]) - 1.0
    
    if z_min > z_max:
        print(rf'No overlap: {required_lines} cannot be observed simultaneously within {lam_min:.0f}-{lam_max:.0f} $\AA$')
        return None
        
    return [z_min, z_max]

# -------------------------------------------------------------------------------------------------------
def make_set(df, condition, label, set_arr, label_arr, colname='par_obj', silent=False):
    '''
    Applies the given condition on given df and appends the ID list into a set and assigns a label
    '''
    id_list = df[condition][colname].values
    if not silent: print(f'{len(id_list)} objects meet the {label} condition')
    set_arr.append(set(id_list))
    label_arr.append(label)

    return set_arr, label_arr

# --------------------------------------------------------------------------------------------------------------------
def plot_venn(df, args, required_lines=None, obs_z_range=None, outfilename='venn_selection.png', silent=False):
    '''
    To plot Venn diagrams with a given df, for whichever selection cuts are active in args
    (redshift observability from required_lines, S/N cuts, EW cuts, mass cuts, etc.)
    Returns intersecting dataframe and the axis handle
    '''
    df = df.copy()
    n_fields = len(pd.unique(df['field'])) if 'field' in df else 1
    if 'par_obj' not in df: df['par_obj'] = df['field'].astype(str) + '-' + df['id'].astype(str)
    if not silent: print(f'\nOut of the total {len(df)} objects in {n_fields} field(s)..\n')

    set_arr = []
    label_arr = []

    # ------add redshift observability range set (from required_lines)-----------
    if required_lines is None:
        required_lines = getattr(args, 'line_list', None) if getattr(args, 'line_list', None) is not None else ['Ha', 'OIII']
    if obs_z_range is None:
        obs_z_range = get_observable_redshift_range(required_lines)

    if obs_z_range is not None:
        condition = df['z'].between(obs_z_range[0], obs_z_range[1])
        set_arr, label_arr = make_set(df, condition, f'{obs_z_range[0]:.2f} < z < {obs_z_range[1]:.2f} ({",".join(required_lines)})', set_arr, label_arr, silent=silent)

    # ---------add line S/N cut sets------------
    sn_cut_map = {
        'OIII': getattr(args, 'sn_cut_o3', getattr(args, 'snr_cut_o3', None)),
        'Ha': getattr(args, 'sn_cut_ha', getattr(args, 'snr_cut_ha', None)),
        'Hb': getattr(args, 'sn_cut_hb', getattr(args, 'snr_cut_hb', None)),
        'OII': getattr(args, 'sn_cut_o2', getattr(args, 'snr_cut_o2', None)),
    }
    sn_cut_all = getattr(args, 'sn_cut_all', getattr(args, 'snr_cut_all', None))

    lines_for_sn = list(dict.fromkeys(list(required_lines) + [line for line, val in sn_cut_map.items() if val is not None]))
    for line in lines_for_sn:
        thresh_list = [val for val in [sn_cut_all, sn_cut_map.get(line)] if val is not None]
        if len(thresh_list) > 0:
            thresh = max(thresh_list)
            try:
                condition = df[f'{line}_sn'] > thresh
                set_arr, label_arr = make_set(df, condition, f'{line} S/N > {thresh}', set_arr, label_arr, silent=silent)
            except:
                if not silent: print(f'Could not apply the {line} S/N condition')
                pass

    # ---------add line EW cut sets (if provided)------------
    ew_cut_all = None #args.EW_thresh
    if ew_cut_all is not None:
        for line in required_lines:
            try:
                condition = (np.isfinite(df[f'{line}_EW'])) & (df[f'{line}_EW'] > ew_cut_all)
                set_arr, label_arr = make_set(df, condition, f'{line} EW > {ew_cut_all}', set_arr, label_arr, silent=silent)
            except:
                if not silent: print(f'Could not apply the {line} EW condition')
                pass

    # ---------add stellar mass cut sets (if provided)------------
    mass_min = getattr(args, 'log_mass_cut', getattr(args, 'mass_cut', getattr(args, 'log_mass_min', None)))
    mass_max = getattr(args, 'log_mass_max', None)
    if mass_min is not None and 'log_mass' in df:
        condition = df['log_mass'] > mass_min
        set_arr, label_arr = make_set(df, condition, rf'$\log M_*/M_\odot > {mass_min}$', set_arr, label_arr, silent=silent)
    if mass_max is not None and 'log_mass' in df:
        condition = df['log_mass'] < mass_max
        set_arr, label_arr = make_set(df, condition, rf'$\log M_*/M_\odot < {mass_max}$', set_arr, label_arr, silent=silent)

    # ---------add optional magnitude / size cuts (if provided)------------
    if getattr(args, 'mag_lim', None) is not None and 'mag' in df:
        condition = df['mag'] <= args.mag_lim
        set_arr, label_arr = make_set(df, condition, f'mag <= {args.mag_lim}', set_arr, label_arr, silent=silent)

    if getattr(args, 'a_thresh', None) is not None and 'a_image' in df:
        condition = df['a_image'] >= args.a_thresh
        set_arr, label_arr = make_set(df, condition, f'a_image >= {args.a_thresh}', set_arr, label_arr, silent=silent)

    # ----------plot the venn diagrams----------
    if len(set_arr) < 2:
        if not silent: print(f'Need at least 2 active criteria to draw a Venn diagram (found {len(set_arr)}); skipping Venn plot..')
        if len(set_arr) == 1:
            intersecting_set = set_arr[0]
            if 'field' in df and 'objid' in df:
                intersecting_par_obj = np.transpose([item.split('-') for item in list(intersecting_set)])
                df_int = pd.DataFrame({'field': intersecting_par_obj[0], 'objid': intersecting_par_obj[1].astype(int)})
                df_int = df.merge(df_int, on=['field', 'objid'], how='inner')
            else:
                intersecting_ids = [int(item.split('-')[1]) for item in list(intersecting_set)]
                df_int = df[df['id'].isin(intersecting_ids)].copy()
            df_int.drop('par_obj', axis=1, inplace=True)
            return df_int, None
        return pd.DataFrame(), None

    if len(set_arr) > 6:
        if not silent: print(f'Warning: Venn diagram supports up to 6 sets; plotting first 6 out of {len(set_arr)} criteria..')
        set_arr = set_arr[:6]
        label_arr = label_arr[:6]

    cmap = 'plasma'
    dataset_dict = dict(zip(label_arr, set_arr))

    # ---------manually calling draw_venn() so as to modify petal labels (for 0 counts)----------
    petal_labels = generate_petal_labels(dataset_dict.values(), fmt="{size}")
    petal_labels = {logic: value if int(value) > 0 else '' for logic, value in petal_labels.items()}
    colors = generate_colors(cmap=cmap, n_colors=len(label_arr))
    
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    ax = draw_venn(petal_labels=petal_labels, dataset_labels=dataset_dict.keys(), hint_hidden=False, colors=colors, figsize=(8, 6), fontsize=args.fontsize / args.fontfactor, legend_loc='upper left', ax=ax)
    save_fig(fig, args.output_dir / 'plots', outfilename, args)

    # ---------printing the conditions for each non-zero petal label----------
    print('\n')
    for index in range(len(petal_labels.keys())):
        if list(petal_labels.values())[index] != '':
            bool_arr = list(petal_labels.keys())[index]
            cond_arr = []
            for index2 in range(len(bool_arr)):
                if int(bool_arr[index2]) == 1:
                    cond_arr.append(label_arr[index2])
            print(f'{int(list(petal_labels.values())[index])} objects in: {" and ".join(cond_arr)}')

    # ----------deriving the dataframe corresponding to the innermost intersection----------
    intersecting_set = set.intersection(*set_arr)
    if len(intersecting_set) > 0:
        if 'field' in df and 'objid' in df:
            intersecting_par_obj = np.transpose([item.split('-') for item in list(intersecting_set)])
            df_int = pd.DataFrame({'field': intersecting_par_obj[0], 'objid': intersecting_par_obj[1]})
            df_int['objid'] = df_int['objid'].astype(int)
            df_int = df.merge(df_int, on=['field', 'objid'], how='inner')
        else:
            intersecting_ids = [int(item.split('-')[1]) for item in list(intersecting_set)]
            df_int = df[df['id'].isin(intersecting_ids)].copy()
        df_int.drop('par_obj', axis=1, inplace=True)
        if 'NUMBER' in df_int: df_int.drop('NUMBER', axis=1, inplace=True)
    else:
        df_int = pd.DataFrame()

    return df_int, ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.5
    required_lines = ['Ha', 'OIII']
    cols_to_print = ['id', 'spec_id', 'z', 'spec_z', 'Ha_sn', 'Ha_EW', 'OIII_sn', 'log_mass']

    # ----------reading catalog----------
    # ---the following data was downloaded from v1.2/catalogs/ in the official NGDEEP data reduction
    # ----link: https://web.corral.tacc.utexas.edu/ngdeepinternal/NIRISSData/Grizli_Extractions/epoch2/
    # ----username: ngdeepHUDF
    # ----pwd: ngdeep4theW!n
    input_catalog_filename = 'master_unicorn_v0.21_wallgrizli_v1.2_wrefit2_rmdups_warcspec_wcigale_addeazycolor.fits'
    input_catalog_path = args.input_dir / 'NGDEEP' / input_catalog_filename
    print(f'Reading in {input_catalog_path}')
    df_input = Table.read(input_catalog_path).to_pandas()

    # --------prepping catalog-----------
    available_lines = ['Ha', 'Hb', 'OIII', 'OII'] # these lines are present in the catalog, see master_catalog.readme in the same directory as above
    for line in available_lines:
        df_input[f'{line}_sn'] = df_input[f'{line}_flux'] / df_input[f'{line}_ferr']
        df_input[f'{line}_EW'] = df_input[f'EW50_{line}']
        df_input[f'{line}_EW_err'] = (df_input[f'EW84_{line}'] - df_input[f'EW16_{line}']) / 2
        df_input[f'{line}_EW_sn'] = df_input[f'{line}_EW'] / df_input[f'{line}_EW_err']

    conditions = (~np.isnan(df_input['spec_z'])) & \
                (df_input['use_flag'] <= 1) & \
                (df_input['phot_dups_flag'] == 0) & \
                (df_input['agn_flag']!=1)

    base_cols_to_extract = ['id', 'spec_id', 'ra', 'dec', 'z', 'z_source', 'spec_z', 'spec_catalog', 'bayes.stellar.m_star', 'bayes.stellar.m_star_err']
    cols_to_extract = np.hstack([base_cols_to_extract, 
                                 np.hstack([[f'{line}_flux', f'{line}_ferr', f'{line}_sn'] for line in available_lines]),
                                 np.hstack([[f'{line}_EW', f'{line}_EW_err', f'{line}_EW_sn'] for line in available_lines]),
                                 ])

    df_usable = df_input[conditions][cols_to_extract]
    df = df_usable.copy()
    print(f'Total {len(df)} usable objects in the catalog..')
    df['id'] = df['id'].astype(int)
    df['field'] = 'NGDEEP'

    quant = unp.log10(unp.uarray(df['bayes.stellar.m_star'], df['bayes.stellar.m_star_err']))
    df['log_mass'] = unp.nominal_values(quant)
    df['log_mass_err'] = unp.std_devs(quant)
    df = df.drop(columns=['bayes.stellar.m_star', 'bayes.stellar.m_star_err'])

    # ---------plotting usable catalog--------------------
    fig_sn_og = plot_hist(df, [f'{line}_sn' for line in required_lines], 'S/N', 'SN_full.png', args, range=[0, 200], bins=20)
    fig_ew_og = plot_hist(df, [f'{line}_EW' for line in required_lines], 'EW', 'EW_full.png', args, range=[0, 1000], bins=20)
    fig_mass_og = plot_hist(df, 'log_mass', r'$\log$ M$_*$/M$_\odot$', 'log_mass_full.png', args)

    # ----------plotting Venn diagram--------------------------
    obs_z_range = get_observable_redshift_range(required_lines)
    df_int, ax_venn = plot_venn(df, args, required_lines=required_lines, obs_z_range=obs_z_range, silent=True)
    print('\n',df_int[cols_to_print])
    '''
    # ----------discarding based on redshift--------------------------
    if obs_z_range is not None:
        nobj0 = len(df)
        df = df[df['z'].between(obs_z_range[0], obs_z_range[1])]
        print(f'\n\tOut of the {nobj0} objects, {len(df)} had {obs_z_range[0]:.2f} < z < {obs_z_range[1]:.2f} for observability of required lines {required_lines}')

    # ----------discarding based on mass--------------------------
    if args.log_mass_cut is not None:
        nobj0 = len(df)
        df = df[df['log_mass'] > args.log_mass_cut]
        print(f'\n\tOut of the {nobj0} objects, {len(df)} had log_mass > {args.log_mass_cut}')

    # -------filtering catalog-------------
    if args.sn_cut_all is not None:
        for line in required_lines:
            df = df[(df[f'{line}_sn'] > args.sn_cut_all)]
        print(f'\tOut of which, {len(df)} objects make the SNR cut {args.sn_cut_all} for ALL {len(required_lines)} required lines..')
    if args.sn_cut_o3 is not None:
        df = df[(df['OIII_sn'] > args.sn_cut_o3)]
        print(f'\tOut of which, {len(df)} objects make the SNR cut {args.sn_cut_o3} for OIII..')
    if args.sn_cut_ha is not None:
        df = df[(df['Ha_sn'] > args.sn_cut_ha)]
        print(f'\tOut of which, {len(df)} objects make the SNR cut {args.sn_cut_ha} for Ha..')
    if args.sn_cut_hb is not None:
        df = df[(df['Hb_sn'] > args.sn_cut_hb)]
        print(f'\tOut of which, {len(df)} objects make the SNR cut {args.sn_cut_hb} for Hb..')

    # ---------plotting filtered catalog--------------------
    # fig_sn_filt = plot_hist(df, [f'{line}_sn' for line in available_lines], 'S/N', 'SN_filtered.png', args, range=[0, 200], bins=20)
    # fig_ew_filt = plot_hist(df, [f'{line}_EW' for line in available_lines], 'EW', 'EW_full.png', args, range=[0, 1000], bins=20)
    # fig_mass_filt = plot_hist(df, 'log_mass', r'$\log$ M$_*$/M$_\odot$', 'log_mass_filtered.png', args)
    '''
    # --------exploring filtered catalog--------------------


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
