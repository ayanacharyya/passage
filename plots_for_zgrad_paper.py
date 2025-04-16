'''
    Filename: plots_for_zgrad_paper.py
    Notes: Makes ALL plots related to PASSAGE metallicity gradient paper and saves them separately; some intentional hard-coding and reduction of flexibility in this script, for posterity
    Author : Ayan
    Created: 09-04-25
    Example: run plots_for_zgrad_paper.py --do_not_correct_pixel --voronoi_line NeIII-3867 --voronoi_snr 4 --phot_models nb --debug_Zsfr
'''
from header import *
from util import *

from get_field_stats import get_crossmatch_with_cosmos, plot_venn, read_stats_df
from plot_mappings_grid import plot_ratio_grid, plot_ratio_model
from make_passage_plots import break_column_into_uncertainty, plot_SFMS_Popesso22, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14
from make_diagnostic_maps import get_emission_line_map, annotate_PAs, plot_linelist, trim_image, get_EB_V, get_voronoi_bin_IDs, get_AGN_func_methods, AGN_func, take_safe_log_ratio, annotate_BPT_axes, get_distance_map, compute_SFR

plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, figname, args):
    '''
    Saves a given figure handle as a given output filename
    '''
    plot_output_dir = args.root_dir / 'zgrad_paper_plots'
    plot_output_dir.mkdir(exist_ok=True, parents=True)

    figname = plot_output_dir / figname

    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def load_full_df(objlist, args, cosmos_name='web'):
    '''
    Loads and returns a dataframe with ALL objects in the fields listed within objlist
    Returns dataframe
    '''
    fields = np.unique(np.array(objlist)[:,0])
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

    conditions_from_cosmos = ['mass', 'sfr', 'sSFR']
    if len(set(conditions_from_cosmos).intersection(set(args.plot_conditions))) > 0:
        df = get_crossmatch_with_cosmos(df, args, cosmos_name=cosmos_name)

    return df

# --------------------------------------------------------------------------------------------------------------------
def plot_passage_venn(df, args):
    '''
    Plots and saves the Venn diagram given a set of conditions
    '''
    print(f'Plotting Venn diagram..')
    df_int, ax = plot_venn(df, args, silent=True)

    # ----------annotate and save the diagram----------
    fig = ax.figure
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)

    has_fields = [str(int(item[3:])) for item in pd.unique(df['field'])]
    has_fields.sort(key=natural_keys)
    field_text = ','.join(has_fields)

    plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
    plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')

    figname = f'Par{field_text}_venn_diagram_{plot_conditions_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def get_logOH_df(df, args, survey='passage'):
    '''
    Loads and returns a dataframe that holds all the metallicity-related properties for a given list of objects
    Returns dataframe
    '''
    filename = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'logOHgrad_df{args.snr_text}{args.only_seg_text}{args.vorbin_text}.txt'
    df_logOH = pd.read_csv(filename)
    df_logOH = df_logOH.drop_duplicates(subset=['field', 'objid', 'logOH_diagnostic'], keep='last')
    df_logOH = pd.merge(df[['field', 'objid']], df_logOH, on=['field', 'objid'], how='inner')

    logOH_cols = ['logOH_int', 'logOH_slope', 'logOH_cen']
    Zdiag_arr = np.hstack([[item] if item in ['NB', 'P25'] else [item + '_low', item + '_high'] for item in args.Zdiag])

    for Zdiag in Zdiag_arr:
        if 'low' in Zdiag: df_sub = df_logOH[(df_logOH['logOH_diagnostic'] == Zdiag[:-4]) & (df_logOH['logOH_branch'] == 'low')]
        elif 'high' in Zdiag: df_sub = df_logOH[(df_logOH['logOH_diagnostic'] == Zdiag[:-5]) & (df_logOH['logOH_branch'] == 'high')]
        else: df_sub = df_logOH[(df_logOH['logOH_diagnostic'] == Zdiag) & (df_logOH['logOH_branch'] == 'low')]

        df_sub = df_sub.drop(['logOH_diagnostic', 'logOH_branch'], axis=1)
        df_sub = df_sub.rename(columns={'log_OH_int_u':'logOH_int_u'})
        rename_dict = {}
        for col in logOH_cols:
            rename_dict.update({f'{col}': f'{col}_{Zdiag}'})
            rename_dict.update({f'{col}_u': f'{col}_{Zdiag}_u'})
        df_sub = df_sub.rename(columns=rename_dict)
        common_columns = list(set(df.columns) & set(df_sub.columns))
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
def make_master_df(df_all, objlist, args):
    '''
    Creates and returns a dataframe that holds all the relevant global properties for a given list of objects
    Returns dataframe
    '''
    df = pd.DataFrame({'field': np.array(objlist)[:, 0], 'objid': np.array(objlist)[:, 1].astype(int)})
    df_short = pd.merge(df, df_all, on=['field', 'objid'], how='left')
    df_short = df_short.loc[:, ~df_short.columns.duplicated()]

    base_cols_to_extract = ['objid', 'ra', 'dec', 'redshift', 'mag', 'a_image', 'b_image', 'cosmos_id']
    line_cols_to_extract = np.hstack([[f'{line}_int', f'{line}_int_u', f'{line}_EW', f'{line}_EW_u', f'{line}_SNR'] for line in args.line_list])
    measured_cols_to_extract = np.hstack([[f'{item}_int', f'{item}_cen', f'{item}_slope'] for item in ['EB_V', 'SFR']])
    sed_cols_to_extract = np.hstack([[f'lp_{item}', f'lp_{item}_u'] for item in ['mass', 'SFR', 'sSFR']])
    cols_to_extract = np.hstack([base_cols_to_extract, line_cols_to_extract, measured_cols_to_extract, sed_cols_to_extract])

    for col in cols_to_extract:
        if col in df_short: df[col] = df_short[col]
        else: df[col] = np.nan

    df = get_logOH_df(df, args, survey='passage')
    df = get_logOH_df(df, args, survey='glass')

    df['marker'] = df['field'].apply(lambda x: get_marker_type(x))

    print(f'Extracted a dataframe with {len(df)} rows and {len(df.columns)} columns')
    return df

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS(df, args, mass_col='lp_mass', sfr_col='lp_SFR'):
    '''
    Plots and saves the SF-main sequence diagram given a dataframe with list of objects and properties
    '''
    print(f'Plotting SFMS...')
    if 'log_' in sfr_col and sfr_col not in df and sfr_col[4:] in df:
        df = break_column_into_uncertainty(df, sfr_col[4:], make_log=True)

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    # ----------plotting----------
    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub[sfr_col], c=df_sub['redshift'], marker=m, plotnonfinite=True, s=100, lw=1, edgecolor='k', vmin=1.7, vmax=3.1, cmap='viridis')
    if sfr_col + '_u' in df: ax.errorbar(df[mass_col], df[sfr_col], yerr=df[sfr_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[sfr_col], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)

    # ----------making colorbar----------
    cbar = plt.colorbar(p)
    cbar.set_label('Redshift', fontsize=args.fontsize)

    # ----------plotting theoretical diagrams----------
    ax = plot_SFMS_Popesso22(ax, 2.0, color='cornflowerblue')
    ax = plot_SFMS_Shivaei15(ax, color='salmon')
    ax = plot_SFMS_Whitaker14(ax, 2.0, color='yellowgreen')

    # ---------annotate axes and save figure-------
    plt.legend(fontsize=args.fontsize)
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log SFR (M$_{\odot}$/yr)', fontsize=args.fontsize)

    ax.set_xlim(6.5, 11)
    ax.set_ylim(-3, 1)

    figname = f'SFMS_colorby_redshift.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_MEx(df, args, mass_col='lp_mass'):
    '''
    Plots and saves the mass-excitation diagram given a dataframe with list of objects and properties
    '''
    print(f'Plotting mass-excitation diagram...')
    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    OIII_factor = 2.98 / (1 + 2.98)
    df['OIII_int'] *= OIII_factor
    df['OIII_int_u'] *= OIII_factor
    quant =unp.log10(unp.uarray(df['OIII_int'], df['OIII_int_u']) / unp.uarray(df['Hb_int'], df['Hb_int_u']))
    df['O3Hb'] = unp.nominal_values(quant)
    df['O3Hb_u'] = unp.std_devs(quant)

    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub['O3Hb'], c=df_sub['redshift'], s=50, lw=0, cmap='cividis', vmin=1.7, vmax=3.1)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df['O3Hb'], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)

    # ----------making colorbar----------
    cbar = plt.colorbar(p)
    cbar.set_label('Redshift', fontsize=args.fontsize)

    # ---------annotate axes and save figure-------
    plt.legend(fontsize=args.fontsize)
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log O III/H$\beta$', fontsize=args.fontsize)

    ax.set_xlim(6.5, 11)
    ax.set_ylim(-1, 1.5)

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

# --------------------------------------------------------------------------------------------------------------------
def plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col='logOH_slope_NB'):
    '''
    Plots and saves the mass-metallicity gradient plot given a dataframe with list of objects and properties
    '''
    print(f'Plotting MZgrad...')

    # ----------setting up the diagram----------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    # ----------plotting----------
    for m in pd.unique(df['marker']):
        df_sub = df[df['marker'] == m]
        p = ax.scatter(df_sub[mass_col], df_sub[zgrad_col], c=df_sub['redshift'], plotnonfinite=True, s=100, lw=1, edgecolor='k', vmin=1.7, vmax=3.1, cmap='viridis')
    if zgrad_col + '_u' in df: ax.errorbar(df[mass_col], df[zgrad_col], yerr=df[zgrad_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if mass_col + '_u' in df: ax.errorbar(df[mass_col], df[zgrad_col], xerr=df[mass_col + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)

    # ----------making colorbar----------
    cbar = plt.colorbar(p)
    cbar.set_label('Redshift', fontsize=args.fontsize)

    ax.axhline(0, ls='--', c='k', lw=0.5)

    # ---------annotate axes and save figure-------
    plt.legend(fontsize=args.fontsize)
    ax.set_xlabel(r'log M$_*$/M$_{\odot}$', fontsize=args.fontsize)
    ax.set_ylabel(r'log $\nabla$Z$_r$ (dex/kpc)', fontsize=args.fontsize)

    ax.set_xlim(6.5, 11)
    ax.set_ylim(-0.5, 0.5)

    figname = f'MZgrad_colorby_redshift.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def make_latex_table(df, args):
    '''
    Makes and saves a latex table for with metallicity and other columns, given a dataframe with list of objects and properties
    Returns the latex table as pandas dataframe
    '''
    column_dict = {'field':'Field', 'objid':'ID', 'redshift':r'$z$', 'lp_mass':r'$\log$ M$_{\star}$/M$_{\odot}$', 'lp_SFR':r'SFR (M$_{\odot}$/yr)', 'logOH_int_NB':'$\log$ O/H + 12$_{total}$', 'logOH_slope_NB':r'$\nabla Z$ (dex/kpc)'}
    decimal_dict = defaultdict(lambda: 2, redshift=1, lp_mass=1)
    base_cols = ['field', 'objid', 'redshift', 'lp_mass', 'lp_SFR']
    cols_with_errors = ['logOH_int_NB', 'logOH_slope_NB']

    tex_df = pd.DataFrame()
    for thiscol in base_cols + cols_with_errors:
        try:
            if thiscol in cols_with_errors: # special treatment for columns with +/-
                tex_df[thiscol] = [f'{data: .{decimal_dict[thiscol]}} \pm {err: .{decimal_dict[thiscol]}}' for (data, err) in zip(df[thiscol], df[thiscol + '_u'])]
            elif 'logOH_int' in thiscol:
                tex_df[thiscol] = df[thiscol].map(lambda x: '$%.2f$' % x if x < 0 else '$\phantom{-}%.2f$' % x)
            elif thiscol in ['field', 'objid']:
                tex_df[thiscol] = df[thiscol]
            else:
                tex_df[thiscol] = df[thiscol].map(('${:,.' + str(decimal_dict[thiscol]) + 'f}$').format)
        except ValueError: # columns that do not have numbers
            continue

    tex_df = tex_df.rename(columns=column_dict) # change column names to nice ones
    tabname = args.root_dir / 'zgrad_paper_plots' / 'paper_table.tex'
    tex_df.to_latex(tabname, index=None, escape=False)
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
def plot_photoionisation_model_grid(ratio_x, ratio_y, args, fit_y_envelope=False):
    '''
    Plots and saves the ratio vs ratio parameter space and, optionally, a fit to its envelope, for a given photoionisation model grid
    '''
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
    ax.grid(which='both', color='gray', linestyle='solid', linewidth=1, alpha=0.3)
    ax.set_xlabel(f'Log {xratio_name}', fontsize=args.fontsize)
    ax.set_ylabel(f'Log {yratio_name}', fontsize=args.fontsize)
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
        ax.plot(xbins, ybins, c='r', lw=0.5)

        color = 'brown'
        popt, pcov = curve_fit(func, xbins, ybins, p0=p_init)
        ax.plot(xarr, func(xarr, *popt), c=color, lw=2)
        ax.text(0.98, 0.01, f'fit = [{",".join([f"{item:.2f}" for item in popt])}]', c=color, ha='right', va='bottom', fontsize=args.fontsize, transform=ax.transAxes)

    # ---------saving figure--------------
    if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: phot_model_text = 'NebulaBayes'
    figname = f'{phot_model_text}_grid_{yratio_name.replace("/", "-")}_va_{xratio_name.replace("/", "-")}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_photoionisation_models(ratio_y, parameter_x, args):
    '''
    Plots and saves the ratio vs parameter for a given photoionisation model grid and parameter
    '''
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
    ax.set_ylabel(f'Log {ratio_name}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)


    xmin = args.xmin if args.xmin is not None else np.min(np.unique(df[parameter_x])) * 0.9
    xmax = args.xmax if args.xmax is not None else  np.max(np.unique(df[parameter_x])) * 1.1
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
def load_object_specific_args(full_hdu, args):
    '''
    Loads some object specific details into args
    Returns modified args
    '''
    # ----------loading object specific things into args---------------------
    args.z = full_hdu[0].header['REDSHIFT']
    args.id = full_hdu[0].header['ID']
    line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
    args.pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
    args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
    args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines])  # replace 'OIII-5007' with 'OIII'

    # ---------------segmentation map---------------
    segmentation_map = full_hdu['SEG'].data
    args.segmentation_map = trim_image(segmentation_map, args)

    # ---------------dust value---------------
    try: _, args.EB_V, _ = get_EB_V(full_hdu, args, verbose=False, silent=True)
    except: args.EB_V = 0.

    # ---------------voronoi binning stuff---------------
    if args.vorbin and args.voronoi_line is not None:
        line_map, _, _, _, _ = get_emission_line_map(args.voronoi_line, full_hdu, args, for_vorbin=True, silent=True)
        args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, args.voronoi_snr, plot=False, quiet=True, args=args)

    return args

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args, ylim=None, xlim=None):
    '''
    Plots and fits the radial profile from a given 2D image in a given axis
    Returns the axis handle and the linefit
    '''
    # ----------getting the distance map--------
    distance_map = get_distance_map(np.shape(image), args)
    distance_map = np.ma.masked_where(image.mask, distance_map)

    # ----making the dataframe before radial profile plot--------------
    df = pd.DataFrame({'radius': np.ma.compressed(distance_map), 'logOH': unp.nominal_values(np.ma.compressed(image)), 'logOH_u': unp.std_devs(np.ma.compressed(image))})

    # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
    if args.vorbin:
        df['bin_ID'] = np.ma.compressed(np.ma.masked_where(image.mask, args.voronoi_bin_IDs.data))
        df = df.groupby('bin_ID', as_index=False).agg(np.mean)

    if args.radius_max is not None: df = df[df['radius'] <= args.radius_max]
    df = df.sort_values(by='radius').reset_index(drop=True)

    # -------plotting--------
    ax.scatter(df['radius'], df['logOH'], c='grey', s=20, alpha=1)
    ax.errorbar(df['radius'], df['logOH'], yerr=df['logOH_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

    # -------radial fitting-------------
    fit_color = 'salmon'
    linefit, linecov = np.polyfit(df['radius'], df['logOH'], 1, cov=True, w=1. / (df['logOH_u']) ** 2)
    y_fitted = np.poly1d(linefit)(df['radius'])
    ax.plot(df['radius'], y_fitted, color=fit_color, lw=1, ls='dashed')
    linefit = np.array([ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))])
    ax.text(0.05, 0.05, r'$\nabla$Z$_r$' + f' = {linefit[0]: .2f} dex/kpc', c=fit_color, fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)

    # --------annotating axis--------------
    if xlim is not None: ax.set_xlim(xlim[0], xlim[1]) # kpc
    if ylim is not None: ax.set_ylim(ylim[0], ylim[1])
    ax.set_box_aspect(1)
    ax = annotate_axes(ax, 'Radius', 'log (O/H) + 12', args)
    
    return ax, linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_1D_spectra(od_hdu, ax, args):
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

        ax.fill_between(table['rest_wave'], table['norm_flux'] - table['norm_flux_u']/2, table['norm_flux'] + table['norm_flux_u']/2, lw=0, color=color, alpha=0.5, step='pre')#, drawstyle='steps')
        ax.step(table['rest_wave'], table['norm_flux'], lw=1, c=color, alpha=1, where='pre')
        ax.plot(table['rest_wave'], table['norm_cont'], lw=1, c='grey', alpha=1)

    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel(r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)', fontsize=args.fontsize)
    ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    ax = plot_linelist(ax, fontsize=args.fontsize / args.fontfactor, line_list_file=HOME / 'Work/astro/Mappings/labframe.passagelinelist')

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def annotate_axes(ax, xlabel, ylabel, args, label='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''

    ax.text(0.1, 0.9, label, c='k', fontsize=args.fontsize/args.fontfactor, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9), transform=ax.transAxes)

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(xlabel, fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(ylabel, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    if not hide_cbar and p is not None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def annotate_kpc_scale_bar(kpc, ax, args, color='k', loc='lower left'):
    '''
    Annotate existing axis with a scale bar corresponding to a given kpc length
    Returns axis handle
    '''
    pix = kpc * cosmo.arcsec_per_kpc_proper(args.z).value  # converting kpc to arcsec
    scalebar = AnchoredSizeBar(ax.transData, pix, f'{kpc} kpc', loc, pad=0.5, color=color, frameon=False, size_vertical=0.01, fontproperties={'size':args.fontsize / args.fontfactor})
    ax.add_artist(scalebar)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_image(full_hdu, filter, ax, args, cmap='Greens'):
    '''
    Plots the direct image for a given filter for a given object, in a given axis
    Returns the axis handle
    '''
    filter_dummy = 'F140W'
    try:
        hdu = full_hdu['DSCI', filter.upper()]
    except:
        try:
            hdu = full_hdu['DSCI', filter_dummy.upper()]
            print(f'WARNING: Plotting direct image for filter {filter_dummy} instead of {filter}..')
        except:
            sys.exit(f'Neither {filter} nor {filter_dummy} available')

    image = hdu.data
    image = trim_image(image, args)
    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, alpha=1)  # , vmin=0, vmax=0.03)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_rgb_image(full_hdu, filters, ax, args):
    '''
    Plots the direct image as an RGB image, combining the given list of filters, for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting the RGB images with filters {filters}..')
    cmap_arr = ['Blues', 'Greens', 'Reds']

    # -------plot direct image for each filter---------
    for index, filter in enumerate(filters):
        ax = plot_direct_image(full_hdu, filter, ax, args, cmap=cmap_arr[index])
        textcolor = mpl_cm.get_cmap(cmap_arr[index])(0.9)
        ax.text(0.05, 0.98 - index * 0.1, filter, c=textcolor, fontsize=args.fontsize / args.fontfactor, ha='left', va='top', transform=ax.transAxes)

    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.text(0.05, 0.05, f'z={args.z:.2f}', c='k', fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)

    # ----------annotate axis---------------
    pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
    ax = annotate_PAs(pa_arr, ax, fontsize=args.fontsize / args.fontfactor)

    ax.set_xlim(-args.arcsec_limit, args.arcsec_limit)  # arcsec
    ax.set_ylim(-args.arcsec_limit, args.arcsec_limit)  # arcsec

    if 'segmentation_map' in args: ax.contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=0.5)
    ax = annotate_axes(ax, 'arcsec', 'arcsec', args, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)
    ax = annotate_kpc_scale_bar(2, ax, args, loc='lower right')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='cividis', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True):
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

    if 'segmentation_map' in args: ax.contour(args.segmentation_map != args.id, levels=0, colors='w' if args.fortalk else 'k', extent=args.extent, linewidths=0.5) # demarcating the segmentation map zone
    ax = annotate_axes(ax, 'arcsec', 'arcsec', args, label=label, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_maps(full_hdu, ratio_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, hide_cbar=True):
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
        ax = plot_2D_map(ratio_map, ax, ratio, args, cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_line_maps(full_hdu, line_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, hide_cbar=True):
    '''
    Plots the 2D line flux maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''
    print(f'Plotting the line maps {line_labels}..')
    for index, ax in enumerate(axes):
        line = line_labels[index]
        line_map, _, _, _, _ = get_emission_line_map(line, full_hdu, args, dered=True, silent=True)
        ax = plot_2D_map(line_map, ax, line, args, cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_galaxy_example_fig(objid, field, args):
    '''
    Plots and saves a single figure with the direct image, 1D spectra, emission line maps and emission line ratio maps for a given object
    '''
    print(f'\nPlotting example galaxy {field}:{objid}..')
    nrow = 2 # 2 rows: one for direct image + 1D spectra, one for flux maps
    if args.plot_ratio_maps: nrow += 1 # additional row for ratio maps
    ncol = len(args.line_list) # one each for OII, OIII, Hb and NeIII line

    # -------setting up the figure layout-------------
    fig = plt.figure(figsize=(10, 7) if args.plot_ratio_maps else (10, 5))
    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.08, wspace=0.18, hspace=0.3)
    ax_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
    ax_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)

    # ---------setting up emission line map axes----------------
    axes_line_maps = [plt.subplot2grid(shape=(nrow, ncol), loc=(1, item), colspan=1) for item in np.arange(ncol)]
    if args.plot_ratio_maps:
        axes_ratio_maps = [plt.subplot2grid(shape=(nrow, ncol), loc=(2, item), colspan=1) for item in np.arange(ncol)]
        ratios_to_plot = ['OII/NeIII-3867', 'OIII/OII', 'OII/Hb', 'OIII/Hb']

    # ----------loading the full.fits and 1D.fits files--------------
    full_hdu = load_full_fits(objid, field, args)
    od_hdu = load_1d_fits(objid, field, args)
    args = load_object_specific_args(full_hdu, args)

    # ----------plotting direct image--------------
    ax_dirimg = plot_rgb_image(full_hdu, ['F115W', 'F150W', 'F200W'], ax_dirimg, args)

    # ----------plotting 1D spectra--------------
    ax_1dspec = plot_1D_spectra(od_hdu, ax_1dspec, args)

    # ----------plotting line flux maps--------------
    axes_line_maps = plot_line_maps(full_hdu, args.line_list, axes_line_maps, args, cmap='pink', vmin=-20, vmax=-18, hide_xaxis=args.plot_ratio_maps, hide_cbar=True)

    # ----------plotting line ratio image--------------
    if args.plot_ratio_maps:
        axes_ratio_maps = plot_line_ratio_maps(full_hdu, ratios_to_plot, axes_ratio_maps, args, cmap='bone', vmin=-1, vmax=1, hide_xaxis=False, hide_cbar=True)

    # ----------annotating and saving figure--------------
    fig.text(0.05, 0.98, f'{field}: ID {objid}', fontsize=args.fontsize, c='k', ha='left', va='top')
    figname = f'{field}_{objid}_example.png'
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
    df = pd.DataFrame({'xnum': unp.nominal_values(np.atleast_1d(x_num)).flatten(), \
                       'xnum_u': unp.std_devs(np.atleast_1d(x_num)).flatten(), \
                       'xden': unp.nominal_values(np.atleast_1d(x_den)).flatten(), \
                       'xden_u': unp.std_devs(np.atleast_1d(x_den)).flatten(), \
                       'ynum': unp.nominal_values(np.atleast_1d(y_num)).flatten(), \
                       'ynum_u': unp.std_devs(np.atleast_1d(y_num)).flatten(), \
                       'yden': unp.nominal_values(np.atleast_1d(y_den)).flatten(), \
                       'yden_u': unp.std_devs(np.atleast_1d(y_den)).flatten(), \
                       })
    df = df.drop_duplicates().reset_index(drop=True)
    df = df[(df['xnum'] > 0) & (df['xden'] > 0) & (df['ynum'] > 0) & (df['yden'] > 0)]

    y_ratio = unp.log10(unp.uarray(df['ynum'], df['ynum_u']) / unp.uarray(df['yden'], df['yden_u']))
    x_ratio = unp.log10(unp.uarray(df['xnum'], df['xnum_u']) / unp.uarray(df['xden'], df['xden_u']))

    if color is None:
        color = get_distance_map_from_AGN_line(df['xnum'], df['xden'], df['ynum'], df['yden'], args).data
        dist_lim = np.max(np.abs(color))

    p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, cmap=args.diverging_cmap, vmin=-dist_lim, vmax=dist_lim, marker=marker, s=size, lw=lw, edgecolor='w' if args.fortalk else 'k', zorder=10)
    ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c='gray', fmt='none', lw=lw, alpha=0.5)

    return ax, p

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_figure(full_hdu, args, marker='o'):
    '''
    Plots and saves the spatially resolved AGN demarcation for a given object
    '''
    print(f'Plotting AGN demarcation diagram for ID: {full_hdu[0].header["ID"]}..')
    theoretical_lines, line_labels = get_AGN_func_methods(args)
    dist_method = theoretical_lines[0]

    # -----------getting the fluxes------------------
    ynum_map, _, ynum_int, _, _ = get_emission_line_map(args.ynum_line, full_hdu, args, silent=True)
    yden_map, _, yden_int, _, _ = get_emission_line_map(args.yden_line, full_hdu, args, silent=True)

    xnum_map, _, xnum_int, _, _ = get_emission_line_map(args.xnum_line, full_hdu, args, silent=True)
    xden_map, _, xden_int, _, _ = get_emission_line_map(args.xden_line, full_hdu, args, silent=True)

    if not args.do_not_correct_flux and args.AGN_diag in ['H21', 'B22'] and args.xden_line == 'Ha': # special treatment for H-alpha line, in order to add the NII 6584 component back
        factor = 0.823  # from grizli source code
        print(f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
        xden_map = np.ma.masked_where(xden_map.mask, xden_map.data / factor)
        xden_int = xden_int / factor

    # -------setting up the figure--------------------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.1, top=0.98)
    args.fontsize *= args.fontfactor

    # -----------integrated-----------------------
    ax, _ = plot_AGN_demarcation_ax(xnum_int, xden_int, ynum_int, yden_int, ax, args, marker=marker, size=200, lw=2)

    # -----------spatially_resolved-----------------------
    ax, scatter_plot_handle = plot_AGN_demarcation_ax(xnum_map, xden_map, ynum_map, yden_map, ax, args, marker=marker, size=50, lw=0.5)

    # -----------2D map inset-----------------------
    ax_inset = ax.inset_axes([0.7, 0.68, 0.3, 0.3])
    distance_from_AGN_line_map = get_distance_map_from_AGN_line(xnum_map, xden_map, ynum_map, yden_map, args)
    dist_lim = np.max(np.abs(np.ma.compressed(distance_from_AGN_line_map)))
    ax_inset = plot_2D_map(distance_from_AGN_line_map, ax_inset, '', args, takelog=False, cmap=args.diverging_cmap, vmin=-dist_lim, vmax=dist_lim, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

    # -----------annotating axes-----------------------
    ax = annotate_BPT_axes(scatter_plot_handle, ax, args, color_label='Distance from ' + theoretical_lines[0], theoretical_lines=theoretical_lines, line_labels=line_labels)

    # -----------saving figure------------
    ax.text(0.05, 0.9, f'ID #{full_hdu[0].header["ID"]}', fontsize=args.fontsize, c='k', ha='left', va='top', transform=ax.transAxes)
    figname = f'BPT_{full_hdu[0].header["ID"]:05d}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_map(field, objid, Zdiag, args):
    '''
    Loads the pre-computed 2D metallicity map for a given object
    Returns 2D image
    '''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    output_fitsname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'{field}_{objid:05d}_logOH_map{args.snr_text}{args.only_seg_text}{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}_AGNdiag_{args.AGN_diag}{exclude_text}.fits'
    ####################################
    if 'glass' in field:
        output_fitsname = Path(str(output_fitsname).replace('SNR_4.0', 'SNR_2.0'))
        print(f'\nWARNING: Actually choosing metallicity map corresponding to vorbin SNR=2 for {field}-{objid}') ##
    ####################################

    # ------reading in existing outputfile--------------
    print(f'\nReading metallicity map from existing fits file {output_fitsname}')
    hdul = fits.open(output_fitsname)
    hdu = hdul['log_OH']
    logOH_map = np.ma.masked_where(np.isnan(hdu.data), unp.uarray(hdu.data, hdul['log_OH_u'].data))
    if hdu.header['LOG_OH_INT'] is None: logOH_int = ufloat(np.nan, np.nan)
    else: logOH_int = ufloat(hdu.header['LOG_OH_INT'], hdu.header['LOG_OH_INT_ERR'])
    if hdu.header['log_oh_sum'] is None: logOH_sum = ufloat(np.nan, np.nan)
    else: logOH_sum = ufloat(hdu.header['log_oh_sum'], hdu.header['log_oh_sum_err'])

    return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_df(field, objid, Zdiag, args):
    '''
    Loads the pre-computed list of spatially resolved metallicity values for a given object as a dataframe
    Returns dataframe
    '''

    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''
    output_fitsname = args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / 'catalogs' / f'{field}_{objid:05d}_logOH_map{args.snr_text}{args.only_seg_text}{args.vorbin_text}_Zdiag_{Zdiag}{Zbranch_text}_AGNdiag_{args.AGN_diag}{exclude_text}.fits'
    ####################################
    if 'glass' in field:
        output_fitsname = Path(str(output_fitsname).replace('SNR_4.0', 'SNR_2.0'))
        print(f'\nWARNING: Actually choosing metallicity map corresponding to vorbin SNR=2 for {field}-{objid}') ##
    ####################################

    # ------reading in existing outputfile--------------
    print(f'\nReading metallicity map from existing fits file {output_fitsname}')
    hdul = fits.open(output_fitsname)
    logOH_df = Table(hdul['tab'].data).to_pandas()

    hdu = hdul['log_OH']
    if hdu.header['LOG_OH_INT'] is None: logOH_int = ufloat(np.nan, np.nan)
    else: logOH_int = ufloat(hdu.header['LOG_OH_INT'], hdu.header['LOG_OH_INT_ERR'])

    return logOH_df, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fig(full_hdu, field, Zdiag, args):
    '''
    Plots and saves a single figure with the direct image, 2D metallicity map and metallicity radial profile for a given object
    '''
    objid = full_hdu[0].header["ID"]
    print(f'Plotting metallicity ({Zdiag}) figure for {objid}..')

    # -----------loading the data---------------
    logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)

    # --------setting up the figure------------
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.5, hspace=0.3)

    # -----plotting direct image-----------------
    axes[0] = plot_rgb_image(full_hdu, ['F115W', 'F150W', 'F200W'], axes[0], args)

    # -----plotting 2D metallicity map-----------
    Zlim = [7.1, 8.1]
    axes[1] = plot_2D_map(logOH_map, axes[1], f'log O/H + 12 ({Zdiag})', args, takelog=False, cmap='cividis', vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=False, hide_yaxis=True, hide_cbar=False)

    # ------plotting metallicity radial profile-----------
    axes[2], _ = plot_radial_profile(logOH_map, axes[2], args, ylim=Zlim, xlim=[0, 5])
    
    # -----------saving figure------------
    axes[2].text(0.05, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[2].transAxes)
    figname = f'metallicity_{objid:05d}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_fig(full_hdu, field, Zdiag, args):
    '''
    Plots and saves a single figure with the direct image, 2D metallicity map, metallicity radial profile, SFR map, metallicity vs SFR plot for a given object
    '''
    objid = full_hdu[0].header["ID"]
    print(f'Plotting metallicity-SFR figure for {objid}..')

    # -----------loading the data---------------
    logOH_map, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)
    distance = cosmo.comoving_distance(args.z)
    Zlim = [7.1, 8.1]
    log_sfr_lim = [-3, -2]

    # --------setting up the figure------------
    fig, axes = plt.subplots(2 if args.debug_Zsfr else 1, 3, figsize=(9, 5) if args.debug_Zsfr else (9, 3))
    fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.5, hspace=0.3)

    # -------deriving H-alpha map-----------
    N2_plus_Ha_map, _, _, _, _ = get_emission_line_map('Ha', full_hdu, args, silent=True)

    # ------plotting native Ha and SFR maps------------
    if args.debug_Zsfr:
        axes[0][0] = plot_2D_map(N2_plus_Ha_map, axes[0][0], r'H$\alpha$ (orig)', args, takelog=True, cmap='RdPu_r', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=False, hide_cbar=False)
        sfr_map = compute_SFR(N2_plus_Ha_map, distance) # N2_plus_Ha_map here is really Ha_map, because the correction has not been undone yet
        axes[0][1] = plot_2D_map(sfr_map, axes[0][1], f'log SFR (orig)', args, takelog=True, cmap='winter', vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=True, hide_yaxis=True, hide_cbar=False)

    # ----------correcting Ha map------------
    if not args.do_not_correct_flux:
        factor = 0.823 # from James et al. 2023?
        N2_plus_Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask, N2_plus_Ha_map.data / factor)

    log_N2Ha_logOH_poly_coeff = [1, -10] # from approx fit to MAPPINGS models

    logOH_arr = logOH_map.data.flatten()
    log_N2Ha_arr = []
    for logOH in logOH_arr:
        try:
            log_N2Ha = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH)
            log_N2Ha_arr.append(log_N2Ha)
        except:
            log_N2Ha_array.append(np.nan)
    log_N2Ha_map = np.ma.masked_where(logOH_map.mask, np.reshape(log_N2Ha_arr, np.shape(logOH_map)))
    N2Ha_map = np.ma.masked_where(log_N2Ha_map.mask, 10 ** log_N2Ha_map.data)

    Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask | N2Ha_map.mask, N2_plus_Ha_map.data / (1 + N2Ha_map.data))

    if args.debug_Zsfr:
        axes[0][2] = plot_2D_map(Ha_map, axes[0][2], r'H$\alpha$ (corrected)', args, takelog=True, cmap='RdPu_r', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=True, hide_cbar=False)
        axes = axes[1]

    # -------deriving SFR map-----------
    sfr_map = compute_SFR(Ha_map, distance)
    log_sfr_map = np.ma.masked_where(sfr_map.mask, unp.log10(sfr_map.data))

    # -----plotting 2D metallicity map-----------
    axes[0] = plot_2D_map(logOH_map, axes[0], f'log O/H + 12 ({Zdiag})', args, takelog=False, cmap='cividis', vmin=Zlim[0], vmax=Zlim[1], hide_xaxis=False, hide_yaxis=False, hide_cbar=False)

    # -----plotting 2D SFR map-----------
    log_sfr_lim = [-3, -2]
    axes[1] = plot_2D_map(log_sfr_map, axes[1], f'log SFR (corrected)' if args.debug_Zsfr else f'log SFR', args, takelog=False, cmap='winter', vmin=log_sfr_lim[0], vmax=log_sfr_lim[1], hide_xaxis=False, hide_yaxis=True, hide_cbar=False)

    # ------plotting metallicity vs SFR-----------
    df = pd.DataFrame({'logOH': unp.nominal_values(np.ma.compressed(logOH_map)), 'logOH_u': unp.std_devs(np.ma.compressed(logOH_map)), 'log_sfr':unp.nominal_values(np.ma.compressed(log_sfr_map)), 'log_sfr_u': unp.std_devs(np.ma.compressed(log_sfr_map))})
    if args.vorbin:
        df['bin_ID'] = np.ma.compressed(np.ma.masked_where(log_sfr_map.mask, args.voronoi_bin_IDs.data))
        df = df.groupby('bin_ID', as_index=False).agg(np.mean)

    axes[2].scatter(df['log_sfr'], df['logOH'], c='grey', s=20, alpha=1)
    axes[2].errorbar(df['log_sfr'], df['logOH'], xerr=df['log_sfr_u'], yerr=df['logOH_u'], c='grey', fmt='none', lw=0.5, alpha=0.2)

    # -------radial fitting-------------
    fit_color = 'salmon'
    linefit, linecov = np.polyfit(df['log_sfr'], df['logOH'], 1, cov=True, w=1. / (df['logOH_u']) ** 2)
    y_fitted = np.poly1d(linefit)(df['log_sfr'])
    axes[2].plot(df['log_sfr'], y_fitted, color=fit_color, lw=1, ls='dashed')
    linefit = np.array([ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))])
    axes[2].text(0.05, 0.05, f'Slope = {linefit[0]: .2f} dex/kpc', c=fit_color, fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=axes[2].transAxes)

    # --------annotating axis--------------
    axes[2].set_xlim(log_sfr_lim[0], log_sfr_lim[1]) # kpc
    axes[2].set_ylim(Zlim[0], Zlim[1])
    axes[2].set_box_aspect(1)
    axes[2] = annotate_axes(axes[2], 'log SFR (Msun/yr)', 'log (O/H) + 12', args)

    # -----------saving figure------------
    axes[2].text(0.05, 0.9, f'ID #{objid}', fontsize=args.fontsize / args.fontfactor, c='k', ha='left', va='top', transform=axes[2].transAxes)
    debug_text = '_debug' if args.debug_Zsfr else ''
    figname = f'metallicity-sfr_{objid:05d}{debug_text}.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_comparison_fig(objlist, Zdiag_arr, args, Zbranch='low'):
    '''
    Plots and saves a single figure with corner plots for comparisons across a given list of different metallicity diagnostics for a given list of objects
    Repeats that for both high and low metallicity branch solutions
    '''
    print(f'Plotting metallicity diagnostics comparison for {Zdiag_arr}..')
    args.Zbranch = Zbranch

    # -------setting limits and colors-----------
    Z_limits = [7.1, 9.1]
    color_lim_dict = {'color':[None, None, '', ''], 'bin_ID':[None, None, 'Voronoi bin ID', 'rainbow'], 'radius':[0, 5, 'Galactocentric distance (kpc)', 'cividis'], 'agn_dist':[-1, 1, f'Distance from {args.AGN_diag} SF line', args.diverging_cmap]}
    color = 'goldenrod'

    # --------setting up full figure----------------------------------
    Zbranch_text = '' if np.array([item in ['NB', 'P25'] for item in Zdiag_arr]).all() else f'-{args.Zbranch}'
    nrow, ncol = len(Zdiag_arr) - 1, len(Zdiag_arr) - 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(8, 7), layout='tight')  # layout = 'tight' or 'constrained'
    axes = np.atleast_2d(axes)

    # --------looping over objects---------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        markersize = 100
        marker='o' if 'Par' in field else 's'
        log_int_array = []

        # --------looping over diagnostics---------------------
        for index2, Zdiag in enumerate(Zdiag_arr):
            this_df, this_logOH_int = load_metallicity_df(field, objid, Zdiag, args)
            if index2 == 0:
                df = this_df
            else:
                df = pd.merge(df, this_df, on = ['radius', 'bin_ID', 'agn_dist'])
            df = df.rename(columns={'log_OH': f'log_OH_{Zdiag}', 'log_OH_u': f'log_OH_{Zdiag}_err'})
            log_int_array.append(this_logOH_int)

       # ------now plotting for every diag combination--------------
        df['color'] = color
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

                    if f'log_OH_{Zdiag1}' in df and f'log_OH_{Zdiag2}' in df:
                        # ---plotting integrated--
                        ax.scatter(log_int_array[Z1_index].n, log_int_array[Z2_index].n, s=markersize, c=color, lw=1, edgecolor='k', marker=marker)
                        ax.errorbar(log_int_array[Z1_index].n, log_int_array[Z2_index].n, xerr=log_int_array[Z1_index].s, yerr=log_int_array[Z2_index].s, c='grey', fmt='none', lw=0.5, alpha=0.5)

                        # ---plotting spatially resolved--
                        p = ax.scatter(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], s=markersize/4, c=df[args.colorcol], lw=0, marker=marker, cmap=color_lim_dict[args.colorcol][3], vmin=color_lim_dict[args.colorcol][0], vmax=color_lim_dict[args.colorcol][1])
                        ax.errorbar(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], xerr=df[f'log_OH_{Zdiag1}_err'], yerr=df[f'log_OH_{Zdiag2}_err'], c='grey', fmt='none', lw=0.1, alpha=0.5)

                        ax.plot(Z_limits, Z_limits, ls='dotted', lw=0.1, c='k')

                        # ----annotate axis----------
                        if index == len(objlist) - 1:
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

    # ------------saving the full figure--------------------------
    colorby_text = f'_colorby_{args.colorcol}' if args.colorcol != 'color' else ''
    figname = f'Zdiag{Zbranch_text}_comparison{colorby_text}{args.vorbin_text}.png'
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
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.colorcol == 'ez_z_phot': args.colorcol = 'color'

    # -----------setting up hard-coded values-------------------
    args.fontsize = 10
    args.fontfactor = 1.5
    args.only_seg = True
    args.vorbin = True
    args.plot_ratio_maps = True
    args.AGN_diag = 'Ne3O2'
    args.exclude_lines = 'SII'
    primary_Zdiag = 'NB'
    cosmos_name = '2020' # choose between '2020' (i.e. COSMOS2020 catalog) or 'web' (i.e. COSMOSWeb catalog))
    args.version_dict = {'passage': 'v0.5', 'glass': 'orig'}

    args.plot_conditions = 'SNR,mass'.split(',')
    args.line_list = 'OII,NeIII-3867,Hb,OIII'.split(',')
    args.SNR_thresh = 2
    args.Zdiag = 'O3O2,R2,R3,R23,NB'.split(',')
    args.colorcol = 'radius'

    # -------setting up objects to plot--------------
    Par28_objects = [300, 1303, 1849, 2171, 2727, 2867]
    glass_objects = [1333, 1721, 1983, 2128]

    passage_objlist = [['Par028', item] for item in Par28_objects]
    glass_objlist = [['glass-a2744', item] for item in glass_objects]
    objlist = passage_objlist + glass_objlist

    # -----------setting up global properties-------------------
    args.snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    args.only_seg_text = '_onlyseg' if args.only_seg else ''
    args.vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'

    # ---------loading full dataframe for all relevant PASSAGE fields------------
    #df_all = load_full_df(passage_objlist, args, cosmos_name=cosmos_name)

    # ---------venn diagram plot----------------------
    #plot_passage_venn(df_all, args)

    # ---------loading master dataframe with only objects in objlist------------
    #df = make_master_df(df_all, objlist, args)

    # ---------photoionisation model plots----------------------
    #plot_photoionisation_model_grid('NeIII/OII', 'OIII/Hb', args, fit_y_envelope=True)
    #plot_photoionisation_models('OIII/Hb', 'Z', args)

    # ---------full population plots----------------------
    #plot_SFMS(df, args, mass_col='lp_mass', sfr_col='lp_SFR')
    #plot_MEx(df, args, mass_col='lp_mass')
    #plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col='logOH_slope_NB')

    # ---------metallicity latex table for paper----------------------
    #df_latex = make_latex_table(df, args)
    
    # ---------individual galaxy plot: example galaxy----------------------
    #plot_galaxy_example_fig(1303, 'Par028', args)
    '''
    # ---------individual galaxy plots: looping over objects----------------------
    objlist = objlist[:1] ##
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'Doing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')

        full_hdu = load_full_fits(objid, field, args)
        args = load_object_specific_args(full_hdu, args)

        #plot_AGN_demarcation_figure(full_hdu, args, marker='o' if 'Par' in field else 's')
        #plot_metallicity_fig(full_hdu, field, primary_Zdiag, args, Zbranch=args.Zbranch)
        try:
            plot_metallicity_sfr_fig(full_hdu, field, primary_Zdiag, args)
        except:
            print(f'Could not plot SFR figure for ID# {objid}, so skipping this object..')
            pass
    '''
    # ---------metallicity comparison plots----------------------
    plot_metallicity_comparison_fig(objlist, args.Zdiag, args, Zbranch='low')
    plot_metallicity_comparison_fig(objlist, args.Zdiag, args, Zbranch='high')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
