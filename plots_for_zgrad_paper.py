'''
    Filename: plots_for_zgrad_paper.py
    Notes: Makes ALL plots related to PASSAGE metallicity gradient paper and saves them separately; some intentional hard-coding and reduction of flexibility in this script, for posterity
    Author : Ayan
    Created: 09-04-25
    Example: run plots_for_zgrad_paper.py --voronoi_line NeIII-3867 --voronoi_snr 4 --colorcol radius --phot_models nb
'''
from header import *
from util import *

from get_field_stats import get_crossmatch_with_cosmos, plot_venn, read_stats_df
from plot_mappings_grid import plot_ratio_grid, plot_ratio_model
from make_passage_plots import break_column_into_uncertainty, plot_SFMS_Popesso22, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14
from make_diagnostic_maps import get_emission_line_map, annotate_PAs, plot_linelist, trim_image, get_EB_V_int, get_voronoi_bin_IDs


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
        args.filters = available_filters_for_field_dict[args.field]
        output_dir = get_output_path(args.field, args)
        df_filename = output_dir / f'{args.field}_all_diag_results.csv'
        thisdf = read_stats_df(df_filename, args)  # read in the stats dataframe
        df = pd.concat([df, thisdf], ignore_index=True)

    df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype( str)  # making a unique combination of field and object id
    df = df.drop_duplicates('par_obj', keep='last')
    df['filters'] = df['field'].map(lambda x: ', '.join(available_filters_for_field_dict[x]))
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
        df = df[quant_names].sort_values(by=quant_names)
        df = df.rename(columns={v: k for k, v in quant_names_dict.items()})

    print(f'Reading in existing grid from {grid_filename}')

    return df

# --------------------------------------------------------------------------------------------------------------------
def plot_photoionisation_model_grid(ratio_x, ratio_y, args, fit_y_envelope=False):
    '''
    Plots and saves the ratio vs ratio parameter space and, optionally, a fit to its envelope, for a given photoionisation model grid
    '''
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
def annotate_ax(ax, xlabel, ylabel, args):
    '''
    Annotates the axis of a given x-y plot
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args):
    '''
    Plots and fits the radial profile from a given 2D image in a given axis
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_1D_spectra(od_hdu, ax, args):
    '''
    Plots the 1D spectra for a given object, in a given axis
    Returns the axis handle
    '''
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
    ax = plot_linelist(ax, fontsize=args.fontsize / args.fontfactor)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def annotate_2d_plot(ax, xlabel, ylabel, args, label='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None):
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
    cmap_arr = ['Blues', 'Greens', 'Reds']

    # -------plot direct image for each filter---------
    for index, filter in enumerate(filters):
        ax = plot_direct_image(full_hdu, filter, ax, args, cmap=cmap_arr[index])
        textcolor = mpl_cm.get_cmap(cmap_arr[index])(0.9)
        ax.text(0.05, 0.98 - index * 0.1, filter, c=textcolor, fontsize=args.fontsize / args.fontfactor, ha='left', va='top', transform=ax.transAxes)

    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.text(0.05, 0.05, f'z={args.z:.2f}', c='k', fontsize=args.fontsize / args.fontfactor, ha='left', va='bottom', transform=ax.transAxes)

    pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
    ax = annotate_PAs(pa_arr, ax, fontsize=args.fontsize / args.fontfactor)

    ax.set_xlim(-args.arcsec_limit, args.arcsec_limit)  # arcsec
    ax.set_ylim(-args.arcsec_limit, args.arcsec_limit)  # arcsec

    if 'segmentation_map' in args: ax.contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=0.5)
    ax = annotate_2d_plot(ax, 'arcsec', 'arcsec', args, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)

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
    ax = annotate_2d_plot(ax, 'arcsec', 'arcsec', args, label=label, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_maps(full_hdu, line_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, hide_cbar=True):
    '''
    Plots the 2D line flux maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''
    for index, ax in enumerate(axes):
        line = line_labels[index]
        line_map, _, _, _ = get_emission_line_map(line, full_hdu, args, dered=True)
        ax = plot_2D_map(line_map, ax, line, args, cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_maps(full_hdu, ratio_labels, axes, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, hide_cbar=True):
    '''
    Plots the 2D line ratio maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''
    for index, ax in enumerate(axes):
        ratio = ratio_labels[index]
        num_line, den_line = ratio.split('/')
        line_map_num, _, _, _ = get_emission_line_map(num_line, full_hdu, args, dered=True)
        line_map_den, _, _, _ = get_emission_line_map(den_line, full_hdu, args, dered=True)

        bad_mask = unp.nominal_values(line_map_den.data) == 0
        line_map_den[bad_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        line_map_den = np.ma.masked_where(line_map_den.mask | bad_mask, line_map_den.data)

        ratio_map = np.ma.masked_where(line_map_num.mask | line_map_den.mask, line_map_num.data / line_map_den.data)
        ax = plot_2D_map(ratio_map, ax, ratio, args, cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=index, hide_cbar=index < len(axes) - 1)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_galaxy_example_fig(objid, field, args):
    '''
    Plots and saves a single figure with the direct image, 1D spectra, emission line maps and emission line ratio maps for a given object
    '''
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
    try: args.EB_V = get_EB_V_int(full_hdu, args, verbose=False)
    except: args.EB_V = 0.

    # ---------------voronoi binning stuff---------------
    if args.vorbin and args.voronoi_line is not None:
        line_map, _, _, _ = get_emission_line_map(args.voronoi_line, full_hdu, args, for_vorbin=True)
        args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, args.voronoi_snr, plot=False, quiet=True, args=args)

    # ----------plotting direct image--------------
    ax_dirimg = plot_rgb_image(full_hdu, ['F115W', 'F150W', 'F200W'], ax_dirimg, args)

    # ----------plotting 1D spectra--------------
    ax_1dspec = plot_1D_spectra(od_hdu, ax_1dspec, args)

    # ----------plotting line flux maps--------------
    axes_line_maps = plot_line_maps(full_hdu, args.line_list, axes_line_maps, args, cmap='RdBu_r', vmin=None, vmax=None, hide_xaxis=args.plot_ratio_maps, hide_cbar=True)

    # ----------plotting line ratio image--------------
    if args.plot_ratio_maps:
        axes_ratio_maps = plot_line_ratio_maps(full_hdu, ratios_to_plot, axes_ratio_maps, args, cmap='GnBu_r', vmin=-1, vmax=1, hide_xaxis=False, hide_cbar=True)

    # ----------annotating and saving figure--------------
    fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
    figname = f'{args.field}_{args.id}_example.png'
    save_fig(fig, figname, args)

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation(full_hdu, args):
    '''
    Plots and saves the spatially resolved AGN demarcation for a given object
    '''

    return

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_map(field, objid, Zdiag, args):
    '''
    Loads the pre-computed 2D metallicity map for a given object
    Returns 2D image
    '''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''

    return logOH_map

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_df(field, objid, Zdiag, args):
    '''
    Loads the pre-computed list of spatially resolved metallicity values for a given object as a dataframe
    Returns dataframe
    '''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and Zdiag == 'NB' else ''

    return logOH_df

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fig(full_hdu, Zdiag, args):
    '''
    Plots and saves a single figure with the direct image, 2D metallicity map and metallicity radial profile for a given object
    '''

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_sfr_fig(full_hdu, Zdiag, args):
    '''
    Plots and saves a single figure with the direct image, 2D metallicity map, metallicity radial profile, SFR map, metallicity vs SFR plot for a given object
    '''

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_comparison_fig(objlist, Zdiag_list, args):
    '''
    Plots and saves a single figure with corner plots for comparisons across a given list of different metallicity diagnostics for a given list of objects
    Repeats that for both high and low metallicity branch solutions
    '''

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
    args.AGN_diag = 'Ne3O2'
    args.exclude_lines = 'SII'
    primary_Zdiag = 'NB'
    cosmos_name = '2020' # choose between '2020' (i.e. COSMOS2020 catalog) or 'web' (i.e. COSMOSWeb catalog))
    args.version_dict = {'passage': 'v0.5', 'glass': 'orig'}

    args.plot_conditions = 'SNR,mass'.split(',')
    args.line_list = 'OII,NeIII-3867,Hb,OIII'.split(',')
    args.Zdiag = 'O3O2,R3,R23,NB'.split(',')
    args.SNR_thresh = 2

    # -------setting up objects to plot--------------
    Par28_objects = [300, 1303, 1634, 1849, 2171, 2727, 2867]
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
    plot_galaxy_example_fig(300, 'Par028', args)
    '''
    # ---------individual galaxy plots: looping over objects----------------------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'Doing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')

        full_hdu = load_full_fits(objid, field, args)
        plot_AGN_demarcation(full_hdu, args)
        plot_metallicity_fig(full_hdu, primary_Zdiag, args)
        plot_metallicity_sfr_fig(full_hdu, primary_Zdiag, args)

    # ---------metallicity comparison plots----------------------
    plot_metallicity_comparison_fig(objlist, args.Zdiag, args)
    '''

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
