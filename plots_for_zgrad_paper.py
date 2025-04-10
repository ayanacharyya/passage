'''
    Filename: plots_for_zgrad_paper.py
    Notes: Makes ALL plots related to PASSAGE metallicity gradient paper and saves them separately; some intentional hard-coding and reduction of flexibility in this script, for posterity
    Author : Ayan
    Created: 09-04-25
    Example: run plots_for_zgrad_paper.py --voronoi_line NeIII-3867 --voronoi_snr 4 --Zdiag O3O2,R3,R23,NB --colorcol radius --phot_models nb
'''
from header import *
from util import *

from get_field_stats import get_crossmatch_with_cosmos, plot_venn, read_stats_df
from plot_mappings_grid import plot_ratio_grid, plot_ratio_model
from make_passage_plots import break_column_into_uncertainty, plot_SFMS_Popesso22, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14, plot_mass_excitation

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
    Zdiag_arr = np.hstack([[item] if item in ['NB', 'P25'] else [item + '_low', item + '_high'] for item in args.Zdiag_arr])

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
    p = ax.scatter(df[mass_col], df[sfr_col], c=df['redshift'], plotnonfinite=True, s=100, lw=1, edgecolor='k', vmin=1.7, vmax=3.1, cmap='viridis')
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

    p = ax.scatter(df[mass_col], df['O3Hb'], c=df['redshift'], s=50, lw=0, cmap='cividis', vmin=1.7, vmax=3.1)
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
    p = ax.scatter(df[mass_col], df[zgrad_col], c=df['redshift'], plotnonfinite=True, s=100, lw=1, edgecolor='k', vmin=1.7, vmax=3.1, cmap='viridis')
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

    return df_latex

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

    return full_hdu

# --------------------------------------------------------------------------------------------------------------------
def load_1d_fits(objid, field, args):
    '''
    Loads the 1D.fits file for a given object and returns the hdulist
    '''

    return oned_hdu

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
def annotate_2d_map(ax, xlabel, ylabel, args):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, args):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_image(full_hdu, filter, ax, args):
    '''
    Plots the direct image for a given filter for a given object, in a given axis
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_rgb_image(full_hdu, filters, ax, args):
    '''
    Plots the direct image as an RGB image, combining the given list of filters, for a given object, in a given axis
    Returns the axis handle
    '''

    return ax


# --------------------------------------------------------------------------------------------------------------------
def plot_1D_spectra(oned_hdu, filter,ax, args):
    '''
    Plots the 1D spectra for a given object, in a given axis
    Returns the axis handle
    '''

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_maps(full_hdu, line_labels, axes, args):
    '''
    Plots the 2D line flux maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_maps(full_hdu, ratio_labels, axes, args):
    '''
    Plots the 2D line ratio maps for a given list of lines for a given object, in a given axis
    Returns the axis handle
    '''

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_galaxy_example_fig(objid, field, args):
    '''
    Plots and saves a single figure with the direct image, 1D spectra, emission line maps and emission line ratio maps for a given object
    '''

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
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and args.Zdiag == 'NB' else ''

    return logOH_map

# --------------------------------------------------------------------------------------------------------------------
def load_metallicity_df(field, objid, Zdiag, args):
    '''
    Loads the pre-computed list of spatially resolved metallicity values for a given object as a dataframe
    Returns dataframe
    '''
    Zbranch_text = '' if Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 and args.Zdiag == 'NB' else ''

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
    args.only_seg = True
    args.vorbin = True
    args.AGN_diag = 'Ne3O2'
    args.exclude_lines = 'SII'
    primary_Zdiag = 'NB'
    cosmos_name = '2020' # choose between '2020' (i.e. COSMOS2020 catalog) or 'web' (i.e. COSMOSWeb catalog))
    args.version_dict = {'passage': 'v0.5', 'glass': 'orig'}

    args.plot_conditions = 'SNR,mass'.split(',')
    args.line_list = 'OIII,Hb,OII,NeIII-3867'.split(',')
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
    args.Zdiag_arr = args.Zdiag.split(',')

    # ---------loading full dataframe for all relevant PASSAGE fields------------
    df_all = load_full_df(passage_objlist, args, cosmos_name=cosmos_name)

    # ---------venn diagram plot----------------------
    #plot_passage_venn(df_all, args)

    # ---------loading master dataframe with only objects in objlist------------
    df = make_master_df(df_all, objlist, args)

    # ---------photoionisation model plots----------------------
    #plot_photoionisation_model_grid('NeIII/OII', 'OIII/Hb', args, fit_y_envelope=True)
    #plot_photoionisation_models('OIII/Hb', 'Z', args)

    # ---------full population plots----------------------
    #plot_SFMS(df, args, mass_col='lp_mass', sfr_col='lp_SFR')
    #plot_MEx(df, args, mass_col='lp_mass')
    #plot_MZgrad(df, args, mass_col='lp_mass', zgrad_col='logOH_slope_NB')

    # ---------metallicity latex table for paper----------------------
    df_latex = make_latex_table(df, args)
    '''
    # ---------individual galaxy plot: example galaxy----------------------
    plot_galaxy_example_fig(300, 'Par028', args)

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
    Zdiag_list = args.Zdiag.split(',')
    plot_metallicity_comparison_fig(objlist, Zdiag_list, args)
    '''

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
