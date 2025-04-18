'''
    Filename: integrated_vs_summed.py
    Notes: Comparing grizli-reported integrated line fluxes, ratios, metallicities etc. vs those derived from summing the 2D emission line maps
    Author : Ayan
    Created: 09-04-25
    Example: run integrated_vs_summed.py --do_not_correct_pixel --Zbranch low
'''
from header import *
from util import *

from plot_mappings_grid import plot_ratio_grid, plot_ratio_model
from make_diagnostic_maps import get_emission_line_maps
from plots_for_zgrad_paper import load_full_fits, load_metallicity_map, get_photoionisation_model_grid, load_object_specific_args

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_photoionisation_models(ratio_y, parameter_x, args):
    '''
    Plots the ratio vs parameter for a given photoionisation model grid and parameter
    Returns axis handle
    '''
    print(f'Plotting photoionisation models for {ratio_y} vs {parameter_x}..')
    if args.phot_models.lower() in ['mappings', 'map']: line_label_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,9', 'NeIII':'[NeIII]3869', 'NeIII-3867':'[NeIII]3869', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31'}) # to map between user input line labels and line labels used in ratio_list.txt file
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})
    args.ynum_line, args.yden_line = ratio_y.split('/')
    if args.ynum_line == 'NeIII-3867': args.ynum_line = 'NeIII'
    if args.yden_line == 'NeIII-3867': args.yden_line = 'NeIII'
    df = get_photoionisation_model_grid(args)

    ratio_y = f'{line_label_dict[args.ynum_line]}/{line_label_dict[args.yden_line]}'
    if ratio_y not in df: df[ratio_y] = df[line_label_dict[args.ynum_line]] / df[line_label_dict[args.yden_line]]

    # ------declare the figure-----------------------
    fig, ax = plt.subplots(figsize=(14, 4))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.15, top=0.95, wspace=0.2)

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

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-2, 2)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def overplot_on_photoionisation_models(df, ratio_list, parameter_x, args, col_int='sandybrown', col_sum='cornflowerblue'):
    '''
    Plots and saves photoionisation models overplotted with given observed ratios from a dataframe
    Returns figure handle
    '''
    fig_arr = []
    # ---------------looping over different ratios-------------
    for ratio in ratio_list:
        ratio_name_int = 'log_' + ratio.replace('/', '_') + '_int'
        ratio_name_sum = 'log_' + ratio.replace('/', '_') + '_sum'

        ax = plot_photoionisation_models(ratio, parameter_x, args)
        for index, row in df.iterrows():
            ax.axhline(row[ratio_name_int], ls='solid', c=col_int, lw=0.5, alpha=1)
            ax.fill_between([ax.get_xlim()[0], ax.get_xlim()[1]], row[ratio_name_int] - row[ratio_name_int + '_u'] / 2, row[ratio_name_int] + row[ratio_name_int + '_u'] / 2, lw=0, color=col_int, alpha=0.1)

            ax.axhline(row[ratio_name_sum], ls='dashed', c=col_sum, lw=0.5, alpha=1)
            ax.fill_between([ax.get_xlim()[0], ax.get_xlim()[1]], row[ratio_name_sum] - row[ratio_name_sum + '_u'] / 2, row[ratio_name_sum] + row[ratio_name_sum + '_u'] / 2, lw=0, color=col_sum, alpha=0.1)

        ax.text(0.98, 0.05, 'Integrated (Grizli) ratio', c=col_int, fontsize=args.fontsize, va='bottom', ha='right', transform=ax.transAxes)
        ax.text(0.98, 0.15, 'Summed ratio', c=col_sum, fontsize=args.fontsize, va='bottom', ha='right', transform=ax.transAxes)

        if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
        elif args.phot_models.lower() in ['nebulabayes', 'nb', 'neb']: phot_model_text = 'NebulaBayes'
        figname = args.output_dir / 'plots' / f'{phot_model_text}_model_{ratio.replace("/", "-")}_vs_Z_withobs.png'

        fig = ax.figure
        fig.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure as {figname}')
        plt.show(block=False)

        fig_arr.append(fig)

    return fig_arr

# --------------------------------------------------------------------------------------------------------------------
def compare_line_fluxes(df, line_list, args, lim=[-17, -15.5]):
    '''
    Plots and saves the integrated vs summed line fluxes from a given dataframe
    Returns the figure handle
    '''
    # ------setting up the figure-------------
    fig, axes = plt.subplots(1, len(ratio_list), figsize=(14, 4), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.15, top=0.95, wspace=0.2)

    # -----------plotting the line ratios-------------
    for line, ax in zip(line_list, axes):
        line_name_int = f'log_{line}_int'
        line_name_sum = f'log_{line}_sum'

        for marker in pd.unique(df['marker']):
            df_sub = df[df['marker'] == marker]
            ax.scatter(df_sub[line_name_int], df_sub[line_name_sum], c=df_sub['color'], marker=marker, s=50, edgecolor='k', lw=0.5, cmap=args.cmap, vmin=args.clim[0], vmax=args.clim[1])
        ax.errorbar(df[line_name_int], df[line_name_sum], xerr=df[line_name_int + '_u'], yerr=df[line_name_sum + '_u'], fmt='none', c='gray', lw=0.5, zorder=-2)

        # --------annotating the figure-----------
        ax.set_xlabel(f'log {line} (int)', fontsize=args.fontsize)
        ax.set_ylabel(f'log {line} (sum)', fontsize=args.fontsize)

        ax.plot([lim[0], lim[1]], [lim[0], lim[1]], c='k', ls='dotted', lw=1)

        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[0], lim[1])

    # --------making colorbar if needed-------------
    if not args.nocolorcoding and args.colorcol != 'ez_z_phot':
        cax = fig.add_axes([0.8, 0.27, 0.15, 0.01])
        cbar = matplotlib.colorbar.ColorbarBase(cax, orientation='horizontal', cmap=args.cmap, norm=mplcolors.Normalize(args.clim[0], args.clim[1]), label=args.colorcol)

    # ----------saving the figure-----------------
    figname = args.output_dir / 'plots' / f'sum_vs_int_line_flux_comparison_{args.snr_text}{args.only_seg_text}_{",".join(line_list)}{args.colorby_text}.png'
    ax.figure.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def compare_line_ratios(df, ratio_list, args, lim=[-2, 2]):
    '''
    Plots and saves the integrated vs summed line ratio from a given dataframe
    Returns the figure handle
    '''
    # -------------photoionisation models-----------
    if args.phot_models.lower() in ['mappings', 'map']:
        line_label_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,9', 'NeIII':'[NeIII]3869', 'NeIII-3867':'[NeIII]3869', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31'}) # to map between user input line labels and line labels used in ratio_list.txt file
        phot_model_text = 'Mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']:
        line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})
        phot_model_text = 'NebulaBayes'

    df_model = get_photoionisation_model_grid(args)
    model_color = 'salmon'

    # ------setting up the figure-------------
    fig, axes = plt.subplots(1, len(ratio_list), figsize=(14, 4), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.15, top=0.95, wspace=0.2)

    # -----------plotting the line ratios-------------
    for ratio, ax in zip(ratio_list, axes):

        # -------------gettinug max and min model ratios-----------
        ynum_line, yden_line = ratio.split('/')
        if ynum_line == 'NeIII-3867': ynum_line = 'NeIII'
        if yden_line == 'NeIII-3867': yden_line = 'NeIII'
        ratio_name = f'{line_label_dict[ynum_line]}/{line_label_dict[yden_line]}'
        if ratio_name not in df_model: df_model[ratio_name] = df_model[line_label_dict[ynum_line]] / df_model[line_label_dict[yden_line]]
        log_model_ratio = np.log10(df_model[ratio_name])
        max_model_ratio = np.max(log_model_ratio)
        min_model_ratio = np.min(log_model_ratio)

        # --------starting the plot-----------
        ratio_name_int = 'log_' + ratio.replace('/', '_') + '_int'
        ratio_name_sum = 'log_' + ratio.replace('/', '_') + '_sum'

        for marker in pd.unique(df['marker']):
            df_sub = df[df['marker'] == marker]
            ax.scatter(df_sub[ratio_name_int], df_sub[ratio_name_sum], c=df_sub['color'], marker=marker, s=50, edgecolor='k', lw=0.5, cmap=args.cmap, vmin=args.clim[0], vmax=args.clim[1])
        ax.errorbar(df[ratio_name_int], df[ratio_name_sum], xerr=df[ratio_name_int + '_u'], yerr=df[ratio_name_sum + '_u'], fmt='none', c='gray', lw=0.5, zorder=-2)

        # --------annotating the figure-----------
        ax.set_xlabel(f'log {ratio} (int)', fontsize=args.fontsize)
        ax.set_ylabel(f'log {ratio} (sum)', fontsize=args.fontsize)

        ax.plot([lim[0], lim[1]], [lim[0], lim[1]], c='k', ls='dotted', lw=1, zorder=-5)

        # -------adding the model limits-----------
        ax.add_patch(plt.Rectangle((min_model_ratio, min_model_ratio), max_model_ratio - min_model_ratio, max_model_ratio - min_model_ratio, lw=0, fc=model_color, alpha=0.2, zorder=-5))
        #ax.scatter(log_model_ratio, log_model_ratio, c=df_model['Z'], s=5, lw=0, marker='d', cmap='cividis', vmin=8, vmax=9, zorder=-4) ##
        ax.text(0.05, 0.05, f'{phot_model_text}\nmodel coverage', c=model_color, fontsize=args.fontsize, va='bottom', ha='left', transform=ax.transAxes)

        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[0], lim[1])

    # --------making colorbar if needed-------------
    if not args.nocolorcoding and args.colorcol != 'ez_z_phot':
        cax = fig.add_axes([0.93, 0.3, 0.005, 0.5])
        cbar = matplotlib.colorbar.ColorbarBase(cax, orientation='vertical', cmap=args.cmap, norm=mplcolors.Normalize(args.clim[0], args.clim[1]), label=args.colorcol)

    # ----------saving the figure-----------------
    colorby_text = '' if args.nocolorcoding or args.colorcol == 'ez_z_phot' else f'_colby_{args.colorcol}'
    figname = args.output_dir / 'plots' / f'sum_vs_int_line_ratio_comparison_{args.snr_text}{args.only_seg_text}_{",".join(ratio_list).replace("/", "-")}{args.colorby_text}.png'
    ax.figure.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def compare_metallicities(df, Zdiag_arr, args, lim=[7, 9]):
    '''
    Plots and saves the integrated vs summed metallicities from a given dataframe
    Returns the figure handle
    '''
    # ------setting up the figure-------------
    fig, axes = plt.subplots(1, len(Zdiag_arr), figsize=(14, 4), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.15, top=0.95, wspace=0.2)

    # -----------plotting the line ratios-------------
    for Zdiag, ax in zip(Zdiag_arr, axes):
        Zdiag_name_int = f'Z_{Zdiag}_int'
        Zdiag_name_sum = f'Z_{Zdiag}_sum'

        for marker in pd.unique(df['marker']):
            df_sub = df[df['marker'] == marker]
            ax.scatter(df_sub[Zdiag_name_int], df_sub[Zdiag_name_sum], c=df_sub['color'], marker=marker, s=50, edgecolor='k', lw=0.5, cmap=args.cmap, vmin=args.clim[0], vmax=args.clim[1])
        ax.errorbar(df[Zdiag_name_int], df[Zdiag_name_sum], xerr=df[Zdiag_name_int + '_u'], yerr=df[Zdiag_name_sum + '_u'], fmt='none', c='gray', lw=0.5, zorder=-2)

        # --------annotating the figure-----------
        ax.set_xlabel(f'log O/H + 12: {Zdiag} (int)', fontsize=args.fontsize)
        ax.set_ylabel(f'log O/H + 12: {Zdiag} (sum)', fontsize=args.fontsize)

        ax.plot([lim[0], lim[1]], [lim[0], lim[1]], c='k', ls='dotted', lw=1)

        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[0], lim[1])

    # --------making colorbar if needed-------------
    if not args.nocolorcoding and args.colorcol != 'ez_z_phot':
        cax = fig.add_axes([0.832, 0.92, 0.15, 0.01])
        cbar = matplotlib.colorbar.ColorbarBase(cax, orientation='horizontal', cmap=args.cmap, norm=mplcolors.Normalize(args.clim[0], args.clim[1]), label=args.colorcol)

    # ----------saving the figure-----------------
    figname = args.output_dir / 'plots' / f'sum_vs_int_metallicity_comparison_{args.snr_text}{args.only_seg_text}_{",".join(Zdiag_arr)}_Zbranch-{args.Zbranch}_AGNdiag_{args.AGN_diag}{args.exclude_text}{args.colorby_text}.png'
    ax.figure.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -----------setting up hard-coded values-------------------
    args.fontsize = 10
    args.only_seg = True
    args.vorbin = True
    args.AGN_diag = 'Ne3O2'
    args.voronoi_line = 'NeIII-3867'
    args.voronoi_snr = float(4)
    args.exclude_lines = 'SII'
    args.version_dict = {'passage': 'v0.5', 'glass': 'orig'}

    args.line_list = 'OII,NeIII-3867,Hb,OIII'.split(',')
    args.Zdiag = 'O3O2,R2,R3,R23,NB'.split(',')
    ratio_list = ['OIII/Hb', 'OII/Hb', 'OIII/OII', 'NeIII-3867/OII']

    col_arr = np.tile(['salmon', 'sienna', 'cornflowerblue', 'darkolivegreen', 'darkgoldenrod', 'darkorchid', 'darkcyan', 'hotpink', 'rosybrown', 'peru', 'darkorange'], 30)
    args.cmap = 'cividis'
    args.clim = [7, 9] if 'Z' in args.colorcol else [None, None]

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
    args.colorby_text = '' if args.nocolorcoding or args.colorcol == 'ez_z_phot' else f'_colby_{args.colorcol}'
    args.exclude_text = f'_without_{args.exclude_lines}' if len(args.exclude_lines) > 0 else ''

    # -------loading the dataframe------------------
    output_dfname = args.output_dir / 'catalogs' / f'sum_vs_int_comparison_{args.snr_text}{args.only_seg_text}{args.vorbin_text}_Zbranch-{args.Zbranch}_AGNdiag_{args.AGN_diag}{args.exclude_text}.csv'

    if not os.path.exists(output_dfname) or args.clobber:
        # ---------setting up the dataframe---------------
        df = pd.DataFrame(columns=np.hstack([['field', 'objid', 'redshift'], \
                                             np.hstack([[f'{item}_int', f'{item}_int_u'] for item in args.line_list]), \
                                             np.hstack([[f'{item}_sum', f'{item}_sum_u'] for item in args.line_list]),
                                             np.hstack([[f'Z_{item}_int', f'Z_{item}_int_u'] for item in args.Zdiag]), \
                                             np.hstack([[f'Z_{item}_sum', f'Z_{item}_sum_u'] for item in args.Zdiag]), \
                                             ]))

        # -----------gathering emission line fluxes and metallicties for all objects into a dataframe-----------
        for index, obj in enumerate(objlist):
            field = obj[0]
            objid = obj[1]
            print(f'\nDoing object {field}:{objid} which is {index + 1} of {len(objlist)}..')

            full_hdu = load_full_fits(objid, field, args)
            args = load_object_specific_args(full_hdu, args)
            z = full_hdu[0].header['redshift']

            _, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, args.line_list, args, silent=True)

            logOH_int_arr, logOH_sum_arr = [], []
            for Zdiag in args.Zdiag:
                _, logOH_int, logOH_sum = load_metallicity_map(field, objid, Zdiag, args)
                logOH_int_arr.append(logOH_int)
                logOH_sum_arr.append(logOH_sum)

            this_row = np.hstack([[field, objid, z], \
                                  np.hstack([[item.n, item.s] for item in line_int_arr]), \
                                  np.hstack([[item.n, item.s] for item in line_sum_arr]), \
                                  np.hstack([[item.n, item.s] for item in logOH_int_arr]), \
                                  np.hstack([[item.n, item.s] for item in logOH_sum_arr]), \
                                  ])
            this_df = pd.DataFrame(dict(map(lambda i, j: (i, [j]), df.columns, this_row)))
            df = pd.concat([df, this_df])

        # --------saving dataframe------------
        df.to_csv(output_dfname, index=None)
        print(f'Saved comparison df as {output_dfname}')

    else:
        print(f'Reading comparison df from existing {output_dfname}')
        df = pd.read_csv(output_dfname)

    # --------------prepare dataframe---------------
    df['marker'] = df['field'].map(lambda x: 's' if 'glass' in x else 'o')
    if args.nocolorcoding: df['color'] = np.tile(['brown'], len(df))
    elif args.colorcol == 'ez_z_phot': df['color'] = col_arr[: len(df)]
    else: df['color'] = df[args.colorcol]

    for line in args.line_list:
        line_name_int = f'{line}_int'
        line_name_sum = f'{line}_sum'

        line_int = unp.uarray(df[line_name_int].astype(float), df[line_name_int + '_u'].astype(float))
        line_sum = unp.uarray(df[line_name_sum].astype(float), df[line_name_sum + '_u'].astype(float))

        line_int[line_int <= 0] = ufloat(np.nan, np.nan)
        line_sum[line_sum <= 0] = ufloat(np.nan, np.nan)

        log_line_int = unp.log10(line_int)
        log_line_sum = unp.log10(line_sum)

        df['log_' + line_name_int] = unp.nominal_values(log_line_int)
        df['log_' + line_name_int + '_u'] = unp.std_devs(log_line_int)

        df['log_' + line_name_sum] = unp.nominal_values(log_line_sum)
        df['log_' + line_name_sum + '_u'] = unp.std_devs(log_line_sum)

    #df = df[df['objid'].isin([2171])]

    # -----------taking the ratios-------------------------
    for ratio in ratio_list:
        num, den = ratio.split('/')

        ratio_name_int = 'log_' + ratio.replace('/', '_') + '_int'
        ratio_int = unp.uarray(df[num+ '_int'].astype(float), df[num + '_int_u'].astype(float)) / unp.uarray(df[den+ '_int'].astype(float), df[den + '_int_u'].astype(float))
        ratio_int[ratio_int <= 0] = ufloat(np.nan, np.nan)
        quant = unp.log10(ratio_int)
        df[ratio_name_int] = unp.nominal_values(quant)
        df[ratio_name_int + '_u'] = unp.std_devs(quant)

        ratio_name_sum = 'log_' + ratio.replace('/', '_') + '_sum'
        ratio_sum = unp.uarray(df[num+ '_sum'].astype(float), df[num + '_sum_u'].astype(float)) / unp.uarray(df[den+ '_sum'].astype(float), df[den + '_sum_u'].astype(float))
        ratio_sum[ratio_sum <= 0] = ufloat(np.nan, np.nan)
        quant = unp.log10(ratio_sum)
        df[ratio_name_sum] = unp.nominal_values(quant)
        df[ratio_name_sum + '_u'] = unp.std_devs(quant)

    # -------making the comparison figurse--------------
    #fig_models_arr = overplot_on_photoionisation_models(df, ratio_list, 'Z', args, col_int='sandybrown', col_sum='cornflowerblue')
    #fig_flux = compare_line_fluxes(df, args.line_list, args, lim=[-19, -15.5])
    #fig_ratio = compare_line_ratios(df, ratio_list, args, lim=[-2, 2])
    fig_Z = compare_metallicities(df, args.Zdiag, args, lim=[7, 9])

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
