'''
    Filename: plots_for_stacked_paper.py
    Notes: Plots stacked mass-metallicity and mass-metallicity gradient relations vs literature plots
    Author : Ayan
    Created: 16-07-26
    Example: run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag R23 --use_C25 --adaptive_bins --bin_by_distance_mass --fold_maps --skip_deproject --cut_z_flag 4
             run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag NB --adaptive_bins --bin_by_distance_mass --fold_maps --skip_deproject --cut_z_flag 4
             run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag NB --adaptive_bins --bin_by_sfh_mass --fold_maps --skip_deproject
'''

from header import *
from util import *
setup_plot_style()
from make_sfms_bins import log_mass_bins, log_sfr_bins, get_stacking_sample, get_binned_df, get_sfms_func, sfms, required_lines, passage_catalog
from make_passage_plots import plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14, plot_SFMS_PASSAGE
from plot_stacked_gradients import read_stacked_df, label_dict, fix_interval_precision

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_MZR_literature(ax):
    '''
    Overplots literature values of MZR on a given axis handle
    Returns axis handle
    '''
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

    zahid_data = zahid_data[zahid_data['Redshift'] > 1].reset_index(drop=True) # removing low-z surveys

    col_ar = ['orange', 'black', 'cyan', 'seagreen', 'blue', 'salmon', 'gray']
    for i in range(len(zahid_data)):
        xarr = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)
        ax.plot(xarr, zahid_func(xarr, zahid_data['Z0'][i], zahid_data['M0'][i], zahid_data['gamma'][i]), color=col_ar[i], lw=2, ls=zahid_data['linestyle'][i], label= 'z = '+str(zahid_data['Redshift'][i]) + '; ' + zahid_data['Sample'][i], zorder=-5)

    # Nedkova+26 MZR
    N26_coeff = unp.uarray([-0.046, 0.992, 3.085], [0.012, 0.204, 0.863]) # from eq 1 of Nedkova+26
    ax.plot(xarr, np.poly1d(unp.nominal_values(N26_coeff))(xarr), color='red', lw=2, ls='--', label='1.7 < z < 3.4; N26')
    #ax.fill_between(xarr, np.poly1d(unp.nominal_values(N26_coeff) - unp.std_devs(N26_coeff))(xarr), np.poly1d(unp.nominal_values(N26_coeff) + unp.std_devs(N26_coeff))(xarr), color='salmon', alpha=0.5)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_MZR(df, args, xcol='log_mass_median', ycol='logOH_int', colorcol=None, cmap='RdBu', qualifiers=''):
    '''
    Plots the stacked mass-metallicity relation, overplotted with relations from the literature
    Saves the figure
    Returns figure handle
    '''
    # ------setup figure----------
    fig, ax = plt.subplots(1, 1, figsize = (6., 5))
    fig.subplots_adjust(left=0.14, right=0.82, top=0.95, bottom=0.12, wspace=0., hspace=0.)
    
    print_sr_corr(df, xcol, ycol, xcol2=colorcol)
    log_mass_cut_low, log_mass_cut_high = 8.5, 8.5
    df_low = df[df[xcol] < log_mass_cut_low]
    df_high = df[df[xcol] > log_mass_cut_high]
    print(f'\nAfter log_mass < {log_mass_cut_low}..')
    print_sr_corr(df_low, xcol, ycol, xcol2=colorcol)
    print(f'\nAfter log_mass > {log_mass_cut_high}..')
    print_sr_corr(df_high, xcol, ycol, xcol2=colorcol)
    
    # ------plot data-----------
    if colorcol is None:
        color = 'cornflowerblue'
        cmin, cmax = None, None
        cmap = None
    else:
        color = df[colorcol]
        cmin, cmax = lim_dict[colorcol][0], lim_dict[colorcol][1]
        cmap = cmap
            
    p = ax.scatter(df[xcol], df[ycol], s=100, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker='o')
    if f'{ycol}_u' in df:
        ax.errorbar(df[xcol], df[ycol], yerr=df[f'{ycol}_u'], c='grey', lw=0.7, fmt='none', alpha=1)

    # -----plot literature---------
    ax = plot_MZR_literature(ax)
    args.fontfactor *= 1.3
    ax.legend(fontsize=args.fontsize / args.fontfactor, loc='best')

    # -------annotate and save fig--------
    ax = annotate_axes(ax, label_dict[xcol], label_dict[ycol], xlim=lim_dict[xcol], ylim=lim_dict[ycol], args=args, 
                       clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, cbar_width=2,
                       p=p, hide_cbar_ticks=False, cticks_integer=False)    
    figname = f'MZR_{qualifiers}.png'
    save_fig(fig, args.fig_dir, figname, args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def plot_MZGR_literature(ax, this_work_legend=[], skip_legend=False):
    '''
    Overplots literature values of MZGR (in units of dex/re only) on a given axis handle
    Returns axis handle
    '''
    # ---------for plotting other observed data from literature------------
    legend_dict = {'sami': 'SAMI', 'manga': 'MaNGA', 'califa': 'CALIFA', 'sharda_scaling1': 'S21 scaling 1', 'sharda_scaling2': 'S21 scaling 2', 'mingozzi2020_izi': 'Mingozzi+20 (IZI)', 'wang17': 'Wang+17', 'jones15': 'Jones+15', 'venturi24': 'Venturi+24', 'li25': 'Li+25', 'ju25': 'Ju+25', 'khoram25': r'Khoram+25 (T$_e$-based)', 'acharyya26': 'Acharyya+26'}
    marker_dict = {'sami': '+', 'manga': 'x', 'califa': '1', 'sharda_scaling1': 'v', 'sharda_scaling2': '^', 'mingozzi2020_izi': '<', 'wang17': 'v', 'jones15': '^', 'venturi24': '>', 'li25': 'D', 'ju25': 'd', 'khoram25': 'o', 'acharyya26': 'X'}
    ls_dict = {'sami': 'dotted', 'manga': 'dotted', 'califa': 'dotted', 'sharda_scaling1': 'solid', 'sharda_scaling2': 'dashed', 'mingozzi2020_izi': 'dotted', 'wang17': 'dotted', 'jones15': 'dotted', 'venturi24': 'dotted', 'li25': 'dotted', 'ju25': 'dotted', 'khoram25': 'dotted', 'acharyya26': 'dotted'}
    color_dict = {'sami': 'firebrick', 'manga': 'chocolate', 'califa': 'darkgoldenrod', 'sharda_scaling1': 'k', 'sharda_scaling2': 'k', 'mingozzi2020_izi': 'peru', 'wang17': 'brown', 'jones15': 'sandybrown', 'venturi24': 'bisque', 'li25': 'grey', 'ju25': 'lightskyblue', 'khoram25': 'cornflowerblue', 'acharyya26': 'goldenrod'}
 
   # --------plotting Sharda+21 data: for dex/re----------
    s21 = []
    literature_dir = args.root_dir / 'zgrad_paper_plots' / 'literature'
    search_text = 'mzgr_*re.csv'
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
    m20 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])

     # --------plotting GLASS-HST data: for dex/re----------
    sample = 'jones15'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv')
    j15 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')
    
    sample = 'wang17'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv')
    w17 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')

    # --------plotting Venturi+24 data: for dex/re and dex/kpc----------
    sample = 'venturi24'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', comment='#', delim_whitespace=True)
    v24 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad_re'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], marker=marker_dict[sample], ec='k')
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad_re'], xerr=df_lit['log_mass_u'], yerr=df_lit['Zgrad_re_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')

    # --------plotting Ju+25 data: for dex/re and dex/kpc----------
    sample = 'ju25'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', comment='#', delim_whitespace=True)
    j25 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad_re'], color=color_dict[sample], lw=0.5, label=legend_dict[sample] + ' (z~1)', marker=marker_dict[sample], ec='k')
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad_re'], yerr=df_lit['Zgrad_re_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')

    # --------plotting direct Te-metallicity gradient using stacked MaNGA data from Khoram+2025: for dex/re----------
    sample = 'khoram25'
    df_lit = pd.read_csv(literature_dir / f'mzgr_{sample}.csv', comment='#')
    k25 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, s=50, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')

    # -------plotting Acharyys+26 data: for dex/re--------------------
    sample = 'acharyya26'
    df_lit = pd.DataFrame({'ID':[300, 1303, 1849, 2867, 1721, 1983, 1991, 1333], \
                            'Redshift':[1.9, 1.9, 3.1, 2.0, 2.2, 1.9, 2.2, 2.0], \
                            'log_mass':[9.5, 8.8, 8.7, 8.7, 8.4, 8.8, 8.0, 8.5], \
                            'log_mass_u':[0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1], \
                            'Zgrad':[-0.53, -0.01, 0.07, 0.45, -0.05, 0.26, 0.08, -0.53], \
                            'Zgrad_u':[0.52, 0.28, 0.14, 0.24, 0.12, 0.36, 0.20, 0.33], \
                            }) # these are taken directly from Table 2 of Acharyya+2026
    a26 = ax.scatter(df_lit['log_mass'], df_lit['Zgrad'], color=color_dict[sample], lw=0.5, label=legend_dict[sample], ec='k', marker=marker_dict[sample])
    ax.errorbar(df_lit['log_mass'], df_lit['Zgrad'], yerr=df_lit['Zgrad_u'], color=color_dict[sample], lw=1, alpha=1, fmt='none')

    # ---------annotate axes and save figure-------
    if not skip_legend:
        handles = this_work_legend + [j15, w17, m20] + s21 + [v24, j25, k25, a26]
        labels = [h.get_label() for h in handles]
        fig = ax.figure
        fig.legend(handles, labels, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.99), fontsize=args.fontsize / args.fontfactor / 1.22)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def print_sr_corr(df, xcol, ycol, xcol2=None):
    '''
    Prints out the partial (if ycol2 is not None) Spearman Rank correlation coefficients
    '''
    p_corr_xcol = pg.corr(df[xcol], df[ycol], method='spearman')
    print(f'\nSpearman Rank correlation of {ycol} vs {xcol} is r={p_corr_xcol.loc["spearman", "r"]:.2f}, p-val={p_corr_xcol.loc["spearman", "p-val"]:.2f}')
    
    if xcol2 is not None:
        p_corr_ccol = pg.corr(df[xcol2], df[ycol], method='spearman')
        print(f'Spearman Rank correlation of {ycol} vs {xcol2} is r={p_corr_ccol.loc["spearman", "r"]:.2f}, p-val={p_corr_ccol.loc["spearman", "p-val"]:.2f}')

        p_corr_xcol = pg.partial_corr(data=df, x=xcol, y=ycol, covar=xcol2, method='spearman')
        p_corr_ccol = pg.partial_corr(data=df, x=xcol2, y=ycol, covar=xcol, method='spearman')
        print(f'Parial Spearman Rank correlation of {ycol} vs {xcol} (keeping {xcol2} fixed) is r={p_corr_xcol.loc["spearman", "r"]:.2f}, p-val={p_corr_xcol.loc["spearman", "p-val"]:.2f}')
        print(f'Parial Spearman Rank correlation of {ycol} vs {xcol2} (keeping {xcol} fixed) is r={p_corr_ccol.loc["spearman", "r"]:.2f}, p-val={p_corr_ccol.loc["spearman", "p-val"]:.2f}')

    return

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_MZGR(df, args, xcol='log_mass_median', ycol='radial_logOH_grad', colorcol=None, cmap='RdBu', qualifiers=''):
    '''
    Plots the stacked mass-metallicity gradient relation (both radial and minor-major gradients), overplotted with relations from the literature
    Saves the figure
    Returns figure handle
    '''
    # ------setup figure----------
    fig, axes = plt.subplots(2, 1, figsize = (10, 7.6), sharex=True)
    fig.subplots_adjust(left=0.1, right=0.87, top=0.88, bottom=0.08, wspace=0., hspace=0.04)

    # ------prepare plotting attributes-----------
    df = df.sort_values(by=xcol)
    print_sr_corr(df, xcol, ycol, xcol2=colorcol)

    log_mass_cut_low, log_mass_cut_high = 8.5, 8.5
    df_low = df[df[xcol] < log_mass_cut_low]
    df_high = df[df[xcol] > log_mass_cut_high]
    print(f'\nAfter log_mass < {log_mass_cut_low}..')
    print_sr_corr(df_low, xcol, ycol, xcol2=colorcol)
    print(f'\nAfter log_mass > {log_mass_cut_high}..')
    print_sr_corr(df_high, xcol, ycol, xcol2=colorcol)

    if colorcol is None:
        color = 'cornflowerblue'
        cmin, cmax = None, None
        cmap = None
    else:
        color = df[colorcol]
        cmin, cmax = lim_dict[colorcol][0], lim_dict[colorcol][1]
        cmap = cmap

    # -------plot radial gradient------------
    axes[0].plot(df[xcol], df[ycol], lw=0.7, c='k', ls='dashed')
    p = axes[0].scatter(df[xcol], df[ycol], s=100, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker='o', zorder=20, label=f'This Work (azimuthally averaged)')
    if f'{ycol}_u' in df:
        axes[0].errorbar(df[xcol], df[ycol], yerr=df[f'{ycol}_u'], c='grey', lw=0.7, fmt='none', alpha=1)
    this_work_legend = [p]

    axes[0].axhline(0, ls='dashed', lw=0.5, c='k')

    # -----plot literature---------
    axes[0] = plot_MZGR_literature(axes[0], skip_legend=True)

    # -------annotate axis--------
    axes[0] = annotate_axes(axes[0], label_dict[xcol], label_dict[ycol], xlim=lim_dict[xcol], ylim=lim_dict[ycol], args=args, 
                       clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, cbar_width=2,
                       p=p, hide_cbar_ticks=False, cticks_integer=False, hide_xaxis=True)    

    vline_col = 'sienna'
    log_mass_cut = log_mass_cut_low
    axes[0].axvline(log_mass_cut, ls='dotted', lw=1., c=vline_col)
    axes[0].text(log_mass_cut - 0.1, axes[0].get_ylim()[1] * 0.95, 'Lower mass regime', c=vline_col, fontsize=args.fontsize / args.fontfactor, ha='right', va='top')
    axes[0].text(log_mass_cut + 0.1, axes[0].get_ylim()[1] * 0.95, 'Higher mass regime', c=vline_col, fontsize=args.fontsize / args.fontfactor, ha='left', va='top')
    axes[1].axvline(log_mass_cut, ls='dotted', lw=1., c=vline_col)

    # -------plot minor major gradient------------
    ycol_arr = [ycol.replace('radial', 'minor'), ycol.replace('radial', 'major')]
    marker_arr = ['s', 'D']
    ls_arr = ['solid', 'dashed']

    for index, ycol in enumerate(ycol_arr):
        print_sr_corr(df, xcol, ycol, xcol2=colorcol)
        print(f'\nAfter log_mass < {log_mass_cut_low}..')
        print_sr_corr(df_low, xcol, ycol, xcol2=colorcol)
        print(f'\nAfter log_mass > {log_mass_cut_high}..')
        print_sr_corr(df_high, xcol, ycol, xcol2=colorcol)
        
        axes[1].plot(df[xcol], df[ycol], lw=0.7, c='k', ls=ls_arr[index])
        p = axes[1].scatter(df[xcol], df[ycol], s=70, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker=marker_arr[index], zorder=20, label=f'This work ({ycol.split("_")[0]} axis)')
        if f'{ycol}_u' in df:
            axes[1].errorbar(df[xcol], df[ycol], yerr=df[f'{ycol}_u'], c='grey', lw=0.7, fmt='none', alpha=1)
        this_work_legend.append(p)

    axes[1].axhline(0, ls='dashed', lw=0.5, c='k')

    # -----plot literature---------
    axes[1] = plot_MZGR_literature(axes[1], this_work_legend=this_work_legend)

    # -------annotate axis--------
    ycol = ycol.replace('minor', 'radial').replace('major', 'radial')
    axes[1] = annotate_axes(axes[1], label_dict[xcol], label_dict[ycol], xlim=lim_dict[xcol], ylim=lim_dict[ycol], args=args, 
                       clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, cbar_width=2,
                       p=p, hide_cbar_ticks=False, cticks_integer=False)    

    # -------save fig--------
    figname = f'MZGR_{qualifiers}.png'
    save_fig(fig, args.fig_dir, figname, args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
lim_dict = {'minor_logOH_grad': [-1.2, 1.2],\
                'major_logOH_grad': [-1.2, 1.2],\
                'radial_logOH_grad': [-1.2, 1.2],\
                'logOH_int': [7.0, 9.0],\
                'delta_sfms_median': [-0.6, 0.6],\
                'log_mass_median': [7.0, 10.0],\
                'tform_ratio_median': [0, 1],\
                }  
# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.2
    
    # ---------reading in the master SED catalog----------------
    passage_catalog_filename = args.output_dir / 'catalogs' / passage_catalog
    df_input = get_stacking_sample(passage_catalog_filename, args, required_lines=required_lines, sfms=sfms)

    # ------------reading and binning dataframe-------------
    interval_cols = ['delta_sfms_bin', 'log_mass_bin', 'log_sfr_bin', 'mass_interval', 'mass_intervals', 'bin_intervals', 'sfr_interval', 'sfr_intervals']
    df = df_input.copy()
    df, bin_list, args = get_binned_df(args, df=df, skip_stacking=True, required_lines=required_lines, sfms=sfms)
    for col in interval_cols:
        if col in df: 
            try: df[col] = fix_interval_precision(df[col], precision=3)
            except: pass

    # -------------reading in stacked gradient dataframe-----------------------
    df_grad = read_stacked_df(args.grad_filename)
    df_grad['nobj'] = df_grad['nobj'].astype(int)
    for col in interval_cols:
        if col in df_grad: df_grad[col] = fix_interval_precision(df_grad[col], precision=3)

    # -------------merging the two dataframes-----------------------
    if args.bin_by_distance_mass:
        df2 = df.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'bin_intervals':'delta_sfms_bin'})
        df2 = df2.groupby(['delta_sfms_bin', 'log_mass_bin']).agg(log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
        df_grad = pd.merge(df_grad, df2, on=['delta_sfms_bin', 'log_mass_bin'], how='left').reset_index(drop=True)
    elif args.bin_by_distance:
        df2 = df.rename(columns={'bin_intervals':'delta_sfms_bin'})
        df2 = df2.groupby(['delta_sfms_bin']).agg(log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
        df_grad = pd.merge(df_grad, df2, on=['delta_sfms_bin'], how='left').reset_index(drop=True)
    elif args.bin_by_sfh_mass:
        df2 = df.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'bin_intervals':'tform_ratio_bin'})
        df2 = df2.groupby(['tform_ratio_bin', 'log_mass_bin']).agg(log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median'), tform_ratio_median=('delta_tform_ratio', 'median')).reset_index()
        df_grad = pd.merge(df_grad, df2, on=['tform_ratio_bin', 'log_mass_bin'], how='left').reset_index(drop=True)
    else:
        df2 = df.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'sfr_intervals':'log_sfr_bin', 'sfr_interval':'log_sfr_bin'})
        df2 = df2.groupby(['log_mass_bin', 'log_sfr_bin']).agg(log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
        df_grad = pd.merge(df_grad, df2, on=['log_mass_bin', 'log_sfr_bin'], how='left').reset_index(drop=True)

    # ------------plotting stacked MZR and MZGR--------------------------
    qualifiers = f'{args.binby_text}{args.fold_text}_Zdiag_{args.Zdiag}{args.C25_text}{args.deproject_text}{args.rescale_text}'
    if args.bin_by_distance_mass:
        colorcol, cmap = 'delta_sfms_median', 'PRGn' # diverging cmap
    else:
        colorcol, cmap = 'tform_ratio_median', 'viridis' # sequential cmap
    
    fig_mzr = plot_stacked_MZR(df_grad, args, xcol='log_mass_median', ycol='logOH_int', colorcol=colorcol, qualifiers=qualifiers, cmap=cmap)
    #fig_mzgr = plot_stacked_MZGR(df_grad, args, xcol='log_mass_median', ycol='radial_logOH_grad', colorcol=colorcol, qualifiers=qualifiers, cmap=cmap)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
