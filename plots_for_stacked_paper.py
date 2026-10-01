'''
    Filename: plots_for_stacked_paper.py
    Notes: Plots stacked mass-metallicity and mass-metallicity gradient relations vs literature plots
    Author : Ayan
    Created: 16-07-26
    Example: run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag R23 --use_C25 --adaptive_bins --bin_by_distance_mass --fold_maps --skip_deproject --cut_z_flag 4
             run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag NB --adaptive_bins --bin_by_distance_mass --fold_maps --skip_deproject --cut_z_flag 4
             run plots_for_stacked_paper.py --system ssd --do_all_fields --Zdiag NB --adaptive_bins --bin_by_sfh_mass --fold_maps --skip_deproject --cut_z_flag 4
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

    # -------------Nedkova+26 MZR--------------------------
    #N26_coeff = unp.uarray([-0.046, 0.992, 3.085], [0.012, 0.204, 0.863]) # from eq 1 of Nedkova+26
    N26_coeff = unp.uarray([-0.038, 0.856, 3.669], [0.021, 0.372, 1.613]) # from Kalina's slack msg
    logOH_mid = np.poly1d(unp.nominal_values(N26_coeff))(xarr)
    ax.plot(xarr, logOH_mid, color='red', lw=2, ls='--', label='1.7 < z < 3.4; N26')
    '''
    # -------this is my attempt to propagate errors---------------
    logOH_low = np.poly1d(unp.nominal_values(N26_coeff) - unp.std_devs(N26_coeff))(xarr)
    logOH_high = np.poly1d(unp.nominal_values(N26_coeff) + unp.std_devs(N26_coeff))(xarr)
    '''
    # --------the following is code from Kalina--------------
    # Covariance matrix extracted from the MCMC chains (accounting for correlation); Order of parameters: [a, b, c]
    cov_matrix = np.array([
        [ 0.00046, -0.00790,  0.03360],
        [-0.00790,  0.13840, -0.59600],
        [ 0.03360, -0.59600,  2.60100]
    ])

    # Correctly propagate errors to find the 1-sigma range *around* the curve; Using: sigma_y^2 = J^T * Cov * J  where J is the Jacobian matrix [x^2, x, 1]
    y_fit_err = np.zeros_like(xarr)
    for i, x in enumerate(xarr):
        jacobian = np.array([x**2, x, 1.0])
        # Compute local model variance at this specific x coordinate
        variance = np.dot(jacobian, np.dot(cov_matrix, jacobian))
        y_fit_err[i] = np.sqrt(variance)

    # Define the true 1-sigma upper and lower bounds flanking the best fit
    logOH_low = logOH_mid - y_fit_err
    logOH_high = logOH_mid + y_fit_err
    
    # -------plotting shaded region--------
    ax.fill_between(xarr, logOH_low, logOH_high, color='salmon', alpha=0.5)

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
    Returns a dictionary of the results
    '''
    corr_xcol = pg.corr(df[xcol], df[ycol], method='spearman')
    corr_xcol_r, corr_xcol_p = corr_xcol.loc["spearman", "r"], corr_xcol.loc["spearman", "p-val"]
    print(f'\nSpearman Rank correlation of {ycol} vs {xcol} is r={corr_xcol_r:.2f}, p-val={corr_xcol_p:.2f}')
    
    if xcol2 is not None:
        corr_ccol = pg.corr(df[xcol2], df[ycol], method='spearman')
        corr_ccol_r, corr_ccol_p = corr_ccol.loc["spearman", "r"], corr_ccol.loc["spearman", "p-val"]
        print(f'Spearman Rank correlation of {ycol} vs {xcol2} is r={corr_ccol_r:.2f}, p-val={corr_ccol_p:.2f}')

        pcorr_xcol = pg.partial_corr(data=df, x=xcol, y=ycol, covar=xcol2, method='spearman')
        pcorr_xcol_r, pcorr_xcol_p = pcorr_xcol.loc["spearman", "r"], pcorr_xcol.loc["spearman", "p-val"]
        print(f'Partial Spearman Rank correlation of {ycol} vs {xcol} (keeping {xcol2} fixed) is r={pcorr_xcol_r:.2f}, p-val={pcorr_xcol_p:.2f}')

        pcorr_ccol = pg.partial_corr(data=df, x=xcol2, y=ycol, covar=xcol, method='spearman')
        pcorr_ccol_r, pcorr_ccol_p = pcorr_ccol.loc["spearman", "r"], pcorr_ccol.loc["spearman", "p-val"]
        print(f'Partial Spearman Rank correlation of {ycol} vs {xcol2} (keeping {xcol} fixed) is r={pcorr_ccol_r:.2f}, p-val={pcorr_ccol_p:.2f}')
    else:
        corr_ccol_r, corr_ccol_p, pcorr_xcol_r, pcorr_xcol_p, pcorr_ccol_r, pcorr_ccol_p = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    results = {
        'corr_x_r': corr_xcol_r, 'corr_x_p': corr_xcol_p,
        'corr_c_r': corr_ccol_r, 'corr_c_p': corr_ccol_p,
        'pcorr_x_r': pcorr_xcol_r, 'pcorr_x_p': pcorr_xcol_p,
        'pcorr_c_r': pcorr_ccol_r, 'pcorr_c_p': pcorr_ccol_p,
    }

    return results

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_MZGR_old(df, args, xcol='log_mass_median', ycol='radial_logOH_grad', colorcol=None, cmap='RdBu', qualifiers=''):
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

    log_mass_cut_low, log_mass_cut_high = log_mass_cut, log_mass_cut
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
def format_latex_cell(r, p):
    '''
    Formats correlation coefficient r and p-value into a LaTeX shortstack string with significance bolding.
    '''
    if np.isnan(r) or np.isnan(p):
        return r"\shortstack{N/A}"
    
    direction = "positive" if r > 0 else "negative"
    abs_r = abs(r)
    
    # Qualitative classification
    if np.round(p, 2) > 0.10:
        qual = "No sig. correlation"
    elif 0.05 < np.round(p, 2) <= 0.10:
        qual = f"Marginal {direction}"
    else:
        if abs_r >= 0.65:
            qual = f"Strong {direction}"
        elif abs_r >= 0.30:
            qual = f"Moderate {direction}"
        else:
            qual = f"Weak {direction}"
            
    p_str = "p < 0.001" if p < 0.001 else (f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.3f}")
    
    # Bold significant results (p <= 0.05)
    if np.round(p, 2) <= 0.05:
        return "\\shortstack{\\textbf{" + qual + "}\\\\ (\\textbf{$r = " + f"{r:.2f}, {p_str}" + "$})}"
    else:
        return "\\shortstack{" + qual + "\\\\ ($r = " + f"{r:.2f}, {p_str}" + "$)}"

# --------------------------------------------------------------------------------------------------------------------
def generate_latex_table(results_dict, outfilename):
    '''
    Constructs LaTeX tabular code from collected correlation statistics and prints/saves it
    '''
    
    header = (
        "\\begin{tabular}{l|c|c||c|c}\n"
        "\\toprule\n"
        "& \\multicolumn{2}{c||}{Stellar mass ($\\log{(M_*/M_\\odot)}$)} & \\multicolumn{2}{c}{Offset from SFMS ($\\delta_{\\rm SFMS}$)} \\\\\n"
        "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n"
        "Quantity & Zero-Order & Partial\\tnote{a} & Zero-Order & Partial\\tnote{a} \\\\\n"
        "\\midrule\n"
    )
    
    body = ""
    sections = [
        ('full', r"\multicolumn{5}{l}{\textit{Full mass range ($7 \lesssim \log{(M_*/M_\odot)} \lesssim 11$)}} \\[5pt]"),
        ('low', r"\multicolumn{5}{l}{\textit{Lower-mass regime ($\log{(M_*/M_\odot)} < 8.5$)}} \\[5pt]"),
        ('high', r"\multicolumn{5}{l}{\textit{Higher-mass regime ($\log{(M_*/M_\odot)} > 8.5$)}} \\[5pt]")
    ]
    
    labels = {
        'logOH_int': 'Integrated metallicity',
        'radial_logOH_grad': 'Radial gradient',
        'minor_logOH_grad': 'Minor-axis gradient',
        'major_logOH_grad': 'Major-axis gradient'
    }
    
    for sec_key, sec_title in sections:
        if sec_key not in results_dict or not results_dict[sec_key]:
            continue
            
        body += f"{sec_title}\n"
        rows = results_dict[sec_key]
        
        if sec_key == 'full': quants_in_this_section = ['logOH_int']
        else: quants_in_this_section = ['radial_logOH_grad', 'minor_logOH_grad', 'major_logOH_grad']

        for q_key in quants_in_this_section:
            if q_key not in rows:
                continue
                
            q_label = labels[q_key]
            stats = rows[q_key]
            
            c1 = format_latex_cell(stats['x_zero'][0], stats['x_zero'][1])
            c2 = format_latex_cell(stats['x_partial'][0], stats['x_partial'][1])
            c3 = format_latex_cell(stats['c_zero'][0], stats['c_zero'][1])
            c4 = format_latex_cell(stats['c_partial'][0], stats['c_partial'][1])
            
            spacing = "\\\\[12pt]\n" if q_key != list(rows.keys())[-1] else "\\\\[8pt]\n"
            body += f"{q_label}\n  & {c1}\n  & {c2}\n  & {c3}\n  & {c4} {spacing}"
            
        body += "\\midrule\n" if sec_key != 'high' else ""

    footer = "\\bottomrule\n\\end{tabular}\n"
    latex_str = header + body + footer
        
    with open(outfilename, 'w') as f:
        f.write(latex_str)
    print(f'LaTeX table saved as {outfilename}\n')
        
    return latex_str

# ----------------------------------------------------------------------------------------------------------------
def record_and_print_stats(results_dict, df, xcol, ycol, colorcol=None, log_xcol_cut=8.5):
    '''
    Computes SR correlations, prints them, and adds the result to a given dictionary for later LaTeX table generation
    Returns the dictionary
    '''
    print(f"\nDoing correlations for {ycol}...")

    df_low = df[df[xcol] < log_xcol_cut]
    df_high = df[df[xcol] > log_xcol_cut]

    # Full mass range
    print(f'\nFor full range of {xcol}..')
    results = print_sr_corr(df, xcol, ycol, xcol2=colorcol)
    results_dict['full'][ycol] = {'x_zero': [results['corr_x_r'], results['corr_x_p']],
                                    'c_zero': [results['corr_c_r'], results['corr_c_p']],
                                    'x_partial': [results['pcorr_x_r'], results['pcorr_x_p']],
                                    'c_partial': [results['pcorr_c_r'], results['pcorr_c_p']]
                                    }
    
    # Lower mass range
    print(f'\nAfter slicing to {xcol} < {log_xcol_cut}..')
    results = print_sr_corr(df_low, xcol, ycol, xcol2=colorcol)
    results_dict['low'][ycol] = {'x_zero': [results['corr_x_r'], results['corr_x_p']], 
                                    'c_zero': [results['corr_c_r'], results['corr_c_p']],
                                    'x_partial': [results['pcorr_x_r'], results['pcorr_x_p']],
                                    'c_partial': [results['pcorr_c_r'], results['pcorr_c_p']]
                                    }
    
    # Higher mass range
    print(f'\nAfter slicing to {xcol} > {log_xcol_cut}..')
    results = print_sr_corr(df_high, xcol, ycol, xcol2=colorcol)
    results_dict['high'][ycol] = {'x_zero': [results['corr_x_r'], results['corr_x_p']],
                                    'c_zero': [results['corr_c_r'], results['corr_c_p']],
                                    'x_partial': [results['pcorr_x_r'], results['pcorr_x_p']],
                                    'c_partial': [results['pcorr_c_r'], results['pcorr_c_p']]
                                    }

    return results_dict

# --------------------------------------------------------------------------------------------------------------------
def compute_correlations(df, args, xcol='log_mass_median', ycol='radial_logOH_grad', qualifiers=''):
    '''
    Computes Spearman correlations and automatically outputs a LaTeX table.
    '''
    # ------storing correlations-----------
    df = df.sort_values(by=xcol)
    results_dict = {'full': {}, 'low': {}, 'high': {}}

    # -------for integrated metallicity------------
    results_dict = record_and_print_stats(results_dict, df, xcol, 'logOH_int', colorcol=colorcol, log_xcol_cut=log_mass_cut) # just getting in the stats for integrated metallicity too

    # -------for radial gradient------------
    results_dict = record_and_print_stats(results_dict, df, xcol, ycol, colorcol=colorcol, log_xcol_cut=log_mass_cut)

    # -------for minor & major gradient------------
    results_dict = record_and_print_stats(results_dict, df, xcol, ycol.replace('radial', 'minor'), colorcol=colorcol, log_xcol_cut=log_mass_cut)
    results_dict = record_and_print_stats(results_dict, df, xcol, ycol.replace('radial', 'major'), colorcol=colorcol, log_xcol_cut=log_mass_cut)

    # -------generate LaTeX Table----------
    tex_filename = os.path.join(args.fig_dir, f'MZGR_{qualifiers}_stats_table.tex')
    generate_latex_table(results_dict, tex_filename)

    return results_dict

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_MZGR(df, args, xcol='log_mass_median', ycol='radial_logOH_grad', colorcol=None, cmap='RdBu', qualifiers=''):
    '''
    Plots the stacked mass-metallicity gradient relation (both radial and minor-major gradients), overplotted with relations from the literature.
    Returns figure handle
    '''
    # ------setup figure----------
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.6), sharex=True)
    fig.subplots_adjust(left=0.1, right=0.87, top=0.88, bottom=0.08, wspace=0., hspace=0.04)

    df = df.sort_values(by=xcol)

    # ------prepare plotting attributes-----------
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
    p = axes[0].scatter(df[xcol], df[ycol], s=100, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker='o', zorder=20, label='This Work (azimuthally averaged)')
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
    axes[0].axvline(log_mass_cut, ls='dotted', lw=1., c=vline_col)
    axes[0].text(log_mass_cut - 0.1, axes[0].get_ylim()[1] * 0.95, 'Lower mass regime', c=vline_col, fontsize=args.fontsize / args.fontfactor, ha='right', va='top')
    axes[0].text(log_mass_cut + 0.1, axes[0].get_ylim()[1] * 0.95, 'Higher mass regime', c=vline_col, fontsize=args.fontsize / args.fontfactor, ha='left', va='top')
    axes[1].axvline(log_mass_cut, ls='dotted', lw=1., c=vline_col)

    # -------plot minor & major gradient------------
    ycol_arr = [ycol.replace('radial', 'minor'), ycol.replace('radial', 'major')]
    marker_arr = ['s', 'D']
    ls_arr = ['solid', 'dashed']

    for index, curr_ycol in enumerate(ycol_arr):        
        axes[1].plot(df[xcol], df[curr_ycol], lw=0.7, c='k', ls=ls_arr[index])
        p = axes[1].scatter(df[xcol], df[curr_ycol], s=70, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker=marker_arr[index], zorder=20, label=f'This work ({curr_ycol.split("_")[0]} axis)')
        if f'{curr_ycol}_u' in df:
            axes[1].errorbar(df[xcol], df[curr_ycol], yerr=df[f'{curr_ycol}_u'], c='grey', lw=0.7, fmt='none', alpha=1)
        this_work_legend.append(p)

    axes[1].axhline(0, ls='dashed', lw=0.5, c='k')

    # -----plot literature---------
    axes[1] = plot_MZGR_literature(axes[1], this_work_legend=this_work_legend)

    # -------annotate axis--------
    ycol_rad = ycol.replace('minor', 'radial').replace('major', 'radial')
    axes[1] = annotate_axes(axes[1], label_dict[xcol], label_dict[ycol_rad], xlim=lim_dict[xcol], ylim=lim_dict[ycol_rad], args=args, 
                       clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, cbar_width=2,
                       p=p, hide_cbar_ticks=False, cticks_integer=False)    

    # -------save fig--------
    figname = f'MZGR_{qualifiers}.png'
    save_fig(fig, args.fig_dir, figname, args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def plot_correlation_summary_hatched_old(results_dict, args, qualifiers=''):
    '''
    Plots Spearman rank correlation coefficients (r) across integrated metallicity
    and metallicity gradient quantities in a single-panel grouped bar chart.

    - Bar Fill Color: Red (p <= 0.05) vs. Grey (p > 0.05)
    - Bar Hatching: Solid (Full Mass), '//' (Low Mass), '\\\\' (High Mass)
    - Annotations: Displays r and p values above/below each bar.
    '''
    # Define quantities and correlation types for x-axis items
    # Format: (results_dict ykey, correlation_type, x_tick_label)
    items = [
        ('logOH_int', 'x_zero', r'Integrated $Z$' + '\n' + r'vs $\log{M_*}$'),
        ('logOH_int', 'c_zero', r'Integrated $Z$' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('radial_logOH_grad', 'x_zero', r'Radial Grad' + '\n' + r'vs $\log{M_*}$'),
        ('radial_logOH_grad', 'c_zero', r'Radial Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('minor_logOH_grad', 'x_zero', r'Minor Grad' + '\n' + r'vs $\log{M_*}$'),
        ('minor_logOH_grad', 'c_zero', r'Minor Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('major_logOH_grad', 'x_zero', r'Major Grad' + '\n' + r'vs $\log{M_*}$'),
        ('major_logOH_grad', 'c_zero', r'Major Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
    ]

    regimes = ['full', 'low', 'high']
    regime_labels = ['Full Mass Range', rf'Low Mass ($\log M_* <{log_mass_cut}$)', rf'High Mass ($\log M_* > {log_mass_cut}$)']
    regime_hatches = ['', '//', '\\\\']

    # -------color palette for statistical significance---------
    sig_color = 'brown'    # Coral Red for p <= 0.05
    nonsig_color = 'lightgrey' # Light Grey for p > 0.05

    # ------setting up figure--------------
    fig, ax = plt.subplots(figsize=(13, 5.5))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.98, bottom=0.10)

    x = np.arange(len(items))
    width = 0.25  # Bar width

    # -----------looping through each mass regime to plot grouped bars-----------
    for i, (reg, r_label, hatch) in enumerate(zip(regimes, regime_labels, regime_hatches)):
        r_vals, p_vals = [], []
        
        for q_key, corr_type, _ in items:
            stat = results_dict.get(reg, {}).get(q_key, {}).get(corr_type, [np.nan, np.nan])
            r_vals.append(stat[0])
            p_vals.append(stat[1])
        
        pos = x + (i - 1) * width
        
        # ------color per bar depending on p-value--------
        colors = [sig_color if (p is not np.nan and not np.isnan(p) and p <= 0.05) else nonsig_color for p in p_vals]
        bars = ax.bar(pos, r_vals, width, color=colors, hatch=hatch, edgecolor='black', linewidth=0.8, alpha=0.9, zorder=3)
        
        # --------annotating r and p values above/below bars-------
        for bar, r, p in zip(bars, r_vals, p_vals):
            #if np.isnan(r) or np.isnan(p):
            if np.isnan(r) or np.isnan(p) or p > 0.05: # only print p-values when significant
                continue
            height = bar.get_height()
            va = 'bottom' if height >= 0 else 'top'
            offset = 0.03 if height >= 0 else -0.03
            
            p_str = f"p<{0.01:.2f}" if p < 0.01 else f"p={p:.2f}"
            fontweight = 'bold' if p <= 0.05 else 'normal'
            
            ax.text(bar.get_x() + bar.get_width() / 2.0, height + offset, f"r={r:.2f}\n({p_str})", ha='center', va=va, fontsize=7.5, fontweight=fontweight, color='black', zorder=4)

    # --------annotate plot------------------
    ax.axhline(0, color='black', linewidth=1.0, zorder=2)

    ax.set_xticks(x)
    ax.set_xticklabels([item[2] for item in items], fontsize=args.fontsize / args.fontfactor)
    ax.set_ylabel(r'Spearman Rank Correlation Coefficient ($r$)', fontsize=args.fontsize / args.fontfactor)
    ax.set_ylim(-1.2, 1.2)
    ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)

    # --------making significance legend (left)---------
    legend_sig = [patches.Patch(facecolor=sig_color, edgecolor='black', label=r'$p \leq 0.05$ (Significant)'), patches.Patch(facecolor=nonsig_color, edgecolor='black', label=r'$p > 0.05$ (Not Significant)')]
    leg1 = ax.legend(handles=legend_sig, loc='upper left', frameon=True, fontsize=args.fontsize / args.fontfactor)
    ax.add_artist(leg1)

    # -------making mass legend (right)-------------------
    legend_regime = [patches.Patch(facecolor='white', edgecolor='black', hatch=h, label=lbl)for h, lbl in zip(regime_hatches, regime_labels)]
    ax.legend(handles=legend_regime, loc='upper right', frameon=True, fontsize=args.fontsize / args.fontfactor)

    # -------save fig--------
    figname = f'correlation_summary_{qualifiers}.png'
    save_fig(fig, args.fig_dir, figname, args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def get_p_color_binned(p):
    '''
    Determine the color of bars for the correlation plot; Idea is to have reddish colors for significant,
    dark grey for kinda significant and light gray for non-significant correlations
    '''
    if np.isnan(p):
        return p_value_color_config[-1]['color']
    
    for bin_info in p_value_color_config:
        if p <= bin_info['max_p']:
            return bin_info['color']
            
    return p_value_color_config[-1]['color']
  
# --------------------------------------------------------------------------------------------------------------------
def plot_correlation_summary_hatched(results_dict, args, qualifiers=''):
    '''
    Plots Spearman rank correlation coefficients (r) across integrated metallicity
    and metallicity gradient quantities in a 2-panel grouped bar chart:
    - Top Panel: Zeroth-order correlations
    - Bottom Panel: Partial correlations (controlling for the second variable)

    - Bar Fill Color: Brown (p <= 0.05) vs. Light Grey (p > 0.05)
    - Bar Hatching: Solid (Full Mass), '//' (Low Mass), '\\\\' (High Mass)
    - Annotations: Displays r and p values above/below significant bars.
    '''
    # Define quantities and correlation types for zeroth-order (top) and partial (bottom)
    items_zero = [
        ('logOH_int', 'x_zero', r'Integrated $Z$' + '\n' + r'vs $\log{M_*}$'),
        ('logOH_int', 'c_zero', r'Integrated $Z$' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('radial_logOH_grad', 'x_zero', r'Radial Grad' + '\n' + r'vs $\log{M_*}$'),
        ('radial_logOH_grad', 'c_zero', r'Radial Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('minor_logOH_grad', 'x_zero', r'Minor Grad' + '\n' + r'vs $\log{M_*}$'),
        ('minor_logOH_grad', 'c_zero', r'Minor Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('major_logOH_grad', 'x_zero', r'Major Grad' + '\n' + r'vs $\log{M_*}$'),
        ('major_logOH_grad', 'c_zero', r'Major Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
    ]

    items_partial = [
        ('logOH_int', 'x_partial', r'Integrated $Z$' + '\n' + r'vs $\log{M_*}$'),
        ('logOH_int', 'c_partial', r'Integrated $Z$' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('radial_logOH_grad', 'x_partial', r'Radial Grad' + '\n' + r'vs $\log{M_*}$'),
        ('radial_logOH_grad', 'c_partial', r'Radial Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('minor_logOH_grad', 'x_partial', r'Minor Grad' + '\n' + r'vs $\log{M_*}$'),
        ('minor_logOH_grad', 'c_partial', r'Minor Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
        ('major_logOH_grad', 'x_partial', r'Major Grad' + '\n' + r'vs $\log{M_*}$'),
        ('major_logOH_grad', 'c_partial', r'Major Grad' + '\n' + r'vs $\delta_{\rm SFMS}$'),
    ]

    regimes = ['full', 'low', 'high']
    regime_labels = ['Full Mass Range', rf'Low Mass ($\log M_*/M_\odot <{log_mass_cut}$)', rf'High Mass ($\log M_*/M_\odot > {log_mass_cut}$)']
    regime_hatches = ['', '//', '\\\\']

    # -------color palette for statistical significance---------
    sig_color = 'brown'    # for p <= 0.05
    nonsig_color = 'lightgrey' # for p > 0.05

    # ------setting up figure with shared x-axis and small panel spacing--------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.99, bottom=0.07, hspace=0.02)

    x = np.arange(len(items_zero))
    width = 0.25  # bar width

    panels = [
        (ax1, items_zero, r'Zeroth-Order Spearman ($r$)', 'Zeroth-Order Correlations'),
        (ax2, items_partial, r'Partial Spearman ($r_{\partial}$)', r'Partial Correlations (controlling for $\delta_{\rm SFMS}$ or $\log M_*$)')
    ]

    for ax, items, ylabel_str, panel_title in panels:
        # -----------looping through each mass regime to plot grouped bars-----------
        for i, (reg, r_label, hatch) in enumerate(zip(regimes, regime_labels, regime_hatches)):
            r_vals, p_vals = [], []
            
            for q_key, corr_type, _ in items:
                stat = results_dict.get(reg, {}).get(q_key, {}).get(corr_type, [np.nan, np.nan])
                r_vals.append(stat[0])
                p_vals.append(stat[1])
            
            pos = x + (i - 1) * width
            
            # ------color per bar depending on p-value--------
            if args.nocolorcoding:
                colors = [sig_color if (p is not np.nan and not np.isnan(p) and p <= 0.05) else nonsig_color for p in p_vals] # fixed colors
            else:
                colors = [get_p_color_binned(p) for p in p_vals] # discrete color-map
            bars = ax.bar(pos, r_vals, width, color=colors, hatch=hatch, edgecolor='black', linewidth=0.8, alpha=0.9, zorder=3)
            
            # --------annotating r and p values above/below bars-------
            for bar, r, p in zip(bars, r_vals, p_vals):
                if np.isnan(r) or np.isnan(p) or (p > 0.05 and args.nocolorcoding) or p > 0.1: # only print p-values when significant
                    continue
                height = bar.get_height()
                va = 'bottom' if height >= 0 else 'top'
                offset = 0.03 if height >= 0 else -0.03
                
                p_str = f"p<{0.01:.2f}" if p < 0.01 else f"p={p:.2f}"
                fontweight = 'bold' if p <= 0.05 else 'normal'
                
                ax.text(bar.get_x() + bar.get_width() / 2.0, height + offset, f"r={r:.2f}\n({p_str})", ha='center', va=va, fontsize=7.5, fontweight=fontweight, color='black', zorder=4)

        # --------annotating plot------------------
        ax.axhline(0, color='black', linewidth=1.0, zorder=2)
        ax.set_ylabel(ylabel_str, fontsize=args.fontsize)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, prune='lower'))
        ax.tick_params(axis='y', which='both', labelsize=args.fontsize)
        ax.set_ylim(-1.1, 1.1)
        ax.grid(axis='y', linestyle=':', alpha=0.6, zorder=0)

        ax.text(0.5, 0.98, panel_title, transform=ax.transAxes, ha='center', va='top', fontsize=args.fontsize / args.fontfactor, fontweight='bold', zorder=5)

    # -------shared x-axis tick marks & labels---------
    ax1.tick_params(axis='x', which='both', bottom=True, labelbottom=False) # Top panel
    ax2.tick_params(axis='x', which='both', bottom=True, labelbottom=True) # Bottom panel
    ax2.set_xticks(x)
    ax2.set_xticklabels([item[2] for item in items_partial], fontsize=args.fontsize)

    # --------making significance legend (left on top panel)---------
    if args.nocolorcoding:
        legend_sig = [patches.Patch(facecolor=sig_color, edgecolor='black', label=r'$p \leq 0.05$ (Significant)'), patches.Patch(facecolor=nonsig_color, edgecolor='black', label=r'$p > 0.05$ (Not Significant)')]
    else:
        legend_sig = [patches.Patch(facecolor=cfg['color'], edgecolor='black', label=cfg['label'])for cfg in p_value_color_config]
    leg1 = ax1.legend(handles=legend_sig, loc='best', frameon=True, fontsize=args.fontsize / args.fontfactor)
    #ax1.add_artist(leg1)

    # -------making mass legend (right on top panel)-------------------
    legend_regime = [patches.Patch(facecolor='white', edgecolor='black', hatch=h, label=lbl) for h, lbl in zip(regime_hatches, regime_labels)]
    ax2.legend(handles=legend_regime, loc='best', frameon=True, fontsize=args.fontsize / args.fontfactor)

    # -------save fig--------
    figname = f'correlation_summary_{qualifiers}.png'
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

p_value_color_config = [
    {'max_p': 0.01, 'color': '#7b241c', 'label': r'$p \leq 0.01$ (Very Significant)'},
    {'max_p': 0.05, 'color': '#cd6155', 'label': r'$0.01 < p \leq 0.05$ (Significant)'},
    {'max_p': 0.10, 'color': "#7b7c7d", 'label': r'$0.05 < p \leq 0.10$ (Marginally Significant)'},
    {'max_p': 1.00, 'color': '#ebedef', 'label': r'$p > 0.10$ (Not Significant)'},
]

log_mass_cut = 9.0

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.4
    
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
    
    results_dict = compute_correlations(df_grad, args, xcol='log_mass_median', ycol='radial_logOH_grad', qualifiers=qualifiers)
    
    #fig_mzr = plot_stacked_MZR(df_grad, args, xcol='log_mass_median', ycol='logOH_int', colorcol=colorcol, qualifiers=qualifiers, cmap=cmap)
    #fig_mzgr = plot_stacked_MZGR(df_grad, args, xcol='log_mass_median', ycol='radial_logOH_grad', colorcol=colorcol, qualifiers=qualifiers, cmap=cmap)
    fig_corr = plot_correlation_summary_hatched(results_dict, args, qualifiers=qualifiers)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
