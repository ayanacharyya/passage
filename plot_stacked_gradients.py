'''
    Filename: plot_stacked_gradients.py
    Notes: Reads stacked (in bins of mass and/or SFR) 2D emission line maps for objects in a given field/s, computes metallicity gradient, and plots them
    Author : Ayan
    Created: 18-01-26
    Example: run plot_stacked_gradients.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run plot_stacked_gradients.py --field Par28 --Zdiag R23 --use_C25 --adaptive_bins --overplot_literature --overplot_passage
             run plot_stacked_gradients.py --field Par28 --Zdiag R23 --use_C25 --adaptive_bins --fold_maps --overplot_literature --overplot_passage
             run plot_stacked_gradients.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_minor_major_profile
             run plot_stacked_gradients.py --field Par28 --Zdiag R23 --use_C25 --fold_maps
             run plot_stacked_gradients.py --system ssd --do_all_fields --Zdiag R23 --use_C25 --adaptive_bins --bin_by_distance --overplot_literature --overplot_passage --fold_maps
'''

from header import *
from util import *
from make_sfms_bins import log_mass_bins, log_sfr_bins, read_passage_sed_catalog, bin_SFMS_linear, bin_SFMS_adaptive, read_passage_sed_catalog, bin_SFMS_distance, get_sfms_func, n_adaptive_bins, sfms
from make_passage_plots import plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def read_stacked_df(filename):
    '''
    Reads the stacked gradient list from a fits Table
    Returns pandaas dataframe
    '''
    df = Table.read(filename).to_pandas()
    print(f'Reading in {filename}..')
    for thiscol in ['log_mass_bin', 'log_sfr_bin', 'bin_intervals']:
        if thiscol in df:
            df[thiscol] = df[thiscol].str.decode('utf-8')
            df[thiscol] = pd.IntervalIndex.from_tuples(df[thiscol].apply(lambda x: pd.to_numeric(x.strip('()[]').split(', '))).map(tuple), closed='right')

    return df

# --------------------------------------------------------------------------------------------------------------------
def rescale_for_seaborn(xarr, yarr, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1):
    '''
    Rescale given stellar mass and SFR arrays to the limits of an exisiting seaborn heatmap, so they can be overplotted
    Returns rescaled x and y arrays (which can then be overplotted directly)
    '''
    # -----transformating to seaborn phase space---
    x_plot = ((xarr - xmin) / (xmax - xmin)) * nx #- 0.5
    y_plot = ((yarr - ymin) / (ymax - ymin)) * ny #- 0.5
    y_plot_reversed = ny - y_plot #- 1 # to reverse the SFR axis, as it is in the heatmap

    return x_plot, y_plot_reversed

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Whitaker14_sns(ax, redshift, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='olivegreen'):
    '''
    Overplots fitted SFMS based on Whitaker+14 (https://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdf) Table 1, eq 2
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    log_mass_low_lim = 9.2 # lower limit of mass they fitted up to
    if redshift >= 0.5 and redshift < 1.0:
        a, b, c, z1, z2 = -27.40, 5.02, -0.22, 0.5, 1.0 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 1.0 and redshift < 1.5:
        a, b, c, z1, z2 = -26.03, 4.62, -0.19, 1.0, 1.5 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 1.5 and redshift < 2.0:
        a, b, c, z1, z2 = -24.04, 4.17, -0.16, 1.5, 2.0 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 2.0 and redshift < 2.5:
        a, b, c, z1, z2 = -19.99, 3.44, -0.13, 2.0, 2.5 # polynomial coefficients from Table 1, for this z range
    else:
        print(f'Provided redshift {redshift} should be between 0.5 and 2.5, otherwise cannot plot')
        return ax

    log_mass1 = np.linspace(xmin, log_mass_low_lim, 20)
    log_SFR1 = np.poly1d([c, b, a])(log_mass1)

    log_mass2 = np.linspace(log_mass_low_lim, xmax, 20)
    log_SFR2 = np.poly1d([c, b, a])(log_mass2)

    log_mass1, log_SFR1 = rescale_for_seaborn(log_mass1, log_SFR1, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)

    log_mass2, log_SFR2 = rescale_for_seaborn(log_mass2, log_SFR2, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Whitaker+14: {z1}<z<{z2}')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Shivaei15_sns(ax, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='salmon'):
    '''
    Overplots fitted SFMS based on Shivaei+15 (https://iopscience.iop.org/article/10.1088/0004-637X/815/2/98/pdf) Table 1
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a1, b1, c1 = 0.65, -5.40, 0.40 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 1.37-2.61
    a2, b2, c2 = 0.58, -4.65, 0.36 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 2.09-2.61
    log_mass_low_lim = 9.5 # lower limit of mass they fitted up to

    log_mass1 = np.linspace(xmin, log_mass_low_lim, 20)
    log_SFR1 = a1 * log_mass1 + b1

    log_mass2 = np.linspace(log_mass_low_lim, xmax, 20)
    log_SFR2 = a1 * log_mass2 + b1

    log_mass1, log_SFR1 = rescale_for_seaborn(log_mass1, log_SFR1, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.fill_between(log_mass1, log_SFR1 - c1/2, log_SFR1 + c1/2, alpha=0.3, facecolor=color)

    log_mass2, log_SFR2 = rescale_for_seaborn(log_mass2, log_SFR2, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Shivaei+15: z~2')
    ax.fill_between(log_mass2, log_SFR2 - c1/2, log_SFR2 + c1/2, alpha=0.3, facecolor=color)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso23_sns(ax, redshift, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a0, a1, b0, b1, b2 = ufloat(0.2,0.02), ufloat(-0.034, 0.002), ufloat(-26.134,0.015), ufloat(4.722, 0.012), ufloat(-0.1925, 0.0011)  # Table 2, Eq 10
    log_mass_low_lim = 8.7 # lower limit of mass they fitted up to

    age_at_z = cosmo.age(redshift).value # Gyr

    log_mass1 = np.linspace(xmin, log_mass_low_lim, 20)
    log_SFR1 = (a1 * age_at_z + b1) * log_mass1 + b2 * (log_mass1) ** 2 + b0 + a0 * age_at_z

    log_mass2 = np.linspace(log_mass_low_lim, xmax, 20)
    log_SFR2 = (a1 * age_at_z + b1) * log_mass2 + b2 * (log_mass2) ** 2 + b0 + a0 * age_at_z

    log_mass1, log_SFR1 = rescale_for_seaborn(log_mass1, log_SFR1, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass1, unp.nominal_values(log_SFR1), ls='dashed', c=color, lw=2)
    ax.fill_between(log_mass1, unp.nominal_values(log_SFR1) - unp.std_devs(log_SFR1)/2, unp.nominal_values(log_SFR1) + unp.std_devs(log_SFR1)/2, alpha=0.3, facecolor=color)

    log_mass2, log_SFR2 = rescale_for_seaborn(log_mass2, log_SFR2, nx, ny, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.plot(log_mass2, unp.nominal_values(log_SFR2), ls='solid', c=color, lw=2, label=f'Popesso+23: z = {redshift}')
    ax.fill_between(log_mass2, unp.nominal_values(log_SFR2) - unp.std_devs(log_SFR2)/2, unp.nominal_values(log_SFR2) + unp.std_devs(log_SFR2)/2, alpha=0.3, facecolor=color)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_sns(ax, pivot_color, pivot_annot, args, cmap='viridis', clabel='', cmin=None, cmax=None, ncbins=4):
    '''
    Makes heatmap from a pivot table, using seaborn, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''    
    log_mass_edges = [int.left for int in pivot_color.columns] + [pivot_color.columns[-1].right]
    log_sfr_edges = [int.right for int in pivot_color.index] + [pivot_color.index[-1].left]

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    
    pivot_mask = pivot_color.isnull()
    hm = sns.heatmap(pivot_color, ax=ax, annot=pivot_annot, fmt='.0f', cmap=cmap, cbar_ax=cax, xticklabels=log_mass_edges, yticklabels=log_sfr_edges, mask=pivot_mask, vmin=cmin, vmax=cmax) # This hides the NaN bins

    # -------to draw lines around cells in use----------
    for y in range(pivot_mask.shape[0]):
        for x in range(pivot_mask.shape[1]):
            if not pivot_mask.iloc[y, x]: # If NOT masked (has data)
                ax.add_patch(plt.Rectangle((x, y), 1, 1, fill=False, edgecolor='k', lw=0.5))
    
    # ---------annotating colorbar------------
    cbar = hm.collections[0].colorbar
    cbar.set_label(clabel, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()
    
    # --------annotating axis borders-----------------
    ax.set_xlabel(r'$\log$ Stellar Mass [M$_\odot$]', fontsize=args.fontsize)
    ax.set_xticks(np.arange(len(log_mass_edges)))
    ax.set_xticklabels(log_mass_edges)
    ax.set_xlim(-0.2, len(pivot_color.columns) + 0.2)
    
    ax.set_ylabel(r'$\log$ SFR [M$_\odot$/yr]', fontsize=args.fontsize)
    ax.set_yticks(np.arange(len(log_sfr_edges)))
    ax.set_yticklabels(log_sfr_edges)
    ax.set_ylim(len(pivot_color.index) + 0.2, -0.2)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize / args.fontfactor, labelbottom=True)

    # -----drawing hatches on the missing cells-------------
    # ax.set_facecolor('k' if args.fortalk else 'whitesmoke')
    # matplotlib.rcParams['hatch.linewidth'] = 0.1
    # missing_data_mask = pivot_color.isnull()
    # ax.pcolor(np.where(missing_data_mask, 0, np.nan), hatch='///', alpha=0.) # 'hatch' options: '/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*'

    # ----seaborn heatmaps hide axis edges by default-------
    for this_ax in [ax, cax]:
        for spine in this_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.)
            spine.set_edgecolor('black')
    
    # ----------over-plotting theoretical diagrams----------
    if args.overplot_literature:
        #ax = plot_SFMS_Whitaker14_sns(ax, 2, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='yellowgreen')
        ax = plot_SFMS_Shivaei15_sns(ax, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='darkgreen')
        ax = plot_SFMS_Popesso23_sns(ax, 2, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='darkgoldenrod')
        #ax = plot_SFMS_Popesso23_sns(ax, 3, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='royalblue')
        ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_heatmap_sns(df, args, quant='logOH'):
    '''
    Makes a nice heatmap (with seaborn) of stacked integrated metallicities and metallicity gradients
    Returns figure handle
    '''
    # -----------------preparing the pivot tables---------------
    log_mass_intervals = [pd.Interval(log_mass_bins[i], log_mass_bins[i+1], closed='left') for i in range(len(log_mass_bins)-1)]
    log_sfr_intervals = [pd.Interval(log_sfr_bins[i], log_sfr_bins[i+1], closed='left') for i in range(len(log_sfr_bins)-1)]
    log_sfr_intervals = log_sfr_intervals[::-1] # reverse the Y-axis so SFR increases upwards

    pivot_int = df.pivot(index='log_sfr_bin', columns='log_mass_bin', values=f'{quant}_int')
    pivot_grad = df.pivot(index='log_sfr_bin', columns='log_mass_bin', values=f'{quant}_grad')
    pivot_nobj = df.pivot(index='log_sfr_bin', columns='log_mass_bin', values='nobj')

    pivot_int = pivot_int.reindex(index=log_sfr_intervals, columns=log_mass_intervals)
    pivot_grad = pivot_grad.reindex(index=log_sfr_intervals, columns=log_mass_intervals)
    pivot_nobj = pivot_nobj.reindex(index=log_sfr_intervals, columns=log_mass_intervals)

    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    fig.subplots_adjust(left=0.07, right=0.92, top=0.95, bottom=0.13, wspace=0.3, hspace=0.)
    
    # -----------------plot integrated metallicity heatmap---------------
    axes[0] = make_heatmap_sns(axes[0], pivot_int, pivot_nobj, args, cmap='plasma', clabel=r'$\log$(O/H) + 12', cmin=7.0, cmax=7.5, ncbins=5)

    # -----------------plot metallicity gradient heatmap---------------
    axes[1] = make_heatmap_sns(axes[1], pivot_grad, pivot_nobj, args, cmap='coolwarm', clabel=r'$\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4)

    # ---------overplot PASSAGE galaxies (integrated stellar mass-SFR)--------------------
    if args.overplot_passage:
        df_photcat = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas() # read the photometric catalog file
        df_photcat['field'] = args.field
        df_photcat = get_passage_masses_from_cosmos(df_photcat, args, id_col='id') # crossmatch with cosmos-web to get stellar mass and SFR

        log_mass, log_sfr = df_photcat['log_mass'], df_photcat['log_sfr']
        log_mass, log_sfr = rescale_for_seaborn(log_mass, log_sfr, len(log_mass_intervals), len(log_sfr_intervals), xmin=log_mass_intervals[0].left, xmax=log_mass_intervals[-1].right, ymin=log_sfr_intervals[-1].left, ymax=log_sfr_intervals[0].right)
        
        col, edgecol = 'w', 'sienna'
        for ax in axes:
            ax.scatter(log_mass, log_sfr, s=5, c=col, lw=1, edgecolors=edgecol, label=f'{args.field}')
            ax.legend(fontsize=args.fontsize / args.fontfactor, loc='upper left')

    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_patches(ax, df, quant, args, xcolname='log_mass_bin', ycolname='log_sfr_bin', cmap='viridis', clabel='', cmin=None, cmax=None, ncbins=4):
    '''
    Makes heatmap from a dataframe, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
    cmap = mplcolormaps[cmap]
    if cmin is None: cmin = df[quant].min()
    if cmax is None: cmax = df[quant].max()
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)

    for index, row in df.iterrows():
        if np.isnan(row[quant]): continue # to prevent blank boxes with borders, for nan values
        m_min, m_max = row[xcolname].left, row[xcolname].right
        s_min, s_max = row[ycolname].left, row[ycolname].right
        width = m_max - m_min
        height = s_max - s_min
        
        color = cmap(norm(row[quant]))
        rect = patches.Rectangle((m_min, s_min), width, height, facecolor=color, edgecolor='w' if args.fortalk else 'k', linewidth=0.5)
        ax.add_patch(rect)
        
        # -----annotating bins with nobjects--------
        if 'int' in quant:
            if row[quant] > cmin + 0.5 * (cmax - cmin): text_color = 'k'
            else: text_color = 'w'
        else:
            if abs(row[quant]) < 0.3: text_color = 'k'
            else: text_color = 'w'
        if row['nobj'] > 0: ax.text(m_min + width/2, s_min + height/2, int(row['nobj']), ha='center', va='center', fontsize=args.fontsize / args.fontfactor, color=text_color)

    # --------annotating axis borders-----------------
    ax.set_xlabel(r'$\log$ Stellar Mass [M$_\odot$]', fontsize=args.fontsize)
    ax.set_xlim(log_mass_bins[0] -0.2, log_mass_bins[-1] + 0.2)
    
    ax.set_ylabel(r'$\log$ SFR [M$_\odot$/yr]', fontsize=args.fontsize)
    ax.set_ylim(log_sfr_bins[0] -0.2, log_sfr_bins[-1] + 0.2)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize / args.fontfactor, labelbottom=True)

    # ---------annotating colorbar------------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, ax=ax, label=clabel)
    cbar.set_label(clabel, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()
    
    return ax

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_distance(ax, df, sfms, quant, args, method_text='', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe that has been binned by distance from the SFMS, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
    # ---------obtaining sfms function and color mappable---------
    sfms_func = get_sfms_func(log_mass_bins, sfms)
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
   
    # -----------defining vertices, for annotating and color-coding the bins-------------
    for index, row in df.iterrows():
        interval = row[f'bin_intervals{method_text}']
        m_grid = np.linspace(row['log_mass_min'], row['log_mass_max'], 50)
        sfms_line = sfms_func(m_grid)
        
        color = sm.to_rgba(row[quant])
        ax.fill_between(m_grid, sfms_line + interval.left, sfms_line + interval.right, color=color, alpha=0.8, edgecolor='w' if args.fortalk else 'k', lw=0.5)
        
        s_center = sfms_func(row['log_mass_median']) + (interval.left + interval.right) / 2
        if ~np.isnan(row['nobj']): ax.text(row['log_mass_median'], s_center, int(row['nobj']), color='w' if args.fortalk else 'k', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold', rotation=45)
 
    # --------annotating axis borders-----------------
    ax.set_xlim(log_mass_bins.min() -0.2, log_mass_bins.max() + 0.2)
    ax.set_ylim(log_sfr_bins.min() -0.2, log_sfr_bins.max()+ 0.2)

    if hide_xaxis:
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize, labelbottom=False)
    else:
        ax.set_xlabel(r'$\log$ Stellar Mass [M$_\odot$]', fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize, labelbottom=True)

    if hide_yaxis:
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize, labelleft=False)
    else:
        ax.set_ylabel(r'$\log$ SFR [M$_\odot$/yr]', fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize, labelleft=True)

    # ---------annotating colorbar------------
    if not hide_cbar:
        if cmin is None: cmin = df[quant].min()
        if cmax is None: cmax = df[quant].max()
        norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ax=ax, label=clabel)
        cbar.set_label(clabel, fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)
        cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
        cbar.update_ticks()
    
    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_heatmap_patches(df, args, quant='logOH'):
    '''
    Makes a nice heatmap (with patches) of stacked integrated metallicities and metallicity gradients
    Returns figure handle
    '''
    if quant == 'logOH': df = df[df[f'radial_{quant}_grad'] > -2] # removing spurious values that are clearly bogus
    
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.13, wspace=0.2, hspace=0.)
    
    # ---------plot the heatmaps-------------------
    if args.bin_by_distance:
        if args.plot_minor_major_profile:
            axes[0] = make_heatmap_distance(axes[0], df, sfms, f'minor_{quant}_grad', args, cmap='coolwarm', clabel=r'Minor $\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4) # plot integrated metallicity heatmap
            axes[1] = make_heatmap_distance(axes[1], df, sfms, f'major_{quant}_grad', args, cmap='coolwarm', clabel=r'Major $\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4) # plot metallicity gradient heatmap
        else:
            axes[0] = make_heatmap_distance(axes[0], df, sfms, f'{quant}_int', args, cmap='plasma', clabel=r'$\log$(O/H) + 12', cmin=7.0, cmax=8., ncbins=5) # plot integrated metallicity heatmap
            axes[1] = make_heatmap_distance(axes[1], df, sfms, f'radial_{quant}_grad', args, cmap='coolwarm', clabel=r'$\nabla$Z$_r$ [dex/R$_e$]', cmin=-1.5, cmax=1.5, ncbins=4) # plot metallicity gradient heatmap
    else:
        if args.plot_minor_major_profile:
            axes[0] = make_heatmap_patches(axes[0], df, f'minor_{quant}_grad', args, cmap='coolwarm', clabel=r'Minor $\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4) # plot minor metallicity gradient heatmap
            axes[1] = make_heatmap_patches(axes[1], df, f'major_{quant}_grad', args, cmap='coolwarm', clabel=r'Major $\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4) # plot major metallicity gradient heatmap
        else:
            axes[0] = make_heatmap_patches(axes[0], df, f'{quant}_int', args, cmap='plasma', clabel=r'$\log$(O/H) + 12', cmin=7.0, cmax=7.5, ncbins=5) # plot integrated metallicity heatmap
            axes[1] = make_heatmap_patches(axes[1], df, f'radial_{quant}_grad', args, cmap='coolwarm', clabel=r'$\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4) # plot metallicity gradient heatmap

    # ---------overplot PASSAGE galaxies (integrated stellar mass-SFR)--------------------
    if args.overplot_passage:
        if args.do_all_fields: # reading in the master SED catalog
            passage_catalog_filename = args.output_dir / 'catalogs' / 'passagepipe_v0.5_SED_fits_cosmosweb_v1.0.0-alpha.fits'
            df_photcat = read_passage_sed_catalog(passage_catalog_filename)
        else: # reading in the single-field phot catalog
            product_dir = args.input_dir / args.field / 'Products'
            df_photcat = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
            df_photcat['field'] = args.field
            df_photcat = get_passage_masses_from_cosmos(df_photcat, args, id_col='id') # crossmatch with cosmos-web to get stellar mass and SFR
            
        col, edgecol = 'w', 'sienna'
        for ax in axes:
            ax.scatter(df_photcat['log_mass'], df_photcat['log_sfr'], s=5, c=col, lw=1, edgecolors=edgecol, label=f'{args.field}')
        axes[0].legend(fontsize=args.fontsize / args.fontfactor, loc='upper left')

    plt.show(block=False)

    # ----------over-plotting theoretical diagrams----------
    if args.overplot_literature:
        for ax in axes:
            if sfms == 'Whitaker14': ax = plot_SFMS_Whitaker14(ax, 2, color='crimson')
            elif sfms == 'Shivaei15': ax = plot_SFMS_Shivaei15(ax, color='royalblue')
            elif sfms == 'Popesso23': ax = plot_SFMS_Popesso23(ax, 2, color='darkgreen')
            #elif sfms == 'Popesso23': ax = plot_SFMS_Popesso23(ax, 3, color='royalblue')
            else: raise ValueError(f'Method {sfms} not found in the list of available SFMS literature methods: Whitaker14, Shivaei15, Popesso23')
        axes[0].legend(fontsize=args.fontsize / args.fontfactor, loc='upper left')

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    fold_text = '_folded' if args.fold_maps else ''
    adapt_text = '_adaptivebins' if args.adaptive_bins else ''
    deproject_text = '_nodeproject' if args.skip_deproject else ''
    rescale_text = '_norescale' if args.skip_re_scaling else ''
    C25_text = '_wC25' if args.use_C25 and 'NB' not in args.Zdiag else ''
    if args.bin_by_distance: binby_text = f'_binby_distance{n_adaptive_bins}'
    elif args.adaptive_bins: binby_text = f'_bin_adaptive_max{args.max_gal_per_bin}'
    else: binby_text = '_bin_linear'

    # ---------determining list of fields----------------
    if args.do_all_fields:
        passage_catalog_filename = args.output_dir / 'catalogs' / 'passagepipe_v0.5_SED_fits_cosmosweb_v1.0.0-alpha.fits'
        df = read_passage_sed_catalog(passage_catalog_filename)
        output_dir = args.output_dir / 'stacking'
    else:
        product_dir = args.input_dir / args.field / 'Products'
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        df['field'] = args.field
        df = get_passage_masses_from_cosmos(df, args, id_col='id') # crossmatch with cosmos-web to get stellar mass and SFR
        output_dir = args.output_dir / args.field / 'stacking'

    # -------------binning the mass-SFR plane-------------
    if args.adaptive_bins:
        if args.bin_by_distance: df, bin_list = bin_SFMS_distance(df, method_text='', n_adaptive_bins=n_adaptive_bins, sfms=sfms)
        else: df, bin_list = bin_SFMS_adaptive(df, method_text='', max_n=args.max_gal_per_bin) # binning the dataframe in an adaptive way
    else:
        if args.bin_by_distance: df, bin_list = bin_SFMS_distance(df, method_text='', delta_bin=0.2, sfms=sfms)
        df, bin_list = bin_SFMS_linear(df, method_text='') # -binning the dataframe uniformly by mass and SFR bins

    if args.bin_by_distance: bin_list = np.sort(bin_list)
    else: bin_list.sort(key=lambda x: (x[0].left, x[1].left))

    # ----------getting directory structure----------
    if args.adaptive_bins: output_dir = Path(str(output_dir).replace('stacking', 'stacking_adaptive'))
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / f'plots{deproject_text}{rescale_text}'
    fig_dir.mkdir(parents=True, exist_ok=True)
    fits_dir = output_dir / f'maps{deproject_text}{rescale_text}'
    fits_dir.mkdir(parents=True, exist_ok=True)

    grad_filename = fits_dir / f'stacked{binby_text}{fold_text}_fits_allbins_Zdiag_{args.Zdiag}{C25_text}{deproject_text}{rescale_text}.fits'

    # -------------reading in stacked gradient dataframe-----------------------
    df_grad = read_stacked_df(grad_filename)
    if args.bin_by_distance:
        df = df.groupby(['bin_intervals']).agg(log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median')).reset_index()
        df_grad = pd.merge(df_grad, df, on=['bin_intervals'], how='left')
    
    # ------------plotting stacked gradients on SFMS--------------------------
    #fig = plot_SFMS_heatmap_sns(df_grad, args)
    fig = plot_SFMS_heatmap_patches(df_grad, args)
    save_fig(fig, fig_dir, f'stacked{adapt_text}{binby_text}{fold_text}_Zdiag_{args.Zdiag}{C25_text}{deproject_text}{rescale_text}_SFMS_heatmap.png', args) # saving the figure


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
