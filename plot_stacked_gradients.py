'''
    Filename: plot_stacked_gradients.py
    Notes: Reads stacked (in bins of mass and/or SFR) 2D emission line maps for objects in a given field/s, computes metallicity gradient, and plots them
    Author : Ayan
    Created: 18-01-26
    Example: run plot_stacked_gradients.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run plot_stacked_gradients.py --field Par28 --Zdiag R23 --use_C25
'''

from header import *
from util import *
from stack_emission_maps import log_mass_bins, log_sfr_bins, get_passage_masses_from_cosmos
from plot_stacked_maps import save_fig

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def read_stacked_df(filename):
    '''
    Reads the stacked gradient list from a fits Table
    Returns pandaas dataframe
    '''
    df = Table.read(filename).to_pandas()
    print(f'Reading in {filename}..')
    for thiscol in ['log_mass_bin', 'log_sfr_bin']:
        df[thiscol] = df[thiscol].str.decode('utf-8')
        df[thiscol] = pd.IntervalIndex.from_tuples(df[thiscol].apply(lambda x: pd.to_numeric(x.strip('()[]').split(', '))).map(tuple), closed='left')

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
def plot_SFMS_Whitaker14(ax, redshift, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='olivegreen'):
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
def plot_SFMS_Shivaei15(ax, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='salmon'):
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
def plot_SFMS_Popesso23(ax, redshift, nx, ny, xmin=7, xmax=11, ymin=0, ymax=1, color='cornflowerblue'):
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
def make_heatmap(ax, pivot_color, pivot_annot, args, cmap='viridis', clabel='', cmin=None, cmax=None, ncbins=4):
    '''
    Makes heatmap from a pivot table, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''    
    log_mass_edges = [int.left for int in pivot_color.columns] + [pivot_color.columns[-1].right]
    log_sfr_edges = [int.right for int in pivot_color.index] + [pivot_color.index[-1].left]

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    
    hm = sns.heatmap(pivot_color, ax=ax, annot=pivot_annot, fmt='.0f', cmap=cmap, cbar_ax=cax, xticklabels=log_mass_edges, yticklabels=log_sfr_edges, mask=pivot_color.isnull(), vmin=cmin, vmax=cmax) # This hides the NaN bins

    cbar = hm.collections[0].colorbar
    cbar.set_label(clabel, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()
    
    ax.set_xlabel(r'$\log$ Stellar Mass [M$_\odot$]', fontsize=args.fontsize)
    ax.set_xticks(np.arange(len(log_mass_edges)))
    ax.set_xticklabels(log_mass_edges)
    ax.set_xlim(0, len(pivot_color.columns))
    
    ax.set_ylabel(r'$\log$ SFR [M$_\odot$/yr]', fontsize=args.fontsize)
    ax.set_yticks(np.arange(len(log_sfr_edges)))
    ax.set_yticklabels(log_sfr_edges)
    ax.set_ylim(len(pivot_color.index), 0)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize / args.fontfactor, labelbottom=True)

    # -----drawing hatches on the missing cells-------------
    ax.set_facecolor('whitesmoke')
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
    #ax = plot_SFMS_Whitaker14(ax, 2, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='yellowgreen')
    ax = plot_SFMS_Shivaei15(ax, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='darkgreen')
    ax = plot_SFMS_Popesso23(ax, 2, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='darkgoldenrod')
    #ax = plot_SFMS_Popesso23(ax, 3, len(pivot_color.columns), len(pivot_color.index), xmin=log_mass_edges[0], xmax=log_mass_edges[-1], ymin=log_sfr_edges[-1], ymax=log_sfr_edges[0], color='royalblue')
    ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_heatmap(df, args, quant='logOH'):
    '''
    Makes a nice heatmap of stacked integrated metallicities and metallicity gradients
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
    axes[0] = make_heatmap(axes[0], pivot_int, pivot_nobj, args, cmap='plasma', clabel=r'$\log$(O/H) + 12', cmin=7.0, cmax=7.5, ncbins=5)

    # -----------------plot metallicity gradient heatmap---------------
    axes[1] = make_heatmap(axes[1], pivot_grad, pivot_nobj, args, cmap='coolwarm', clabel=r'$\nabla$Z$_r$ [dex/R$_e$]', cmin=-1., cmax=1., ncbins=4)

    # ---------overplot PASSAGE galaxies (integrated stellar mass-SFR)--------------------
    df_photcat = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas() # read the photometric catalog file
    df_photcat['field'] = args.field
    df_photcat = get_passage_masses_from_cosmos(df_photcat, args, id_col='id') # crossmatch with cosmos-web to get stellar mass and SFR

    log_mass, log_sfr = df_photcat['log_mass'], df_photcat['log_sfr']
    log_mass, log_sfr = rescale_for_seaborn(log_mass, log_sfr, len(log_mass_intervals), len(log_sfr_intervals), xmin=log_mass_intervals[0].left, xmax=log_mass_intervals[-1].right, ymin=log_sfr_intervals[-1].left, ymax=log_sfr_intervals[0].right)
    
    col, edgecol = 'w', 'sienna'
    for ax in axes:
        ax.scatter(log_mass, log_sfr, s=5, c=col, lw=1, edgecolors=edgecol, label=f'{args.field}')
        ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')
        #ax.set_aspect((log_sfr_intervals[0].right - log_sfr_intervals[-1].left) / (log_mass_intervals[-1].right - log_mass_intervals[0].left))
        #ax.set_aspect(0.2)

    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / 'Par*') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr
    
    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'\nStarting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / args.field / 'stacking'
        output_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir / 'plots'
        fig_dir.mkdir(parents=True, exist_ok=True)
        fits_dir = output_dir / 'maps'
        fits_dir.mkdir(parents=True, exist_ok=True)

        C25_text = '_wC25' if args.use_C25 and 'NB' not in args.Zdiag else ''
        grad_filename = fits_dir / f'stacked_fits_allbins_Zdiag_{args.Zdiag}{C25_text}.fits'

        # -------------reading in stacked gradient dataframe-----------------------
        df_grad = read_stacked_df(grad_filename)
        
        # ------------plotting stacked gradients on SFMS--------------------------
        fig = plot_SFMS_heatmap(df_grad, args)
        save_fig(fig, fig_dir, f'stacked_SFMS_heatmap.png', args) # saving the figure

        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
