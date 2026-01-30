'''
    Filename: make_sfms_bins.py
    Notes: Bins the SFR-mass plane in various ways, for testing to decide which type of binning has the best optimisation between more bins and more galaxies per bin
    Author : Ayan
    Created: 30-01-26
    Example: run make_sfms_bins.py --field Par028 --overplot_literature --overplot_passage
             run make_sfms_bins.py --field Par028
'''

from header import *
from util import *
from stack_emission_maps import log_mass_bins, log_sfr_bins, get_passage_masses_from_cosmos
from plot_stacked_maps import save_fig
from make_passage_plots import plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_patches(ax, df, quant, args, xcolname='log_mass_bin', ycolname='log_sfr_bin', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
    cmap = mplcolormaps[cmap]

    if 'm_min' not in df:
        df['m_min'] = df[xcolname].array.left
        df['m_max'] = df[xcolname].array.right

        df['s_min'] = df[ycolname].array.left
        df['s_max'] = df[ycolname].array.right
    
    # -----------defining vertices for bins-------------
    patches, values = [], []
    for _, row in df.iterrows():
        vertices = [(row['m_min'], row['s_min']), (row['m_max'], row['s_min']), (row['m_max'], row['s_max']), (row['m_min'], row['s_max'])]
        patches.append(Polygon(vertices, closed=True))
        values.append(row[quant])

    p = PatchCollection(patches, cmap=cmap, edgecolors='w', alpha=0.8)
    p.set_array(np.array(values))
    ax.add_collection(p)

    # --------annotate at the center of each patch------------
    for i, row in df.iterrows():
        if row[quant] > 0:
            cx = (row['m_min'] + row['m_max']) / 2
            cy = (row['s_min'] + row['s_max']) / 2
            ax.text(cx, cy, int(row[quant]), color='black', ha='center', va='center', fontsize=9, fontweight='bold')
    
    # --------annotating axis borders-----------------
    ax.set_xlim(df['m_min'].min() -0.2, df['m_max'].max() + 0.2)
    ax.set_ylim(df['s_min'].min() -0.2, df['s_max'].max()+ 0.2)

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
    
    # ----seaborn heatmaps hide axis edges by default-------
    for this_ax in [ax]:
        for spine in this_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.)
            spine.set_edgecolor('black')

    # ----------over-plotting theoretical diagrams----------
    if args.overplot_literature:
        #ax = plot_SFMS_Whitaker14(ax,2,  color='yellowgreen')
        ax = plot_SFMS_Shivaei15(ax, color='darkgreen')
        ax = plot_SFMS_Popesso23(ax, 2, color='darkgoldenrod')
        #ax = plot_SFMS_Popesso23(ax, 3, color='royalblue')
        ax.legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_bins(df, methods, args):
    '''
    Makes a nice heatmap (with patches) of stacked integrated metallicities and metallicity gradients
    Returns figure handle
    '''
    cmap, cmin, cmax, ncbins, clabel ='viridis', 0, 100, 5, 'Number of galaxies'
    ncol = len(methods)
    
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, ncol, figsize=(5.2 * ncol, 5.), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.13, wspace=0., hspace=0.)
    
    # ---------plot the heatmaps-------------------
    for index, method in enumerate(methods):
        df_sub = df.groupby(f'bin_intervals_{method}').size().reset_index(name='n_galaxies')
        df_sub[['log_mass_bin', 'log_sfr_bin']] = pd.DataFrame(df_sub[f'bin_intervals_{method}'].tolist(), index=df_sub.index)
        axes[index] = make_heatmap_patches(axes[index], df_sub, 'n_galaxies', args, xcolname=f'log_mass_bin', ycolname=f'log_sfr_bin', cmap=cmap, hide_cbar=True, hide_yaxis=index)
        axes[index].text(0.05, 0.95, f'{method}', ha='left', va='top', c='k', fontsize=args.fontsize, transform=axes[index].transAxes)
        axes[index].set_aspect('equal')

        if args.overplot_passage: axes[index].scatter(df['log_mass'], df['log_sfr'], s=5, c='w', lw=1, edgecolors='sienna', label=f'{args.field}') # overplot PASSAGE galaxies (integrated stellar mass-SFR)

    # ---------annotating colorbar------------
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    cbar = fig.colorbar(sm, ax=axes, location='top', shrink=0.95, pad=0.01, aspect=60)
    cbar.set_label(clabel, fontsize=args.fontsize, labelpad=5)    
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()

    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def get_adaptive_bins(df_subset, m_range, s_range, max_n=20):
    '''
    m_range: (min, max) of log_mass for this specific tile
    s_range: (min, max) of log_sfr for this specific tile
    Courtesy of this function: Gemini
    '''
    # Count how many galaxies are in this specific rectangular area
    mask = (df_subset['log_mass'] >= m_range[0]) & (df_subset['log_mass'] < m_range[1]) & \
           (df_subset['log_sfr'] >= s_range[0]) & (df_subset['log_sfr'] < s_range[1])
    
    subset = df_subset[mask]
    n_count = len(subset)

    # Base case: if count is small OR area is already very tiny, stop splitting
    if n_count <= max_n or (m_range[1] - m_range[0]) < 0.1:
        if n_count == 0: return []
        # Return the coordinates and the mean value for this leaf node
        return [{'m_min': m_range[0], 'm_max': m_range[1], 's_min': s_range[0], 's_max': s_range[1], 'n_count': n_count}]
    
    # Recursive step: Split into 4 quadrants
    m_mid = (m_range[0] + m_range[1]) / 2
    s_mid = (s_range[0] + s_range[1]) / 2
    
    results = []
    results.extend(get_adaptive_bins(subset, (m_range[0], m_mid), (s_range[0], s_mid))) # Bottom-Left
    results.extend(get_adaptive_bins(subset, (m_mid, m_range[1]), (s_range[0], s_mid))) # Bottom-Right
    results.extend(get_adaptive_bins(subset, (m_range[0], m_mid), (s_mid, s_range[1]))) # Top-Left
    results.extend(get_adaptive_bins(subset, (m_mid, m_range[1]), (s_mid, s_range[1]))) # Top-Right
    
    return results

# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_adaptive(df, method_text = '_adaptive'):
    '''
    Bins in SFMS plane in an adaptive but regular-binned way
    Returns dataframe with additional columns containing '_linear', and list of unique bins
    '''
    final_bins = get_adaptive_bins(df, (log_mass_bins[0],log_mass_bins[-1]), (log_sfr_bins[0],log_sfr_bins[-1]), max_n=args.max_gal_per_bin)

    df[f'bin_id{method_text}'] = -1
    df[f'mass_interval{method_text}'] = None
    df[f'sfr_interval{method_text}'] = None
    
    for i, b in enumerate(final_bins):
        mask = (df['log_mass'] >= b['m_min']) & (df['log_mass'] < b['m_max']) & (df['log_sfr'] >= b['s_min']) & (df['log_sfr'] < b['s_max'])
        df.loc[mask, f'bin_id{method_text}'] = i
        
        df.loc[mask, f'mass_interval{method_text}'] = pd.Interval(left=b['m_min'], right=b['m_max'], closed='left')
        df.loc[mask, f'sfr_interval{method_text}'] = pd.Interval(left=b['s_min'], right=b['s_max'], closed='left')
        df[f'bin_intervals{method_text}'] = list(zip(df[f'mass_interval{method_text}'], df[f'sfr_interval{method_text}']))

    df = df[df[f'bin_id{method_text}'] != -1].copy()
    bin_list = pd.unique(df[f'bin_intervals{method_text}'])

    return df, bin_list

# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_linear(df, method_text = '_linear'):
    '''
    Bins in SFMS plane in a linear, regular-binned way
    Returns dataframe with additional columns containing '_linear', and list of unique bins
    '''
    df[f'mass_interval{method_text}'] = pd.cut(df['log_mass'], bins=log_mass_bins)
    df[f'sfr_interval{method_text}'] = pd.cut(df['log_sfr'], bins=log_sfr_bins)
    df = df.dropna(subset=[f'mass_interval{method_text}', f'sfr_interval{method_text}'])
    df[f'bin_intervals{method_text}'] = list(zip(df[f'mass_interval{method_text}'], df[f'sfr_interval{method_text}']))

    all_mass_intervals = df[f'mass_interval{method_text}'].cat.categories
    all_sfr_intervals = df[f'sfr_interval{method_text}'].cat.categories
    bin_list = list(itertools.product(all_mass_intervals, all_sfr_intervals))

    return df, bin_list

# --------------------------------------------------------------------------------------------------------------------
def bin_by_method(df, method):
    '''
    Decides which function to call for the binning in SFMS plane, based on the input method
    Returns dataframe with additional columns containing '_{method}', and list of unique bins
    '''
    if method == 'linear': df, _ = bin_SFMS_linear(df, method_text=f'_{method}')
    elif method == 'adaptive': df, _ = bin_SFMS_adaptive(df, method_text=f'_{method}')

    return df

# ----------declaring mass and SFR bins-------------------
delta_log_mass, delta_log_sfr = 1, 0.5
log_mass_bins = np.arange(7.5, 11.5 + delta_log_mass/2, delta_log_mass)
log_sfr_bins = np.arange(-0.5, 2.5 + delta_log_sfr/2, delta_log_sfr)
    
# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.5

    methods = ['linear', 'adaptive']

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
        if args.adaptive_bins: output_dir = Path(str(output_dir).replace('stacking', 'stacking_adaptive'))
        output_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir / f'plots'
    
        # ---------read the photometric catalog file--------------------
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        df['field'] = args.field

        # ---------crossmatch with cosmos-web to get stellar mass and SFR--------------------
        df = get_passage_masses_from_cosmos(df, args, id_col='id')

        # ---------merge with effective radius catalog--------------------
        df_re = Table.read(args.output_dir / f'catalogs/{args.field}_re_list.fits').to_pandas()
        if 'redshift' in df_re: df_re.drop(columns=['redshift'], axis=1, inplace=True)
        df_re['id'] = df_re['id'].astype(int)
        df_re = df_re[df_re['re_kpc'] > 0]
        df = pd.merge(df, df_re, on='id', how='inner')

        # -----------making bins in various ways---------------------
        for method in methods: df = bin_by_method(df, method)

        # ------------plotting stacked gradients on SFMS--------------------------
        fig = plot_SFMS_bins(df, methods, args)
        save_fig(fig, fig_dir, f'{args.field}_SFMS_binned.png', args) # saving the figure

        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
