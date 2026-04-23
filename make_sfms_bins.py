'''
    Filename: make_sfms_bins.py
    Notes: Bins the SFR-mass plane in various ways, for testing to decide which type of binning has the best optimisation between more bins and more galaxies per bin
    Author : Ayan
    Created: 30-01-26
    Example: run make_sfms_bins.py --field Par028 --overplot_literature --overplot_passage
             run make_sfms_bins.py --field Par028
             run make_sfms_bins.py --system ssd --do_all_fields
'''

from header import *
from util import *
from make_passage_plots import plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14, get_SFMS_Popesso23, get_SFMS_Shivaei15, get_SFMS_Whitaker14

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_patches(ax, df, quant, args, xcolname='log_mass_bin', ycolname='log_sfr_bin', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
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
    p.set_clim(vmin=cmin, vmax=cmax)
    ax.add_collection(p)

    # --------annotate at the center of each patch------------
    for i, row in df.iterrows():
        if row[quant] > 0:
            cx = (row['m_min'] + row['m_max']) / 2
            cy = (row['s_min'] + row['s_max']) / 2
            ax.text(cx, cy, int(row[quant]), color='k' if int(row[quant]) > 30 else 'w', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold')
    
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
            #spine.set_edgecolor('black')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_vorbin(ax, df, bin_summary, centers_scaled, scaling, quant, args, method_text='_voronoi', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe that has been Voronoi binned, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
    # ---------obtaining voronoi bin properties---------
    mean, std = scaling
    #centers = (centers_scaled * std) + mean
    #vor_orig_units = Voronoi(centers)
    #voronoi_plot_2d(vor_orig_units, ax=ax, show_vertices=False, show_points=False, line_colors='k', line_width=1)

    # --------adding 4 points very far away in scaled space, to handle infinite patches-------
    buffer = 100
    dummy_points = np.array([[buffer, buffer], [buffer, -buffer], [-buffer, buffer], [-buffer, -buffer]])
    points_for_vor = np.vstack([centers_scaled, dummy_points])
    vor = Voronoi(points_for_vor)

    # -----------defining vertices, for annotating and color-coding the bins-------------
    patches, values = [], []

    for index, row in bin_summary.iterrows():
        bin_id = int(row[f'bin_id{method_text}'])
        region_idx = vor.point_region[bin_id] # vor.point_region maps the index of the generator point to the index of the region
        vertex_indices = vor.regions[region_idx]
        
        if -1 not in vertex_indices and len(vertex_indices) > 0:
            verts_scaled = vor.vertices[vertex_indices]
            verts = (verts_scaled * std) + mean # Rescale vertices back to physical units
            
            patches.append(Polygon(verts, closed=True))
            values.append(bin_summary.iloc[index][quant])
            ax.text(row['m_center'], row['s_center'], int(row[quant]), color='k' if int(row[quant]) > 30 else 'w', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold') # Annotate: Place text at the physical center of the bin
        else:
            pass    
    
    p = PatchCollection(patches, cmap=cmap, edgecolors='w', alpha=0.8)
    p.set_array(np.array(values))
    p.set_clim(vmin=cmin, vmax=cmax)
    ax.add_collection(p)
 
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
def make_heatmap_distance(ax, df, sfms, quant, args, method_text='_distance', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
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
        ax.fill_between(m_grid, sfms_line + interval.left, sfms_line + interval.right, color=color, alpha=0.8, edgecolor='k', lw=0.5)
        
        s_center = sfms_func(row['log_mass_median']) + (interval.left + interval.right) / 2
        ax.text(row['log_mass_median'], s_center, int(row[quant]), color='k', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold', rotation=45)
 
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
def make_scattermap_sfh(ax, df, args, method_text='_sfh', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe that has been binned by distance from the SFMS, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''
    bin_list = list(pd.unique(df[f'bin_intervals{method_text}']))
    sorted_bins = sorted(bin_list)
    midpoints = [b.mid for b in sorted_bins]
    norm = mplcolors.Normalize(vmin=min(midpoints), vmax=max(midpoints))
    cmap = plt.get_cmap(cmap)
    bin_colors = [cmap(norm(m)) for m in midpoints]
    
    for index2, this_bin in enumerate(bin_list):
        df_sub = df[df[f'bin_intervals{method_text}'] == this_bin]
        ax.scatter(df_sub['log_mass'], df_sub['log_sfr'], s=10, c=bin_colors[index2], lw=0.1, edgecolors='k')
 
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
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.02)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ax=ax, label=clabel, cax=cax)
        cbar.set_label(clabel, fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)
        cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)
        cbar.update_ticks()

        cbar.ax.yaxis.set_ticks_position('left')
        cbar.ax.yaxis.set_label_position('left')
    
    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_bins(df, methods, method_texts, args, scaling=None, centers_scaled=None, bin_summary=None, sfms='Popesso23'):
    '''
    Makes a nice heatmap (with patches) of stacked integrated metallicities and metallicity gradients
    Returns figure handle
    '''
    cmap, cmin, cmax, ncbins, clabel ='viridis', 10, 50, 4, 'Number of galaxies'
    ncol = len(methods)
    
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, ncol, figsize=(5.2 * ncol, 5.), sharex=True, sharey=True, layout='constrained')
    axes = np.atleast_1d(axes)
    
    # ---------plot the heatmaps-------------------
    for index, method in enumerate(methods):
        if 'voronoi' in method:
            axes[index] = make_heatmap_vorbin(axes[index], df, bin_summary, centers_scaled, scaling, 'n_galaxies', args, method_text=method_texts[index], cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
        elif 'sfh' in method:
            axes[index] = make_scattermap_sfh(axes[index], df, args, method_text=method_texts[index], cmap=cmap, hide_cbar=False, hide_yaxis=index)
        elif 'distance_mass' in method:
            groupby_cols = [f'bin_intervals{method_texts[index]}', f'mass_intervals{method_texts[index]}']
            df_sub = df.groupby(groupby_cols).agg(n_galaxies=('log_mass', 'size'), log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median')).reset_index()
            df_sub = df_sub.dropna(subset=groupby_cols)
            df_sub = df_sub[df_sub['n_galaxies'] > 0]
            axes[index] = make_heatmap_distance(axes[index], df_sub, sfms, 'n_galaxies', args, method_text=method_texts[index], cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
        elif 'distance' in method:
            df_sub = df.groupby(f'bin_intervals{method_texts[index]}').agg(n_galaxies=('log_mass', 'size'), log_mass_min=('log_mass', 'min'), log_mass_max=('log_mass', 'max'), log_mass_median=('log_mass', 'median')).reset_index()
            df_sub = df_sub.dropna(subset=[f'bin_intervals{method_texts[index]}'])
            df_sub = df_sub[df_sub['n_galaxies'] > 0]
            axes[index] = make_heatmap_distance(axes[index], df_sub, sfms, 'n_galaxies', args, method_text=method_texts[index], cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
        else:
            df_sub = df.groupby(f'bin_intervals{method_texts[index]}').size().reset_index(name='n_galaxies')
            df_sub[['log_mass_bin', 'log_sfr_bin']] = pd.DataFrame(df_sub[f'bin_intervals{method_texts[index]}'].tolist(), index=df_sub.index)
            axes[index] = make_heatmap_patches(axes[index], df_sub, 'n_galaxies', args, xcolname=f'log_mass_bin', ycolname=f'log_sfr_bin', cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
        axes[index].text(0.05, 0.95, f'{method}', ha='left', va='top', c='k', fontsize=args.fontsize, transform=axes[index].transAxes)
        axes[index].set_aspect('equal')

        # ----------over-plotting data and theoretical diagrams----------
        if args.overplot_passage:
            axes[index].scatter(df['log_mass'], df['log_sfr'], s=5, c='w', lw=1, edgecolors='sienna', label='Huberty+26' if args.do_all_fields else f'{args.field}') # overplot PASSAGE galaxies (integrated stellar mass-SFR)
            if index == 0: axes[index].legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

        if args.overplot_literature:
            if sfms == 'Whitaker14': axes[index] = plot_SFMS_Whitaker14(axes[index], 2, color='crimson')
            elif sfms == 'Shivaei15': axes[index] = plot_SFMS_Shivaei15(axes[index], color='royalblue')
            elif sfms == 'Popesso23': axes[index] = plot_SFMS_Popesso23(axes[index], 2, color='darkgoldenrod')
            #elif sfms == 'Popesso23': axes[index] = plot_SFMS_Popesso23(axes[index], 3, color='royalblue')
            else: raise ValueError(f'Method {sfms} not found in the list of available SFMS literature methods: Whitaker14, Shivaei15, Popesso23')
            if index == 0: axes[index].legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

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
def get_adaptive_bins_nmax(df_subset, m_range, s_range, max_n=20):
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
    results.extend(get_adaptive_bins_nmax(subset, (m_range[0], m_mid), (s_range[0], s_mid))) # Bottom-Left
    results.extend(get_adaptive_bins_nmax(subset, (m_mid, m_range[1]), (s_range[0], s_mid))) # Bottom-Right
    results.extend(get_adaptive_bins_nmax(subset, (m_range[0], m_mid), (s_mid, s_range[1]))) # Top-Left
    results.extend(get_adaptive_bins_nmax(subset, (m_mid, m_range[1]), (s_mid, s_range[1]))) # Top-Right
    
    return results

# --------------------------------------------------------------------------------------------------------------------
def get_adaptive_bins_nmin(df, m_range, s_range, min_n=20):
    '''
    Recursively bins the Mass-SFR plane ensuring each bin has at least min_n galaxies.
    Courtesy of this function: Gemini
    '''    
    m_min, m_max = m_range
    s_min, s_max = s_range
    
    # Filter galaxies within this current rectangle
    mask = (df['log_mass'] >= m_min) & (df['log_mass'] < m_max) & (df['log_sfr'] >= s_min) & (df['log_sfr'] < s_max)
    subset = df[mask]
    count = len(subset)
    
    # If we are already below min_n, this is a terminal bin
    if count < min_n:
        return []

    # Calculate midpoints for a potential 2x2 split
    m_mid = (m_min + m_max) / 2
    s_mid = (s_min + s_max) / 2
    
    # Define the 4 potential quadrants
    quads = [((m_min, m_mid), (s_min, s_mid)), ((m_mid, m_max), (s_min, s_mid)), ((m_min, m_mid), (s_mid, s_max)), ((m_mid, m_max), (s_mid, s_max))]
    
    # Check if EVERY quadrant in the potential split would have >= min_n
    can_split = True
    for (mq, sq) in quads:
        q_count = df[(df['log_mass'] >= mq[0]) & (df['log_mass'] < mq[1]) & (df['log_sfr'] >= sq[0]) & (df['log_sfr'] < sq[1])].shape[0]
        if q_count < min_n:
            can_split = False
            break
            
    # If we can split, recurse into children
    if can_split:
        bins = []
        for (mq, sq) in quads:
            bins.extend(get_adaptive_bins_nmin(df, mq, sq, min_n))
        return bins
    else:
        # Cannot split further without violating min_n, return current bin
        return [{'m_min': m_min, 'm_max': m_max, 's_min': s_min, 's_max': s_max, 'n_galaxies': count, 'galaxy_ids': subset['id'].tolist()}]
    
# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_adaptive(df, method_text = '_adaptive', max_n=None, min_n=None):
    '''
    Bins in SFMS plane in an adaptive but regular-binned way
    Returns dataframe with additional columns containing '_linear', and list of unique bins
    '''
    if max_n is not None: final_bins = get_adaptive_bins_nmax(df, (log_mass_bins[0],log_mass_bins[-1]), (log_sfr_bins[0],log_sfr_bins[-1]), max_n=max_n)
    elif min_n is not None: final_bins = get_adaptive_bins_nmin(df, (log_mass_bins[0],log_mass_bins[-1]), (log_sfr_bins[0],log_sfr_bins[-1]), min_n=min_n)

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
    bin_list = list(pd.unique(df[f'bin_intervals{method_text}']))

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
def bin_SFMS_voronoi(df, method_text='_voronoi', target_n=20):
    '''
    Bins the Mass-SFR plane ensuring each bin has approximately target_n galaxies.
    Courtesy of this function: Gemini
    '''    
    n_bins = len(df) // target_n
    coords = df[['log_mass', 'log_sfr']].values
    coords_mean = coords.mean(axis=0)
    coords_std = coords.std(axis=0)
    coords_scaled = (coords - coords_mean) / coords_std # we should normalize so one doesn't dominate the distance calculation

    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    #kmeans = KMeans(n_clusters=n_bins, n_init=10, random_state=42) # This finds centers that partition the space into target_n chunks
    kmeans = KMeans(n_clusters=n_bins, n_init='auto', random_state=42)
    df[f'bin_id{method_text}'] = kmeans.fit_predict(coords_scaled)
    
    bin_summary = df.groupby(f'bin_id{method_text}').agg(n_galaxies=('id', 'count'), m_center=('log_mass', 'mean'), s_center=('log_sfr', 'mean'), galaxy_ids=('id', list)).reset_index()
    
    return df, bin_summary, kmeans.cluster_centers_, (coords_mean, coords_std)

# --------------------------------------------------------------------------------------------------------------------
def get_sfms_func(log_mass_arr, method):
    '''
    Compute the log_sfr_arr based on an input log_mass_arr and SFMS method
    Interpolates between the input log_mass_arr and the derived log_sfr_arr
    Returns the interpolated function
    '''
    log_mass_min, log_mass_max = np.min(log_mass_arr), np.max(log_mass_arr)
    redshift = 1.8 # based on the median PASSAGE redshift from Huvery+26

    if method == 'Whitaker14': (log_mass1, log_SFR1), (log_mass2, log_SFR2), (z1, z2) = get_SFMS_Whitaker14(log_mass_min, log_mass_max, redshift)
    elif method == 'Shivaei15': (log_mass1, log_SFR1), (log_mass2, log_SFR2), scatter = get_SFMS_Shivaei15(log_mass_min, log_mass_max)
    elif method == 'Popesso23': (log_mass1, log_SFR1), (log_mass2, log_SFR2) = get_SFMS_Popesso23(log_mass_min, log_mass_max, redshift)
    else: raise ValueError(f'Method {method} not found in the list of available SFMS literature methods: Whitaker14, Shivaei15, Popesso23')

    log_mass_arr = np.hstack([log_mass1, log_mass2])
    log_sfr_arr = np.hstack([log_SFR1, log_SFR2])
    if type(log_sfr_arr[0]) != np.float64: log_sfr_arr = unp.nominal_values(log_sfr_arr) # if it is an uncertainties array then just take the nominal values

    sfms_func = interp1d(log_mass_arr, log_sfr_arr, kind='linear', fill_value='extrapolate')

    return sfms_func

# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_distance(df, method_text = '_distance', delta_bin=0.2, n_adaptive_bins=None, sfms='Popesso23'):
    '''
    Bins in SFMS plane based on the distance from a given SFMS relation in bins of width delta_bin
    Returns dataframe with additional columns containing '_distance', and list of bin edges
    '''
    sfms_func = get_sfms_func(df['log_mass'], sfms)

    df['log_sfr_ms_expected'] = sfms_func(df['log_mass'])
    df['delta_sfms'] = df['log_sfr'] - df['log_sfr_ms_expected']
    df = df.dropna(subset=['delta_sfms'])

    if n_adaptive_bins is None: # make uniform width bins of width delta_bin
        max_dist = np.ceil(df['delta_sfms'].abs().max() / delta_bin) * delta_bin
        bin_edges = np.arange(-max_dist, max_dist + delta_bin, delta_bin)
        df[f'bin_intervals{method_text}'] = pd.cut(df['delta_sfms'], bins=bin_edges)
    else: # make n_adaptive_bins of non-uniform width
        df[f'bin_intervals{method_text}'], bin_edges = pd.qcut(df['delta_sfms'], q=n_adaptive_bins, retbins=True)

    bin_list = pd.unique(df[f'bin_intervals{method_text}'])

    return df, bin_list

# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_distance_mass(df, method_text = '_distance_mass', delta_bin=0.2, n_adaptive_bins=None, n_mass_bins=4, sfms='Popesso23'):
    '''
    Bins in SFMS plane based on the distance from a given SFMS relation in bins of width delta_bin
    Returns dataframe with additional columns containing '_distance', and list of bin edges
    '''
    sfms_func = get_sfms_func(df['log_mass'], sfms)

    df['log_sfr_ms_expected'] = sfms_func(df['log_mass'])
    df['delta_sfms'] = df['log_sfr'] - df['log_sfr_ms_expected']
    df = df.dropna(subset=['delta_sfms'])

    if n_adaptive_bins is None: # make uniform width bins of width delta_bin
        max_dist = np.ceil(df['delta_sfms'].abs().max() / delta_bin) * delta_bin
        bin_edges = np.arange(-max_dist, max_dist + delta_bin, delta_bin)
        df[f'bin_intervals{method_text}'] = pd.cut(df['delta_sfms'], bins=bin_edges)
        df[f'mass_intervals{method_text}'] = pd.cut(df['log_mass'], bins=n_mass_bins) # make n_mass_bins of equal widths
    else: # make n_adaptive_bins of non-uniform width
        df[f'bin_intervals{method_text}'], bin_edges = pd.qcut(df['delta_sfms'], q=n_adaptive_bins, retbins=True)
        df[f'mass_intervals{method_text}'] = df.groupby(f'bin_intervals{method_text}')['log_mass'].apply(lambda x: pd.qcut(x, q=n_mass_bins)).reset_index(level=0, drop=True) # make n_mass_bins of widths decided to have roughly equal galaxies in each bin
    
    unique_bins_df = df[[f'bin_intervals{method_text}', f'mass_intervals{method_text}']].drop_duplicates().sort_values([f'bin_intervals{method_text}', f'mass_intervals{method_text}'])
    bin_list = list(unique_bins_df.to_records(index=False))

    return df, bin_list

# --------------------------------------------------------------------------------------------------------------------
def bin_SFMS_sfh(df, method_text = '_sfh', delta_sfh=0.2, n_sfh_bins=None, bin_by_col='delta_tform_90_10'):
    '''
    Bins in SFMS plane in an adaptive (or, optionally, regular) way based on a SFH parameter
    Returns dataframe with additional columns containing method_text, and list of unique bins
    '''

    if n_sfh_bins is not None:
        df[f'bin_intervals{method_text}'] = pd.qcut(df[bin_by_col], q=n_sfh_bins)
    else:
        custom_edges = np.arange(np.round(df[bin_by_col].min(), 2), np.round(df[bin_by_col].max(), 2) + delta_sfh, delta_sfh)
        df[f'bin_intervals{method_text}'] = pd.cut(df[bin_by_col], bins=custom_edges)

    bin_list = list(pd.unique(df[f'bin_intervals{method_text}']))

    return df, bin_list

# --------------------------------------------------------------------------------------------------------------------
def read_passage_sed_catalog(filename, use_old=True):
    '''
    Read the combined master catalog from PASSAGE SED fits, rename a few columns, and only keep the mass and SFR columns
    Return pandas dataframe
    '''
    if use_old:
        print(f'Reading master PASSAGE SED catalog from {filename}..')
        full_df = Table.read(filename).to_pandas()
    else:
        print(f'Reading PASSAGE line finding catalog from {filename}..')
        full_df = pd.read_csv(filename, header=0, sep='\t')
        full_df = full_df[full_df['stellar_mass_50'] > 0].reset_index(drop=True) # to get only those sources that have stellar mass measured

    full_df.rename(columns={'Par':'field', 'passage_id':'id', 'id_photcat':'id', 'objid':'id', 'cosmoswebid_1':'cosmosid', 'stellar_mass_50':'log_mass', 'ssfr_50':'log_ssfr', 'sfr_50':'sfr', 'ra_obj':'ra', 'dec_obj':'dec'}, inplace=True)
    full_df['log_sfr'] = np.log10(full_df['sfr'])
    
    full_df['delta_tform_50_10'] = full_df['tform50_50'] - full_df['tform10_50']
    full_df['delta_tform_90_50'] = full_df['tform90_50'] - full_df['tform50_50']
    full_df['delta_tform_90_10'] = full_df['tform90_50'] - full_df['tform10_50']

    columns_to_extract = ['field', 'id', 'zbest', 'log_mass', 'log_sfr', 'log_ssfr', 'cosmosid', 'delta_tform_50_10', 'delta_tform_90_50', 'delta_tform_90_10']
    df = full_df[columns_to_extract]
    df['field'] = df['field'].astype(str)
    df.rename(columns={'zbest':'redshift'}, inplace=True)

    df = df.dropna(subset=['log_mass', 'log_sfr'])

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_binned_df(args, skip_binning=False, df=None, method_text='', skip_stacking=False):
    '''
    Read in the full SED dataframe and curtail and bin it as needed to produce the binned dataframe
    Returns binned pandas dataframe, bin_list, args (with added keywords)
    '''
    # ------------add new keywords to args---------------
    # args.line_list = ['OIII', 'OII', 'NeIII-3867', 'Hb', 'OIII-4363', 'Ha', 'SII']
    args.line_list = ['OIII', 'OII', 'NeIII-3867', 'Hb', 'Ha', 'SII']

    args.deproject_text = '_nodeproject' if args.skip_deproject else ''
    args.rescale_text = '_norescale' if args.skip_re_scaling else ''
    args.C25_text = '_wC25' if args.use_C25 and 'NB' not in args.Zdiag else ''
    args.fold_text = '_folded' if args.fold_maps else ''

    # ---------reading in the master SED catalog----------------
    if df is None:
        passage_catalog_filename = args.output_dir / 'catalogs' / 'SED_fits_v1.0.2_cosmosweb.fits'
        df = read_passage_sed_catalog(passage_catalog_filename)

    if args.do_all_fields:
        re_catalog_filename = args.output_dir / f'catalogs/all_fields_re_list.fits'
        root_dir = args.output_dir
    # ---------curtailing to single-field----------------
    else:
        df = df[df['field'] == args.field]        
        re_catalog_filename = args.output_dir / f'catalogs/{args.field}_re_list.fits'
        root_dir = args.output_dir / args.field

    # --------------returning unbinned dataframe and no bin_list if skip_binning=True--------------
    if skip_binning:
        bin_list = None        
        args.fig_dir = root_dir / 'stacking'
        args.fig_dir.mkdir(parents=True, exist_ok=True)
    else:
        # -------------binning the mass-SFR plane-------------
        if args.adaptive_bins:
            if args.voronoi_bins:
                df, bin_summary, centers_scaled, scaling = bin_SFMS_voronoi(df, method_text=method_text, target_n=target_n)
                bin_list = [bin_summary, centers_scaled, scaling]
                args.binby_text = f'adap_binby_voronoi_ngal_{target_n}'
            elif args.bin_by_sfh:
                df, bin_list = bin_SFMS_sfh(df, method_text=method_text, n_sfh_bins=n_sfh_bins) # binning the dataframe in an adaptive way
                args.binby_text = f'adap_binby_sfh_{n_sfh_bins}'
            elif args.bin_by_distance:
                df, bin_list = bin_SFMS_distance(df, method_text=method_text, n_adaptive_bins=n_adaptive_bins, sfms=sfms)
                args.binby_text = f'adap_binby_distance_{n_adaptive_bins}_sfms_{sfms}'
            elif args.bin_by_distance_mass:
                df, bin_list = bin_SFMS_distance_mass(df, method_text=method_text, n_adaptive_bins=n_adaptive_bins, sfms=sfms, n_mass_bins=n_mass_bins)
                args.binby_text = f'adap_binby_distance_{n_adaptive_bins}_sfms_{sfms}_mass_{n_mass_bins}'
            else:
                if args.min_gal_per_bin is not None:
                    df, bin_list = bin_SFMS_adaptive(df, method_text=method_text, min_n=args.min_gal_per_bin) # binning the dataframe in an adaptive way
                    args.binby_text = f'adap_binby_mass_sfr_mingal_{args.min_gal_per_bin}'
                else:
                    df, bin_list = bin_SFMS_adaptive(df, method_text=method_text, max_n=args.max_gal_per_bin) # binning the dataframe in an adaptive way
                    args.binby_text = f'adap_binby_mass_sfr_maxgal_{args.max_gal_per_bin}'
        else:
            if args.bin_by_sfh:
                df, bin_list = bin_SFMS_sfh(df, method_text=method_text, delta_sfh=delta_sfh) # binning the dataframe in an adaptive way
                args.binby_text = f'lin_binby_sfh_delta_sfh_{delta_sfh}'
            elif args.bin_by_distance:
                df, bin_list = bin_SFMS_distance(df, method_text=method_text, delta_bin=delta_sfms, sfms=sfms)
                args.binby_text = f'lin_binby_distance_delta_bin_{delta_sfms}_sfms_{sfms}'
            elif args.bin_by_distance_mass:
                df, bin_list = bin_SFMS_distance_mass(df, method_text=method_text, delta_bin=delta_sfms, sfms=sfms, n_mass_bins=n_mass_bins)
                args.binby_text = f'lin_binby_distance_delta_bin_{delta_sfms}_sfms_{sfms}_mass_{n_mass_bins}'
            else:
                df, bin_list = bin_SFMS_linear(df, method_text=method_text) # -binning the dataframe uniformly by mass and SFR bins
                args.binby_text = f'lin_binby_mass_sfr_delta_mass_{delta_log_mass}_delta_sfr_{delta_log_sfr}'

        if not args.voronoi_bins:
            if args.bin_by_distance or args.bin_by_distance_mass or args.bin_by_sfh:
                bin_list = np.sort(bin_list)
            else:
                bin_list.sort(key=lambda x: (x[0].left, x[1].left))

        # ------determining field-specific paths, etc-----------
        output_dir = root_dir / 'stacking' / f'{args.binby_text}'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        args.fig_dir = output_dir / f'plots{args.deproject_text}{args.rescale_text}'
        Path(args.fig_dir / 'binmembers').mkdir(parents=True, exist_ok=True)
        
        args.fits_dir = output_dir / f'maps{args.deproject_text}{args.rescale_text}'
        args.fits_dir.mkdir(parents=True, exist_ok=True)

        args.grad_filename = args.fits_dir / f'stacked{args.binby_text}{args.fold_text}_fits_allbins_Zdiag_{args.Zdiag}{args.C25_text}{args.deproject_text}{args.rescale_text}.fits'

        # ---------merge with effective radius catalog--------------------
        if not skip_stacking:
            if os.path.exists(re_catalog_filename):
                df_re = Table.read(re_catalog_filename).to_pandas()
                if 'redshift' in df_re: df_re.drop(columns=['redshift'], axis=1, inplace=True)
                df_re['id'] = df_re['id'].astype(int)
                df_re['field'] = df_re['field'].astype(str)
                df_re = df_re[df_re['re_kpc'] > 0]
                df = pd.merge(df, df_re, on=['field', 'id'], how='inner')
            else:
                raise FileNotFoundError(f're catalog not found at {re_catalog_filename}; please create this file first by running compute_re.py. Exiting.')

            # -------merge with all photcats by looping over fields--------------
            print(f'Combining all photcats to one file before merging to df, might take a few seconds..')
            fields = pd.unique(df['field'])
            df_master_photcat = pd.DataFrame()
            cols_to_extract = ['field', 'id', 'a_image', 'b_image', 'theta_image']
            for this_field in fields:
                df_phot_this_field = GTable.read(args.input_dir / this_field / 'Products' / f'{this_field}_photcat.fits').to_pandas()
                df_phot_this_field['field'] = this_field
                df_master_photcat = pd.concat([df_master_photcat, df_phot_this_field[cols_to_extract]])
            
            df = pd.merge(df, df_master_photcat, on=['field', 'id'], how='inner')

        # --------------curtailiug bins for debugging-------------------
        if args.debug_bin: bin_list = bin_list[1:2]
        #if args.debug_bin: bin_list = [item for item in bin_list if (item[0].left == 9.5) & (item[0].right == 10.) & (item[1].left == 2.0) & (item[1].right == 2.5)] # to choose the mass=9.5-10.5, sfr=2-2.5 bin for debugging purposes

        # ------determining field-specific paths, etc-----------
        if args.skip_re_scaling:
            args.pix_size = (2 * args.kpc_limit) / args.npix_side # kpc
        else:
            args.pix_size = (2 * args.re_limit) / args.npix_side # Re
        
        if args.fold_maps:
            if args.skip_re_scaling:
                args.extent = (0, args.re_limit, 0, args.re_limit)
            else:
                args.extent = (0, args.re_limit, 0, args.re_limit)
        else: 
            offset = args.pix_size / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
            if args.skip_re_scaling:
                args.extent = (-args.kpc_limit - offset, args.kpc_limit - offset, -args.kpc_limit - offset, args.kpc_limit - offset)
            else:
                args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)

    return df, bin_list, args

# ----------declaring mass and SFR bins-------------------
delta_log_mass, delta_log_sfr = 0.5, 0.5
log_mass_bins = np.arange(6.5, 11.5 + delta_log_mass/2, delta_log_mass)
log_sfr_bins = np.arange(-2.5, 2.5 + delta_log_sfr/2, delta_log_sfr)

# ---------------choose which binning methods to try---------------
methods = [#'adaptive_nmin', \
            # 'adaptive_nmax', \
            # 'adaptive_voronoi', \
            'adaptive_distance', \
            'adaptive_distance_mass', \
            # 'adaptive_sfh', \
            # 'linear', \
            # 'linear_distance', \
            # 'linear_distance_mass', \
            # 'linear_sfh', \
            ]

target_n = 30 # for voronoi binning
n_adaptive_bins = 8 # for distance (from SFMS) binning
n_mass_bins = 4 # number of mass bins within each distance (from SFMS) bin
delta_sfms = 0.4 # delta in distance from SFMS in which to bin in the distance-from-SFMS method, unless binning adaptively
sfms =  'Popesso23' # from 'Popesso23', 'Shivaei15' and 'Whitaker14'; for binning by distance from SFMS
delta_sfh = 0.5 # delta in tform90 - tform10 parameter (in Gyr) in which to bin, unless binning adaptively
n_sfh_bins = 10 # number of bins in SFH parameter, adaptive, so that each bin contains approximately equal number of objects

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.5
    if args.re_limit is None: args.re_limit = 2.

    # ---------reading in the master SED catalog----------------
    passage_catalog_filename = args.output_dir / 'catalogs' / 'SED_fits_v1.0.2_cosmosweb.fits'
    df = read_passage_sed_catalog(passage_catalog_filename)

    # -----------making bins in various ways---------------------
    for method in methods:
        print(f'\nDoing method {method}')
        args.adaptive_bins, args.voronoi_bins, args.bin_by_sfh, args.bin_by_distance, args.bin_by_distance_mass = [False] * 5
        args.min_gal_per_bin = None

        if method.startswith('adaptive'): args.adaptive_bins = True
        if 'voronoi' in method: args.voronoi_bins = True
        elif 'sfh' in method: args.bin_by_sfh = True
        elif 'distance_mass' in method: args.bin_by_distance_mass = True
        elif 'distance' in method: args.bin_by_distance = True
        if '_nmin' in method and args.min_gal_per_bin is None: args.min_gal_per_bin = 20

        df, bin_list, args = get_binned_df(args, df=df, method_text=f'_{method}', skip_stacking=True)

        if 'voronoi' in method:
            bin_summary, centers_scaled, scaling = bin_list
        else:
            bin_summary, centers_scaled, scaling = [None] * 3

        # ------------plotting stacked gradients on SFMS--------------------------
        fig = plot_SFMS_bins(df, [method], [f'_{method}'], args, centers_scaled=centers_scaled, scaling=scaling, bin_summary=bin_summary, sfms=sfms)
        save_fig(fig, args.fig_dir, f'SFMS_binned.png', args) # saving the figure

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
