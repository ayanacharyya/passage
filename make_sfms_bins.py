'''
    Filename: make_sfms_bins.py
    Notes: Bins the SFR-mass plane in various ways, for testing to decide which type of binning has the best optimisation between more bins and more galaxies per bin
    Author : Ayan
    Created: 30-01-26
    Example: run make_sfms_bins.py --field Par028 --overplot_literature --overplot_passage
             run make_sfms_bins.py --field Par028
             run make_sfms_bins.py --do_all_fields
'''

from header import *
from util import *
from make_passage_plots import plot_SFMS_Popesso23, plot_SFMS_Shivaei15, plot_SFMS_Whitaker14

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
            ax.text(cx, cy, int(row[quant]), color='black', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold')
    
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

    return ax

# --------------------------------------------------------------------------------------------------------------------
def make_heatmap_vorbin(ax, df, bin_summary, centers_scaled, scaling, quant, args, method_text='_voronoi', cmap='viridis', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, cmin=None, cmax=None, clabel='', ncbins=4):
    '''
    Makes heatmap from a dataframe that has been Voronoi binned, using patches, and annotates based on another pivot table, onto a given axis handle
    Returns axis handle
    '''    
    # ---------obtaining voronoi bin properties---------
    mean, std = scaling
    centers = (centers_scaled * std) + mean
    vor = Voronoi(centers)
    voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='k', line_width=1)

    # -----------defining vertices, for annotating and color-coding the bins-------------
    patches, values = [], []
        
    for index, row in bin_summary.iterrows():
        bin_id = int(row[f'bin_id{method_text}'])
        region_idx = vor.point_region[bin_id] # vor.point_region maps the index of the generator point to the index of the region
        region_vertices_indices = vor.regions[region_idx]
        verts = vor.vertices[region_vertices_indices] # Get the actual coordinates of the vertices
        
        #verts_phys = (verts * std) + mean
        patches.append(Polygon(verts, closed=True))
        values.append(row[quant])
        print(f'Deb110: index={index}, patch={Polygon(verts, closed=True)}, value={row[quant]}') ##
        
        ax.text(row['m_center'], row['s_center'], int(row[quant]), color='k', ha='center', va='center', fontsize=args.fontsize / args.fontfactor, fontweight='bold') # Annotate: Place text at the physical center of the bin

    p = PatchCollection(patches, cmap=cmap, edgecolors='w', alpha=0.3, linewidths=1)
    p.set_array(np.array(values))
    p.set_clim(vmin=cmin, vmax=cmax)
    print(f'Deb117: len={len(patches)}, patches={p}, values={values}') ##
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
def plot_SFMS_bins(df, methods, args, scaling=None, centers_scaled=None, bin_summary=None):
    '''
    Makes a nice heatmap (with patches) of stacked integrated metallicities and metallicity gradients
    Returns figure handle
    '''
    cmap, cmin, cmax, ncbins, clabel ='viridis', 10, 50, 4, 'Number of galaxies'
    ncol = len(methods)
    
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, ncol, figsize=(5.2 * ncol, 5.), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    fig.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.13, wspace=0., hspace=0.)
    
    # ---------plot the heatmaps-------------------
    for index, method in enumerate(methods):
        if 'vor' in method:
            axes[index] = make_heatmap_vorbin(axes[index], df, bin_summary, centers_scaled, scaling, 'n_galaxies', args, method_text='_voronoi', cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
        else:
            df_sub = df.groupby(f'bin_intervals_{method}').size().reset_index(name='n_galaxies')
            df_sub[['log_mass_bin', 'log_sfr_bin']] = pd.DataFrame(df_sub[f'bin_intervals_{method}'].tolist(), index=df_sub.index)
            axes[index] = make_heatmap_patches(axes[index], df_sub, 'n_galaxies', args, xcolname=f'log_mass_bin', ycolname=f'log_sfr_bin', cmap=cmap, hide_cbar=True, hide_yaxis=index, cmin=cmin, cmax=cmax)
            axes[index].text(0.05, 0.95, f'{method}', ha='left', va='top', c='k', fontsize=args.fontsize, transform=axes[index].transAxes)

        axes[index].set_aspect('equal')

        # ----------over-plotting data and theoretical diagrams----------
        if args.overplot_passage:
            axes[index].scatter(df['log_mass'], df['log_sfr'], s=5, c='w', lw=1, edgecolors='sienna', label=f'{args.field}') # overplot PASSAGE galaxies (integrated stellar mass-SFR)
            if index == 0: axes[index].legend(fontsize=args.fontsize / args.fontfactor, loc='lower right')

        if args.overplot_literature:
            #axes[index] = plot_SFMS_Whitaker14(axes[index], 2, color='yellowgreen')
            axes[index] = plot_SFMS_Shivaei15(axes[index], color='darkgreen')
            axes[index] = plot_SFMS_Popesso23(axes[index], 2, color='darkgoldenrod')
            #axes[index] = plot_SFMS_Popesso23(axes[index], 3, color='royalblue')
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
def bin_by_method(df, method):
    '''
    Decides which function to call for the binning in SFMS plane, based on the input method
    Returns dataframe with additional columns containing '_{method}', and list of unique bins
    '''
    if method == 'linear': output = bin_SFMS_linear(df, method_text=f'_{method}')
    elif method == 'adaptive_nmax': output= bin_SFMS_adaptive(df, method_text=f'_{method}', max_n=40)
    elif method == 'adaptive_nmin': output = bin_SFMS_adaptive(df, method_text=f'_{method}', min_n=20)
    elif 'vor' in method: output = bin_SFMS_voronoi(df, method_text=f'_{method}', target_n=30)

    return output

# --------------------------------------------------------------------------------------------------------------------
def read_passage_sed_catalog(filename):
    '''
    Read the combined master catalog from PASSAGE SED fits, rename a few columns, and only keep the mass and SFR columns
    Return pandas dataframe
    '''
    print(f'Reading master PASSAGE SED catalog from {filename}..')
    full_df = Table.read(filename).to_pandas()
    full_df.rename(columns={'Par':'field', 'passage_id':'id', 'stellar_mass_50':'log_mass', 'ssfr_50':'log_ssfr', 'sfr_50':'sfr', 'ra_obj':'ra', 'dec_obj':'dec', }, inplace=True)

    full_df['log_sfr'] = np.log10(full_df['sfr'])
    columns_to_extract = ['field', 'id', 'zbest', 'log_mass', 'log_sfr', 'log_ssfr', 'cosmosid']
    df = full_df[columns_to_extract]
    df['field'] = df['field'].astype(str)
    df.rename(columns={'zbest':'redshift'}, inplace=True)

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

    methods = ['linear', 'adaptive_nmin', 'voronoi']# 'adaptive_nmin']

    # ---------reading in the master SED catalog----------------
    if args.do_all_fields:
        passage_catalog_filename = args.output_dir / 'catalogs' / 'passagepipe_v0.5_SED_fits_cosmosweb_v1.0.0-alpha.fits'
        df = read_passage_sed_catalog(passage_catalog_filename)
        output_dir = args.output_dir / 'stacking'
    
    # ---------reading in the single-field phot catalog----------------
    else:
        product_dir = args.input_dir / args.field / 'Products'
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        df['field'] = args.field
        df = get_passage_masses_from_cosmos(df, args, id_col='id') # crossmatch with cosmos-web to get stellar mass and SFR
        output_dir = args.output_dir / args.field / 'stacking'

    # ------determining field-specific paths, etc-----------
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / f'plots'

    # -----------making bins in various ways---------------------
    for method in methods:
        output = bin_by_method(df, method)
        if 'vor' in method:
            df, bin_summary, centers_scaled, scaling = output
        else:
            df, bin_list = output
            bin_summary, centers_scaled, scaling = None, None, None

    # ------------plotting stacked gradients on SFMS--------------------------
    fig = plot_SFMS_bins(df, methods, args, centers_scaled=centers_scaled, scaling=scaling, bin_summary=bin_summary)
    save_fig(fig, fig_dir, f'SFMS_binned.png', args) # saving the figure

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
