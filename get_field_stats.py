'''
    Filename: get_field_stats.py
    Notes: Computes and plots various statistics (redshift distribution, detected line distribution, magnitude distribution etc.) for a given PASSAGE field which has been already run through make_diagnostic_maps.py
    Author : Ayan
    Created: 19-08-24
    Example: run get_field_stats.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --re_extract
             run get_field_stats.py --field Par61 --mag_lim 26 --line_list OII,OIII,Ha
             run get_field_stats.py --mag_lim 26 --line_list OIII --do_all_fields --clobber
             run get_field_stats.py --mag_lim 26 --line_list OIII,Ha --do_all_fields --zmin 1 --zmax 2.5 --merge_visual --plot_conditions EW,z,mag,tail,RQ,strong_OIII,PA
             run get_field_stats.py --mag_lim 26 --line_list OIII,Ha --do_all_fields --plot_conditions EW,mag,compact
             run get_field_stats.py --mag_lim 24 --EW_thresh 300 --log_SFR_thresh 0 --line_list OIII,Ha --do_all_fields --plot_conditions EW,mag,PA,mass
             run get_field_stats.py --line_list OIII,Ha --do_all_fields --plot_conditions EW,mass,PA,a_image --clobber_venn
             run get_field_stats.py --line_list OIII,Ha --do_all_fields --plot_conditions EW,mass,PA --clobber_venn
             run get_field_stats.py --line_list Ha --do_all_fields --plot_conditions has,PA --clobber_venn
             run get_field_stats.py --line_list Ha --do_all_fields --plot_conditions SNR,mass,F115W,F150W,F200W --SNR_thresh 10 --clobber_venn
             run get_field_stats.py --line_list Ha --do_all_fields --plot_pie
             run get_field_stats.py --line_list Ha --do_all_fields --plot_sunburst
             run get_field_stats.py --line_list Ha --do_all_fields --plot_conditions mass --plot_columns --xcol redshift --ycol lp_zBEST --colorcol lp_mass
             run get_field_stats.py --line_list Ha --do_all_fields --plot_conditions mass --plot_columns --xcol redshift --ycol ez_z_phot --colorcol ez_mass
'''
from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# -------------------------------------------------------------------------------------------------------
def make_set(df, condition, label, set_arr, label_arr):
    '''
    Applies the given condition on given df and appends the ID list into a set and assigns a label
    '''
    id_list = df[condition]['par_obj'].values
    print(f'{len(id_list)} objects meet the {label} condition')
    set_arr.append(set(id_list))
    label_arr.append(label)

    return set_arr, label_arr

# -----------------------------------------------------------------------------------------------
def is_line_within_filter(line, redshift, filters=None, field=None, trim_factor=0.):
    '''
    Determines whether a given line is within the "good" part of a given list of filters' wavelength range
    Returns boolean
    '''
    if filters is None: filters = available_filters_for_field_dict[field]

    line_obs_wave = rest_wave_dict[line] * (1 + redshift) / 1e3 # factor 1e3 to convert nm to microns
    all_wave_ranges = [filter_waverange_dict[item] for item in filters] # "good" part of wavelength ranges for each filter, in microns
    trimmed_wave_ranges = [[item[0] * (1 + trim_factor), item[1] * (1 - trim_factor)] for item in all_wave_ranges]
    is_line_available = np.array([item[0] < line_obs_wave and item[1] > line_obs_wave for item in trimmed_wave_ranges]).any()

    return is_line_available

# -------------------------------------------------------------------------------------------------------
def plot_venn(df, args):
    '''
    To plot Venn diagrams with a given df, for a bunch of criteria
    Plots and saves the figure
    Returns intersecting dataframe
    '''
    n_fields = len(pd.unique(df["field"]))
    if 'par_obj' not in df: df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype(str)
    print(f'\nOut of the total {len(df)} objects in {n_fields} fields..\n')

    set_arr = []
    label_arr = []

    # ---------add line sets------------
    line_list = args.line_list

    for line in line_list:
        condition1 = df.apply(lambda x: is_line_within_filter(line, x['redshift'], field=x['field'], trim_factor=args.trim_filter_by_wavelength_factor), axis=1)
        set_arr, label_arr = make_set(df, condition1, f'{line} available', set_arr, label_arr)

        if 'present' in args.plot_conditions or 'SNR' in args.plot_conditions or 'EW' in args.plot_conditions:
            try:
                condition2 = (np.isfinite(df[f'{line}_EW'])) & (df[f'{line}_EW'] > 0)
                set_arr, label_arr = make_set(df, condition1 & condition2, f'{line} present', set_arr, label_arr)
            except:
                print(f'Could not apply the {line} present condition')
                pass
        if 'SNR' in args.plot_conditions or 'EW' in args.plot_conditions:
            try:
                condition3 = df[f'{line}_SNR'] > args.SNR_thresh
                set_arr, label_arr = make_set(df, condition1 & condition2 & condition3, f'{line} SNR > {args.SNR_thresh}', set_arr, label_arr)
            except:
                print(f'Could not apply the {line} SNR condition')
                pass
        if 'EW' in args.plot_conditions:
            try:
                condition4 = df[f'{line}_EW'] > args.EW_thresh
                set_arr, label_arr = make_set(df, condition1 & condition2 & condition3 & condition4, f'{line} EW_r > {args.EW_thresh}; SNR > {args.SNR_thresh}', set_arr, label_arr)
            except:
                print(f'Could not apply the {line} EW condition')
                pass

    # ---------add magnitude set------------
    if args.mag_lim is None: mag_lim = 26
    else: mag_lim = args.mag_lim
    condition = df['mag'] <= mag_lim
    set_arr, label_arr = make_set(df, condition, f'mag <= {mag_lim}', set_arr, label_arr)

    # ---------add semi-major (a_image) axis set------------
    condition = df['a_image'] >= args.a_thresh
    set_arr, label_arr = make_set(df, condition, f'a_image >= {args.a_thresh}', set_arr, label_arr)

    # ------add redshift range set-----------
    condition = df['redshift'].between(args.zmin, args.zmax)
    set_arr, label_arr = make_set(df, condition, f'{args.zmin}<z<{args.zmax}', set_arr, label_arr)

    # ------add number of grism orientations set-----------
    condition = df['nPA'] == 2
    set_arr, label_arr = make_set(df, condition, '#PA = 2', set_arr, label_arr)

    # ------add number of filter sets-----------
    if 'filters' in df:
        filter_arr = ['F200W', 'F150W', 'F115W']

        for index, filter in enumerate(filter_arr):
            condition = df['n_filters'] == index + 1
            set_arr, label_arr = make_set(df, condition, f'n_filters = {index + 1}', set_arr, label_arr)

            condition = df['filters'].str.contains(filter)
            set_arr, label_arr = make_set(df, condition, f'has {filter}', set_arr, label_arr)

    # ---------add sets from visual inspection------------
    if 'Notes' in df:
        print('\n')
        for attribute in ['compact', 'tail', 'merging', 'neighbour', 'clumpy', 'bulge', 'pea', 'bar', 'mg']:
            condition = df['Notes'].str.contains(attribute)
            set_arr, label_arr = make_set(df, condition, attribute, set_arr, label_arr)

        print('\n')
        for strong_line in ['OIII', 'Ha']:
            condition = (df[f'{strong_line} emission'].str.contains('strong')) & (df_visual['OIII emission'].str != np.nan)
            set_arr, label_arr = make_set(df, condition, f'strong_{strong_line}', set_arr, label_arr)

        print('\n')
        condition = df['DQ/RQ'].str.contains('okay')
        set_arr, label_arr = make_set(df, condition, 'RQ = okay', set_arr, label_arr)

    # ---------add sets from cosmos dataset------------
    if 'lp_mass' in df:
        print('\n')
        condition = np.isfinite(df['lp_mass'])
        set_arr, label_arr = make_set(df, condition, 'mass available', set_arr, label_arr)

        condition = df['lp_SFR'] > args.log_SFR_thresh
        set_arr, label_arr = make_set(df, condition, f'log sfr > {args.log_SFR_thresh}', set_arr, label_arr)

        condition = df['lp_sSFR'] > args.log_sSFR_thresh
        set_arr, label_arr = make_set(df, condition, f'log sSFR > {args.log_sSFR_thresh}', set_arr, label_arr)

    # ----------plot the venn diagrams----------
    which_sets_to_plot = [np.array([item1 in item2 for item1 in args.plot_conditions]).any() for item2 in label_arr]
    cmap = 'plasma'
    set_arr = np.array(set_arr)[which_sets_to_plot]
    label_arr = np.array(label_arr)[which_sets_to_plot]
    dataset_dict = dict(zip(label_arr, set_arr))

    # ---------manually calling draw_venn() so as to modify petal labels (for 0 counts)----------
    petal_labels = generate_petal_labels(dataset_dict.values(), fmt="{size}")
    petal_labels = {logic: value if int(value) > 0 else '' for logic, value in petal_labels.items()}
    colors = generate_colors(cmap=cmap, n_colors=len(label_arr))
    ax = draw_venn(petal_labels=petal_labels, dataset_labels=dataset_dict.keys(), hint_hidden=False, colors=colors, figsize=(8, 6), fontsize=args.fontsize, legend_loc='lower left', ax=None)

    # ---------printing the conditions for each non-zero petal label----------
    print('\n')
    for index in range(len(petal_labels.keys())):
        if list(petal_labels.values())[index] != '':
            bool_arr = list(petal_labels.keys())[index]
            cond_arr = []
            for index2 in range(len(bool_arr)):
                if int(bool_arr[index2]) == 1:
                    cond_arr.append(label_arr[index2])
            print(f'{int(list(petal_labels.values())[index])} objects in: {" and ".join(cond_arr)}')

    # -------calling the wrapper function, with automatic petal labelling (as opposed to manual calling above) but then 0 counts are displayed as such-------
    #ax = venn(dataset_dict, cmap=cmap, fmt='{size}', fontsize=8, legend_loc='upper left', ax=None)

    # ----------annotate and save the diagram----------
    fig = ax.figure
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    fig.text(0.99, 0.99, f'Total {n_fields} fields: Par{args.field_text}\nTotal {len(df)} objects', c='k', ha='right', va='top', transform=ax.transAxes)
    figname = args.output_dir / f'Par{args.field_text}_venn_diagram_{",".join(args.plot_conditions)}.png'

    # --------for talk plots--------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    # ----------deriving the dataframe corresponding to the innermost intersection----------
    intersecting_set = set.intersection(*set_arr)
    if len(intersecting_set) > 0:
        intersecting_par_obj = np.transpose([item.split('-') for item in list(intersecting_set)])
        df_int = pd.DataFrame({'field': intersecting_par_obj[0], 'objid':intersecting_par_obj[1]})
        df_int['objid'] = df_int['objid'].astype(int)
        df_int = df.merge(df_int, on=['field', 'objid'], how='inner')
        df_int.drop('par_obj', axis=1, inplace=True)
        if 'NUMBER' in df_int: df_int.drop('NUMBER', axis=1, inplace=True)
    else:
        df_int = pd.DataFrame()

    return df_int

# -------------------------------------------------------------------------------------------------------
def plot_nested_pie(df_subset, args, outer_col='nPA', inner_col='filters'):
    '''
    To plot concentric pie chart given a df, for a bunch of criteria
    Plots and saves the figure
    '''

    # -----------setting up dataframe--------
    n_fields = len(pd.unique(df_subset['field']))
    df_subset['nPA'] = df_subset['nPA'].apply(lambda x: f'{x} PA' if x <= 1 else f'{x} PAs')
    df = df_subset.pivot_table('par_obj', [outer_col, inner_col], aggfunc='count').reset_index()
    n_outer = len(pd.unique(df[outer_col]))
    n_inner = df.groupby(outer_col).size().values

    # ------------determining pie values-----------
    size = 0.3
    vals = df['par_obj'].values
    group_sum = df.groupby(outer_col).sum()['par_obj'].values # Major category values = sum of minor category values

    # -----------determining labels to be displayed------------
    outer_labels = pd.unique(df[outer_col])
    inner_labels = np.hstack([df[df[outer_col] == item][inner_col].values for item in pd.unique(df[outer_col])])

    # ----------determining colorbars------------
    color_stretch_arr = np.linspace(0.5, 0.1, 5)
    cmap_arr = [plt.cm.spring, plt.cm.summer, plt.cm.autumn, plt.cm.winter, plt.cm.pink]
    outer_colors = [item(0.6) for item in cmap_arr[:n_outer]]
    inner_colors = np.vstack([[cmap_arr[item](index) for index in color_stretch_arr[:n_inner[item]]] for item in range(n_outer)])

    # -------------function to format label inside wedges--------------
    def func(pct, allvals):
        absolute = int(pct / 100. * np.sum(allvals))
        return f'{pct :.1f}%\n({absolute:d})'

    # -----------function to plot one ring of the pie chart------------------
    def plot_pie_ring(ax, vals, radius, colors, labels, size):
        wedges, texts, counts = ax.pie(vals, radius=radius, colors=colors, labels=labels, autopct=lambda pct: func(pct, vals), wedgeprops=dict(width=size, edgecolor='w'))

        for text, count, wedge in zip(texts, counts, wedges):
            angle = np.mean([wedge.theta1, wedge.theta2])
            r = wedge.r - wedge.width / 2  # convert polar to cartesian
            x = r * np.cos(np.deg2rad(angle))  #
            y = r * np.sin(np.deg2rad(angle))  #
            text.set_x(x)
            text.set_y(y)
            text.set_rotation(angle - 90)

            r = wedge.r - wedge.width / 2 - 0.07  # convert polar to cartesian
            x = r * np.cos(np.deg2rad(angle))  #
            y = r * np.sin(np.deg2rad(angle))  #
            count.set_x(x)
            count.set_y(y - 0.05)
            count.set_rotation(angle - 90)

        return ax

    # -----------plotting------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    #ax = plot_pie_ring(ax, group_sum, 1, outer_colors, outer_labels, size)
    #ax = plot_pie_ring(ax, vals, 1-size, inner_colors, inner_labels, 1-size)
    ax = plot_pie_ring(ax, group_sum, 1-size, outer_colors, outer_labels, 1-size)
    ax = plot_pie_ring(ax, vals, 1, inner_colors, inner_labels, size)

    # ----------annotate and save the diagram----------
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
    fig.text(0.99, 0.99, f'Total {n_fields} fields: Par{args.field_text}\nTotal {len(df_subset)} objects', c='k', ha='right', va='top', transform=ax.transAxes)
    figname = args.output_dir / f'Par{args.field_text}_pie_diagram.png'

    # --------for talk plots--------------
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

# -------------------------------------------------------------------------------------------------------
def plot_sunburst_pie(df, args, outer_col='nPA', inner_col='filters'):
    '''
    To plot concentric pie chart given a df, for a bunch of criteria
    Plots and saves the figure
    '''
    df['count'] = 1
    df['nPA'] = df['nPA'].apply(lambda x: f'{x} PA' if x <= 1 else f'{x} PAs')
    fig = px.sunburst(df, path=[inner_col, outer_col], values='count', hover_data=['count'])
    fig.update_traces(textinfo='label+value+percent entry', texttemplate='<b>%{label}</b><br>(%{value})<br>%{percentEntry:.1%}')
    if args.fortalk: fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)'})


    # ----------annotate and save the diagram----------
    figname = args.output_dir / f'Par{args.field_text}_sunburst_diagram.png'
    fig.write_image(figname)
    print(f'\nSaved figure as {figname}')
    #fig.show()

    return fig

# -------------------------------------------------------------------------------------------------------
def plot_columns(df, args):
    '''
    To plot any column vs any column colorcoded by any column for a given df
    Plots and saves the figure
    Returns figure handle
    '''
    if args.fontsize == 10: args.fontsize = 15
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    p = ax.scatter(df[args.xcol], df[args.ycol], c=df[args.colorcol], s=5, lw=0, cmap='viridis')

    cbar = plt.colorbar(p)
    cbar.set_label(label_dict[args.colorcol] if args.colorcol in label_dict else args.colorcol, fontsize=args.fontsize)

    ax.set_xlabel(label_dict[args.xcol] if args.xcol in label_dict else args.xcol, fontsize=args.fontsize)
    ax.set_ylabel(label_dict[args.ycol] if args.ycol in label_dict else args.ycol, fontsize=args.fontsize)

    if args.xcol in bounds_dict: ax.set_xlim(bounds_dict[args.xcol][0], bounds_dict[args.xcol][1])
    if args.ycol in bounds_dict: ax.set_ylim(bounds_dict[args.ycol][0], bounds_dict[args.ycol][1])

    ax.set_xticklabels(['%.1f' % item for item in ax.get_xticks()], fontsize=args.fontsize)
    ax.set_yticklabels(['%.1f' % item for item in ax.get_yticks()], fontsize=args.fontsize)

    # ---------diagonal line---------
    if 'redshift' in args.xcol.lower() and 'z' in args.ycol.lower():
        line = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
        ax.plot(line, line, ls='dashed', c='k', lw=1)

    fig.text(0.99, 0.99, f'Total {len(pd.unique(df["field"]))} fields: Par{args.field_text}\nTotal {len(df)} objects', c='k', ha='right', va='top', transform=ax.transAxes)
    figname = args.output_dir / f'all_df_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig


# -------------------------------------------------------------------------------------------------------
def get_detection_fraction(df, line, args):
    '''
    To compute fraction of objects in a given field that have a given line detected (beyond a given EW threshold)
    out of all objects where the line was accessible in the given filter wavelength regime
    Plots and saves EW distribution as histograms
    Returns subset dataframe where line is detected
    '''
    color, mag_color = 'cornflowerblue', 'salmon'

    df = df[(np.isfinite(df[f'{line}_EW'])) & (df[f'{line}_EW'] > 0)]
    df_detected = df[df[f'{line}_EW'] > args.EW_thresh]

    # ------plot histogram of line EWs--------
    fig, ax = plt.subplots()

    hist = ax.hist(np.log10(df[f'{line}_EW']), bins=50, range=(-2, 6), color=color, histtype='step', lw=2, label='all')
    ax.axvline(np.log10(args.EW_thresh), c='k', lw=2, ls='dashed')
    ax.text(0.02, 0.75, f'{len(df_detected)} above EW {args.EW_thresh} A', c=color, ha='left', va='top', transform=ax.transAxes)

    # -----magnitude limit-----------
    if args.mag_lim is not None:
        df_magcut = df[df['mag'] <= args.mag_lim]
        df_magcut_detected = df_magcut[df_magcut[f'{line}_EW'] > args.EW_thresh]
        hist2 = ax.hist(np.log10(df_magcut[f'{line}_EW']), bins=hist[1], color=mag_color, histtype='step', lw=2, label=f'mag <= {args.mag_lim}')
        ax.text(0.02, 0.7, f'{len(df_magcut_detected)} above EW {args.EW_thresh} A', c=mag_color, ha='left', va='top', transform=ax.transAxes)

    ax.legend(loc='upper left')
    ax.set_ylabel('#objects')
    ax.set_xlabel(f'log {line} EW (A)')
    #ax.set_xlim(-2, 6)

    if args.do_all_fields:
        ax.text(0.02, 0.8, f'Par{args.field_text}', c='k', ha='left', va='top', transform=ax.transAxes)
        figname = args.output_dir / f'Par{args.field_text}_{line}_EW_histogram.png'
    else:
        ax.text(0.02, 0.8, f'Par{args.field_text}', c='k', ha='left', va='top', transform=ax.transAxes)
        figname = args.output_dir / f'{args.field}' / f'{args.field}_{line}_EW_histogram.png'

    # --------for talk plots--------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    if args.mag_lim is None: return df_detected
    else: return df_magcut_detected

# -------------------------------------------------------------------------------------------------------
def read_stats_df(df_filename, args):
    '''
    To read in the dataframe produced by make_diagnostic_plots.py
    Returns dataframe
    '''
    extract_dir = args.input_dir / args.field / 'Extractions'
    # -------initiliasing dataframe-------------------------------
    df = pd.read_table(df_filename, delim_whitespace=True)

    # ------------getting magnitudes from phot catalog----------------------
    try:
        phot_catalog_file = extract_dir / f'{args.field}_phot.fits'
        phot_catalog = GTable.read(phot_catalog_file)
    except:
        try:
            phot_catalog_file = args.input_dir / args.field / 'Products' / f'{args.field}_photcat.fits'
            phot_catalog = GTable.read(phot_catalog_file)
        except:
            phot_catalog = None

    if phot_catalog is not None:
        print(f'Merging photcat for {args.field}..')
        phot_catalog.rename_column('id', 'objid')
        phot_catalog.rename_column('mag_auto', 'mag')
        phot_catalog_df = phot_catalog['objid', 'mag', 'a_image', 'b_image'].to_pandas()
        df = df.merge(phot_catalog_df, on='objid', how='inner')
    else:
        df['mag'] = np.nan

    if args.field in fields_with_2PA: df['nPA'] = 2
    else: df['nPA'] = 1

    # ------------getting EW error and SNR from spec catalog----------------------
    grab_for_each_line = ['ewhw', 'sn']
    columns_to_extract = ['objid']
    line_list = [item for item in list(rest_wave_dict.keys()) if '+' not in item]
    for line in line_list: columns_to_extract += [f'{item}_{line}' for item in grab_for_each_line]

    try:
        spec_catalog_file = args.input_dir / args.field / 'Products' / f'{args.field}_speccat.fits'
        spec_catalog = GTable.read(spec_catalog_file)
    except:
        spec_catalog = None

    if spec_catalog is not None:
        print(f'Merging speccat for {args.field}..')
        spec_catalog.rename_column('id', 'objid')
        spec_catalog_df = spec_catalog[columns_to_extract].to_pandas()
        df = df.merge(spec_catalog_df, on='objid', how='inner')
    else:
        for col in columns_to_extract[1:]: df[col] = np.nan

    df = df.rename(columns={f'sn_{item}':f'{item}_SNR' for item in line_list})
    df = df.rename(columns={f'ewhw_{item}':f'{item}_EW_u' for item in line_list})

    return df

# -------------------------------------------------------------------------------------------------------
def read_visual_df(args):
    '''
    To read in the Google Sheet that resides online and was produced by visual inspection
    Returns dataframe
    '''
    sheetname = args.field[:3] + str(int(args.field[3:]))
    gsheet_url = f'https://docs.google.com/spreadsheets/d/1TcIjeP_BgjAQ4ABMz5nVUCb3hZjsSSW4-FcLceLdJ7U/export?&format=xlsx'

    try:
        df = pd.read_excel(gsheet_url, sheetname,  skiprows=4 if '51' in args.field else 2)
        df['field'] = args.field
    except ValueError:
        print(f'Field {args.field} has no visual inspection sheet')
        df = pd.DataFrame()

    if 'z1' in df and 'z2' in df:
        df.rename(columns={'z2':'z'}, inplace=True)
        df.drop('z1', axis=1, inplace=True)

    df.rename(columns={'ID':'objid'}, inplace=True)

    return df

# ----------------------------------------------------------------------------------------------------------
def get_crossmatch_with_cosmos(df, args):
    '''
    Determines crossmatch of a given dataframe with COSMOS2020 catalog
    Returns input dataframe merged with COSMOS stellar masses, etc. wherever available
    '''
    df = df.copy()
    fields = pd.unique(df['field'])
    if 'par_obj' not in df: df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype(str)  # making a unique combination of field and object id
    df = df.drop_duplicates('par_obj', keep='last')

    # -------collating only those COSMOS objects that lie within the FoV of available PASSAGE fields------
    print(f'\nTrying to read in COSMOS catalogs..')
    df_cosmos = pd.DataFrame()

    for index, thisfield in enumerate(fields):
        filename = args.input_dir / 'COSMOS' / f'cosmos2020_objects_in_{thisfield}.fits'
        if os.path.exists(filename):
            print(f'{index+1} of {len(fields)} fields: Reading COSMOS subset table from {filename}')
            df_cosmos_thisfield = read_COSMOS2020_catalog(filename=filename)
            df_cosmos = pd.concat([df_cosmos, df_cosmos_thisfield])
        else:
            print(f'{index+1} of {len(fields)} fields: Could not find COSMOS subset table for {thisfield}, so skipping.')
    if 'passage_id' in df_cosmos: df_cosmos = df_cosmos.drop('passage_id', axis=1)

    # -------cross-matching RA/DEC of both catalogs------
    print(f'\nDoing cross-matching between PASSAGE and COSMOS catalogs..')
    passage_coords = SkyCoord(df['ra'], df['dec'], unit='deg')
    cosmos_coords = SkyCoord(df_cosmos['ra'], df_cosmos['dec'], unit='deg')
    nearest_id_in_cosmos, sep_from_nearest_id_in_cosmos, _ = passage_coords.match_to_catalog_sky(cosmos_coords)

    df_crossmatch = pd.DataFrame({'passage_id': df['par_obj'].values, 'cosmos_id': df_cosmos['id'].iloc[nearest_id_in_cosmos].values, 'sep': sep_from_nearest_id_in_cosmos.arcsec})
    df_crossmatch = df_crossmatch[df_crossmatch['sep'] < 1.]  # separation within 1 arcsecond
    df_crossmatch = df_crossmatch.sort_values('sep').drop_duplicates(subset='cosmos_id', keep='first').reset_index(drop=True)  # to avoid multiple PASSAGE objects being linked to the same COSMOS object
    df_crossmatch = pd.merge(df_crossmatch[['passage_id', 'cosmos_id']], df_cosmos, left_on='cosmos_id', right_on='id', how='inner').drop(['id', 'ra', 'dec'], axis=1)

    try: df_crossmatch_subset = df_crossmatch[['passage_id', 'cosmos_id', 'ez_z_phot', 'lp_MK', 'lp_zBEST']]
    except: df_crossmatch_subset = df_crossmatch[['passage_id', 'cosmos_id', 'ez_z_phot', 'lp_MK']]
    lp_quants = ['lp_SFR', 'lp_mass', 'lp_sSFR']
    ez_quants = ['ez_mass', 'ez_sfr', 'ez_ssfr']

    for quant in lp_quants:
        df_crossmatch_subset[quant] = df_crossmatch[quant + '_med']
        df_crossmatch_subset[quant + '_u'] = (df_crossmatch[quant + '_med_max68'] - df_crossmatch[quant + '_med_min68']) / 2

    for quant in ez_quants:
        df_crossmatch_subset[quant] = df_crossmatch[quant + '_p500']
        df_crossmatch_subset[quant + '_u'] = (df_crossmatch[quant + '_p840'] - df_crossmatch[quant + '_p160']) / 2

    print(f'\nFound a total of {len(df_crossmatch_subset)} matching objects')
    df = pd.merge(df, df_crossmatch_subset, left_on='par_obj', right_on='passage_id', how='outer').drop('passage_id', axis=1)

    return df

# --------------------------------------------------------------------------------------------------------------------
label_dict = {'lp_mass': r'log M$_*$/M$_{\odot}$', 'ez_mass': r'log M$_*$/M$_{\odot}$', 'lp_SFR': r'log SFR (M$_{\odot}$/yr)', 'ez_z_phot': 'Photo-z COSMOS2020 (EAZY)', 'lp_zBEST': 'Photo-z COSMOS2020 (Le Phare)', 'redshift': 'Spec-z (Grizli)', 'logOH_slope':r'log $\nabla$Z$_r$ (dex/kpc)'}
bounds_dict = {'lp_mass': (6, 12), 'ez_mass': (6, 12), 'lp_SFR': (-3, 1), 'log_SFR_int': (-3, 1), 'ez_z_phot': (0, 6), 'lp_zBEST': (0, 6), 'redshift': (0, 6), 'logOH_slope': (-0.4, 0.1)}
colormap_dict = defaultdict(lambda: 'viridis', ez_z_phot='plasma')


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ---------setting up fields and paths----------------
    lines_to_consider = args.line_list # ['OII', 'OIII'] # OR
    if args.do_all_fields:
        available_fields = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.output_dir / 'Par*') + '/')]
        available_fields.sort(key=natural_keys)
    else:
        available_fields = args.field_arr

    df_stats_filename = args.output_dir / f'all_fields_diag_results.txt'
    df_visual_filename = args.output_dir / f'all_fields_visual_inspection_results.txt'
    df_outfilename = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df.txt'

    # -------------------------------------------------------------------
    if not os.path.exists(df_outfilename) or args.clobber_venn_df:
        if args.do_all_fields and os.path.exists(df_stats_filename) and not args.clobber:
            print(f'Reading in existing {df_stats_filename}')
            df_stats = pd.read_table(df_stats_filename, delim_whitespace=True)

            print(f'Reading in existing {df_visual_filename}')
            df_visual = pd.read_csv(df_visual_filename)

        else:
            df_stats = pd.DataFrame()
            df_visual = pd.DataFrame()

            # ---------looping over fields-----------
            for index, args.field in enumerate(available_fields):
                print(f'Doing field {args.field} which is {index+1} of {len(available_fields)}..')
                args.filters = available_filters_for_field_dict[args.field]

                # ---------determining filename suffixes-------------------------------
                output_dir = args.output_dir / args.field
                if args.re_extract: output_dir = output_dir / 're_extracted'
                df_filename = output_dir / f'{args.field}_all_diag_results.txt'

                if os.path.exists(df_filename):
                    thisdf_stat = read_stats_df(df_filename, args)  # read in the stats dataframe
                    df_stats = pd.concat([df_stats, thisdf_stat], ignore_index=True)

                    thisdf_visual = read_visual_df(args)  # read in the visually inspected dataframe
                    df_visual = pd.concat([df_visual, thisdf_visual], ignore_index=True)
                else:
                    print(f'{df_filename} does not exists. Skipping this field.')

            # ------------saving the master dfs--------------------
            if args.do_all_fields:
                df_stats.to_csv(df_stats_filename, sep='\t', index=None, na_rep='NULL')
                print(f'Saved master stats df in {df_stats_filename}')

                df_visual.to_csv(df_visual_filename, index=None, na_rep='NULL')
                print(f'Saved master stats df in {df_visual_filename}')

        # ------------doing the line histograms--------------------
        if args.plot_EW_hist:
            for index, line in enumerate(lines_to_consider):
                print(f'Doing line {line} which is {index+1} of {len(lines_to_consider)}..')
                df_detected = get_detection_fraction(df_stats, line, args)

        # ------------merging visual dataframes for the venn diagrams--------------------
        conditions_from_visual = ['compact', 'tail', 'merging', 'neighbour', 'clumpy', 'bulge', 'pea', 'bar', 'mg', 'RQ']
        if args.merge_visual or len(set(conditions_from_visual).intersection(set(args.plot_conditions))) > 0:
            df_visual.drop('nPA', axis=1, inplace=True)
            df = pd.merge(df_stats, df_visual, on=['field', 'objid'], how='inner')
        else:
            df = df_stats

        has_fields = [str(int(item[3:])) for item in pd.unique(df_stats['field'])]
        has_fields.sort(key=natural_keys)
        args.field_text = ','.join(has_fields)

        # ------------adding some common columns before crossmatching--------------------
        df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype(str)  # making a unique combination of field and object id
        df = df.drop_duplicates('par_obj', keep='last')
        df['filters'] = df['field'].map(lambda x: ', '.join(available_filters_for_field_dict[x]))
        df['n_filters'] = df['filters'].map(lambda x: len(x.split(',')))
        EW_cols = [item for item in df.columns if 'EW' in item]
        for thiscol in EW_cols: df[thiscol] = df[thiscol] / (1 + df['redshift']) # converting ALL EWs to rest-frame

        # ------------merging cosmos datasets for the venn diagrams--------------------
        conditions_from_cosmos = ['mass', 'sfr', 'sSFR']
        if len(set(conditions_from_cosmos).intersection(set(args.plot_conditions))) > 0:
            df = get_crossmatch_with_cosmos(df, args)
            df = df[df['ez_mass'] < 20] # masses cannot be over 10^20 Msun

        # ------------doing the photo-z vs spec-z comparison--------------------
        if args.plot_columns:
            fig = plot_columns(df, args)
        # ------------doing the pie charts--------------------
        elif args.plot_pie:
            plot_nested_pie(df, args, outer_col='filters', inner_col='nPA')
        elif args.plot_sunburst:
            fig = plot_sunburst_pie(df, args, outer_col='nPA', inner_col='filters')
        else:
            # ------------doing the venn diagrams--------------------
            df_int = plot_venn(df, args)

            # ------------saving the resultant intersecting dataframe--------------------
            df_int.to_csv(df_outfilename, index=None)
            print(f'Saved intersecting dataframe as {df_outfilename}')

    else:
        print(f'Reading from existing {df_outfilename}')
        df_int = pd.read_csv(df_outfilename)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
