'''
    Filename: get_field_stats.py
    Notes: Computes and plots various statistics (redshift distribution, detected line distribution, magnitude distribution etc.) for a given PASSAGE field which has been already run through make_diagnostic_maps.py
    Author : Ayan
    Created: 19-08-24
    Example: run get_field_stats.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --re_extract
             run get_field_stats.py --field Par61 --mag_lim 26 --line_list OII,OIII,Ha
             run get_field_stats.py --mag_lim 26 --line_list OIII --do_all_fields --clobber
             run get_field_stats.py --mag_lim 26 --line_list OIII,Ha --do_all_fields --plot_venn --zmin 1 --zmax 2.5 --merge_visual --plot_conditions detected,z,mag,tail
'''
from header import *
from util import *

start_time = datetime.now()

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

# -------------------------------------------------------------------------------------------------------
def plot_venn(df, args):
    '''
    To plot Venn diagrams with a given df, for a bunch of criteria
    Plots and saves the figure
    Returns figure handle
    '''
    df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype(str)

    set_arr = []
    label_arr = []

    # ---------add line sets------------
    line_list = args.line_list

    for line in line_list:
        condition1 = (np.isfinite(df[f'{line}_EW'])) & (df[f'{line}_EW'] > 0)
        df_line = df[condition1]
        set_arr, label_arr = make_set(df, condition1, f'{line}_present', set_arr, label_arr)

        condition2 = df_line[f'{line}_EW'] > args.EW_thresh
        set_arr, label_arr = make_set(df_line, condition2, f'{line}_detected', set_arr, label_arr)

    # ---------add magnitude set------------
    if args.mag_lim is None: mag_lim = 26
    else: mag_lim = args.mag_lim
    condition = df['mag'] <= mag_lim
    set_arr, label_arr = make_set(df, condition, f'mag <= {mag_lim}', set_arr, label_arr)

    # ---------add sets from visual inspection------------
    if 'Notes' in df:
        for attribute in ['compact', 'tail', 'strong', 'merging', 'neighbour', 'clumpy', 'bulge', 'pea', 'bar']:
            condition1 = df['Notes'].str.contains(attribute)
            set_arr, label_arr = make_set(df, condition1, attribute, set_arr, label_arr)

    # ------add redshift range set-----------
    condition = df['redshift'].between(args.zmin, args.zmax)
    set_arr, label_arr = make_set(df, condition, f'{args.zmin}<z<{args.zmax}', set_arr, label_arr)

    # ----------plot the enn diagrams----------
    which_sets_to_plot = [np.array([item1 in item2 for item1 in args.plot_conditions]).any() for item2 in label_arr]
    cmap = 'plasma'
    set_arr = np.array(set_arr)[which_sets_to_plot]
    label_arr = np.array(label_arr)[which_sets_to_plot]
    dataset_dict = dict(zip(label_arr, set_arr))

    fig = plt.figure()
    ax = plt.gca()

    venn(dataset_dict, cmap=cmap, fmt='{size}', fontsize=8, legend_loc='upper left', ax=ax)

    if args.do_all_fields:
        fig.text(0.99, 0.99, f'Par{args.field_text}', c='k', ha='right', va='top', transform=ax.transAxes)
        figname = args.output_dir / f'Par{args.field_text}_venn_diagram.png'
    else:
        fig.text(0.99, 0.99, args.field, c='k', ha='right', va='top', transform=ax.transAxes)
        figname = args.output_dir / f'{args.field}' / f'{args.field}_venn_diagram.png'

    fig.savefig(figname)
    print(f'Saved figure as {figname}')
    plt.show(block=False)


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
        ax.text(0.02, 0.8, args.field, c='k', ha='left', va='top', transform=ax.transAxes)
        figname = args.output_dir / f'{args.field}' / f'{args.field}_{line}_EW_histogram.png'

    fig.savefig(figname)
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
    # ------------getting magnitudes from catalog----------------------
    try:
        catalog_file = extract_dir / f'{args.field}-ir.cat.fits'
        catalog = GTable.read(catalog_file)
    except:
        try:
            catalog_file = args.input_dir / args.field / 'Products' / f'{args.field}_photcat.fits'
            catalog = GTable.read(catalog_file)
        except:
            catalog = None

    # -------initiliasing dataframe-------------------------------
    df = pd.read_table(df_filename, delim_whitespace=True)
    if catalog is not None and 'MAG_AUTO' in catalog.columns:
        catalog_df = catalog['NUMBER', 'MAG_AUTO'].to_pandas()
        df = df.merge(catalog_df, left_on='objid', right_on='NUMBER', how='inner')
        df.rename(columns={'MAG_AUTO':'mag'}, inplace=True)
    else:
        df['mag'] = np.nan

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

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ---------setting up fields and paths----------------
    lines_to_consider = args.line_list # ['OII', 'OIII'] # OR
    if args.do_all_fields:
        available_fields = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.output_dir / 'Par*') + '/')]
    else:
        available_fields = [args.field]

    df_stats_filename = args.output_dir / f'all_fields_diag_results.txt'
    df_visual_filename = args.output_dir / f'all_fields_visual_inspection_results.txt'

    # -------------------------------------------------------------------
    if args.do_all_fields and os.path.exists(df_stats_filename) and not args.clobber:
        print(f'Reading in existing {df_stats_filename}')
        df_stats = pd.read_table(df_stats_filename, delim_whitespace=True)

        print(f'Reading in existing {df_visual_filename}')
        df_visual = pd.read_csv(df_visual_filename)

        args.field_text = ','.join([str(int(item[3:])) for item in pd.unique(df_stats['field'])])
    else:
        df_stats = pd.DataFrame()
        df_visual = pd.DataFrame()
        args.field_text = ''

        # ---------looping over fields-----------
        for index, args.field in enumerate(available_fields):
            print(f'Doing field {args.field} which is {index+1} of {len(available_fields)}..')
            args.filters = filter_dict[args.field]
            if '51' in args.field:
                args.re_extract = True
            else:
                args.re_extract = False

            # ---------determining filename suffixes-------------------------------
            output_dir = args.output_dir / args.field
            if args.re_extract: output_dir = output_dir / 're_extracted'
            df_filename = output_dir / f'{args.field}_all_diag_results.txt'

            if os.path.exists(df_filename):
                thisdf_stat = read_stats_df(df_filename, args)  # read in the stats dataframe
                df_stats = pd.concat([df_stats, thisdf_stat], ignore_index=True)
                args.field_text += f'{int(args.field[4:])},'

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
    for index, line in enumerate(lines_to_consider):
        print(f'Doing line {line} which is {index+1} of {len(lines_to_consider)}..')
        df_detected = get_detection_fraction(df_stats, line, args)

    # ------------doing the venn diagrams--------------------
    if args.merge_visual: df = pd.merge(df_stats, df_visual, on=['field', 'objid'], how='inner')
    else: df = df_stats
    plot_venn(df, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
