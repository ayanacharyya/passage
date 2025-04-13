'''
    Filename: plot_mappings_grid.py
    Notes: Makes line ratio grids from MAPPINGS photoionisation models
    Author : Ayan
    Created: 13-03-25
    Example: run plot_mappings_grid.py --plot_grid --ynum_line OIII --yden_line Hb --xnum_line SII --xden_line Ha,NII --slice_at_quantity3 5 --annotate
             run plot_mappings_grid.py --plot_model --ynum_line OIII --yden_line Hb --slice_at_quantity3 4,5,6
             run plot_mappings_grid.py --plot_grid --ynum_line OIII --yden_line Hb --xnum_line NeIII --xden_line OII --fit_y_envelope
             run plot_mappings_grid.py --plot_grid --ynum_line OIII --yden_line Hb --xnum_line NeIII --xden_line OII --xmin -3.5 --xmax 1 --ymin -4 --ymax 1 --keep --phot_model nb
'''
from header import *
from util import *

start_time = datetime.now()

# ---------plots a grid in 2 ratios--------------------------------------
def plot_ratio_grid(df_ratios, ax, args, color1='salmon', color2='cornflowerblue', color3='sienna'):
    '''
    Plots the grid of a specific given line ratio for a given value of Z, q, P, on an existing axis handle
    Returns the axis handle and ratio names
    '''
    if args.phot_models.lower() in ['mappings', 'map']: line_names_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,9', 'NeIII':'[NeIII]3869', 'NeIII-3867':'[NeIII]3869', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31'}) # to map between user input line labels and line labels used in ratio_list.txt file
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: line_names_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})

    # ------------getting the line ratio names------------------
    x_num_labels = ','.join([line_names_dict[item] for item in args.xnum_line.split(',')])
    x_den_labels = ','.join([line_names_dict[item] for item in args.xden_line.split(',')])
    xratio_name = f'{x_num_labels}/{x_den_labels}'

    y_num_labels = ','.join([line_names_dict[item] for item in args.ynum_line.split(',')])
    y_den_labels = ','.join([line_names_dict[item] for item in args.yden_line.split(',')])
    yratio_name = f'{y_num_labels}/{y_den_labels}'

    if args.slice_at_quantity1 is not None:
        print(f'Slicing the model at {args.quantity1} == {args.slice_at_quantity1}')
        df_ratios = df_ratios[df_ratios[args.quantity1].isin(args.slice_at_quantity1)] # to plot grid for only one value of quantity2, to reduce clutter
    if args.slice_at_quantity2 is not None:
        print(f'Slicing the model at {args.quantity2} == {args.slice_at_quantity3}')
        df_ratios = df_ratios[df_ratios[args.quantity2].isin(args.slice_at_quantity2)] # to plot grid for only one value of quantity2, to reduce clutter
    if args.slice_at_quantity3 is not None:
        print(f'Slicing the model at {args.quantity3} == {args.slice_at_quantity3}')
        df_ratios = df_ratios[df_ratios[args.quantity3].isin(args.slice_at_quantity3)] # to plot grid for only one value of quantity3, to reduce clutter

    for i, quant3 in enumerate(np.unique(df_ratios[args.quantity3])):
        df_sub = df_ratios[df_ratios[args.quantity3] == quant3]

        for j, quant1 in enumerate(np.unique(df_sub[args.quantity1])):
            df_sub_sub = df_sub[df_sub[args.quantity1] == quant1]
            if not args.fit_y_envelope: ax.plot(np.log10(df_sub_sub[xratio_name]), np.log10(df_sub_sub[yratio_name]), color=color1, lw=1, alpha=(j+1)/len(np.unique(df_sub[args.quantity1])))
            else: ax.scatter(np.log10(df_sub_sub[xratio_name]), np.log10(df_sub_sub[yratio_name]), color=color1, s=5, alpha=(j+1)/len(np.unique(df_sub[args.quantity1])))
            if args.annotate and i == len(np.unique(df_ratios[args.quantity3])) - 1: ax.annotate(f'{args.quantity1} = {quant1:.2f}', xy=(np.log10(df_sub_sub[xratio_name].values[0]), np.log10(df_sub_sub[yratio_name].values[0])), xytext=(np.log10(df_sub_sub[xratio_name].values[0]) - 0.2, np.log10(df_sub_sub[yratio_name].values[0]) - 0.1), color=color1, fontsize=args.fontsize/1.5, arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color1))

        for j, quant2 in enumerate(np.unique(df_sub[args.quantity2])):
            df_sub_sub = df_sub[df_sub[args.quantity2] == quant2]
            if not args.fit_y_envelope: ax.plot(np.log10(df_sub_sub[xratio_name]), np.log10(df_sub_sub[yratio_name]), color=color2, lw=1, alpha=(j+1)/len(np.unique(df_sub[args.quantity2])))
            if args.annotate and i == len(np.unique(df_ratios[args.quantity3])) - 1: ax.annotate(f'{args.quantity2} = {quant2:.2f}', xy=(np.log10(df_sub_sub[xratio_name].values[-1]), np.log10(df_sub_sub[yratio_name].values[-1])), xytext=(np.log10(df_sub_sub[xratio_name].values[-1]) - 0.1, np.log10(df_sub_sub[yratio_name].values[-1]) + 0.1), color=color2, fontsize=args.fontsize/1.5, arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color2))

        if args.annotate: ax.annotate(f'{args.quantity3} = {quant3:.2f}', xy=(np.log10(df_sub[xratio_name].values[-1]), np.log10(df_sub[yratio_name].values[-1])), xytext=(np.log10(df_sub[xratio_name].values[-1]) + 0.3, np.log10(df_sub[yratio_name].values[-1]) - 0.2), color=color3, fontsize=args.fontsize / 1.5, arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color3))

    for i, quant1 in enumerate(np.unique(df_ratios[args.quantity1])):
        for j, quant2 in enumerate(np.unique(df_ratios[args.quantity2])):
            df_sub = df_ratios[(df_ratios[args.quantity1] == quant1) & (df_ratios[args.quantity2] == quant2)]
            if not args.fit_y_envelope: ax.plot(np.log10(df_sub[xratio_name]), np.log10(df_sub[yratio_name]), color=color3, lw=0.5, alpha=0.7)

    return ax, xratio_name, yratio_name

# --------------------------------------------------------------------------------------------------------------------
def plot_ratio_grid_fig(df_ratios, args):
    '''
    Wrapper for plot_ratio_grid()
    Returns the figure handle
    '''
    # ------declare the figure-----------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.1, top=0.95, wspace=0.2)

    # --------plot the model ratios---------------
    ax, xratio_name, yratio_name = plot_ratio_grid(df_ratios, ax, args)

    # ------annotate figure-----------------------
    ax.grid(which='both', color='gray', linestyle='solid', linewidth=1, alpha=0.3)
    ax.set_xlabel(f'Log {xratio_name}', fontsize=args.fontsize)
    ax.set_ylabel(f'Log {yratio_name}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    df_ratios = df_ratios[(df_ratios[xratio_name] > 0) & (df_ratios[yratio_name] > 0)] # to avoid math errors later while taking log
    xmin = args.xmin if args.xmin is not None else np.log10(np.min(df_ratios[xratio_name]) * 0.9)
    xmax = args.xmax if args.xmax is not None else np.log10(np.max(df_ratios[xratio_name]) * 1.1)
    ymin = args.ymin if args.ymin is not None else np.log10(np.min(df_ratios[yratio_name]) * 0.9)
    ymax = args.ymax if args.ymax is not None else np.log10(np.max(df_ratios[yratio_name]) * 1.1)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # ------for fitting to upper y-axis envelope of data------------
    if args.fit_y_envelope:
        nbins = 25

        def func(x, *popt): return np.poly1d(popt)(x)
        p_init_arr = [[1, 1, 1, 1]]
        func_arr = [func]

        xarr = np.linspace(np.log10(np.min(df_ratios[xratio_name])), np.log10(np.max(df_ratios[xratio_name])), nbins)
        df_ratios['bin'] = pd.cut(np.log10(df_ratios[xratio_name]), bins=xarr)
        grouped = df_ratios.groupby('bin')
        xbins = np.log10(grouped[xratio_name].mean().values)
        ybins = np.log10(grouped[yratio_name].max().values)
        ax.plot(xbins, ybins, c='r', lw=0.5)

        col_arr = ['brown', 'darkgreen', 'cornflowerblue']
        for index, (p_init, func) in enumerate(zip(p_init_arr, func_arr)):
            popt, pcov = curve_fit(func, xbins, ybins, p0=p_init)
            ax.plot(xarr, func(xarr, *p_init), c=col_arr[index], lw=1, ls='--')
            ax.plot(xarr, func(xarr, *popt), c=col_arr[index], lw=2)
            ax.text(0.98, 0.01 + index * 0.07, f'fit{index+1} = [{",".join([f"{item:.2f}" for item in popt])}]', c=col_arr[index], ha='right', va='bottom', fontsize=args.fontsize, transform=ax.transAxes)

    # --------for talk plots--------------
    if args.fortalk:
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    # ---------saving figure--------------
    if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: phot_model_text = 'NebulaBayes'
    slice_text = ''
    if args.slice_at_quantity1 is not None: slice_text += f'_{args.quantity1}={",".join([str(item) for item in args.slice_at_quantity1])}'.replace('/', '-')
    if args.slice_at_quantity2 is not None: slice_text += f'_{args.quantity2}={",".join([str(item) for item in args.slice_at_quantity2])}'.replace('/', '-')
    if args.slice_at_quantity3 is not None: slice_text += f'_{args.quantity3}={",".join([str(item) for item in args.slice_at_quantity3])}'.replace('/', '-')
    figname = args.mappings_dir / 'plots' / f'{phot_model_text}_grid_{yratio_name.replace("/", "-")}_va_{xratio_name.replace("/", "-")}{slice_text}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def plot_ratio_model(df_ratios, ax, args):
    '''
    Plots the given line ratio as a function of quantity1 for a given value of quantity2 and quantity3, on an existing axis handle
    Returns the axis handle and ratio name
    '''
    col_arr = ['k', 'r', 'b', 'orange', 'g', 'brown', 'darkblue', 'm', 'olive', 'darksalmon', 'teal']  # at least as many as "quantity3" values
    lighter_col_dict = {'k':'gray', 'r':'lightsalmon', 'b':'lightskyblue', 'orange':'wheat', 'g':'lightgreen'}
    style_arr = ['solid', 'dashed', 'dotted']
    thickness_arr = [1.2, 0.7, 0.3]  # len(style_arr) \times len(tickness_arr) = at least as many as "quantity2" values
    thickness_arr = np.tile(np.repeat(thickness_arr, len(thickness_arr)), 2)
    if args.phot_models.lower() in ['mappings', 'map']: line_names_dict = smart_dict({'OIII':'[OIII]5007', 'OII':'[OII]3727,29', 'Hb':'Hbeta', 'Ha':'Halpha', 'NII':'[NII]6584', 'SII':'[SII]6717,31', 'NeIII':'[NeIII]3869'}) # to map between user input line labels and line labels used in ratio_list.txt file
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: line_names_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})

    # ------------getting the line ratio names------------------
    num_labels = ','.join([line_names_dict[item] for item in args.ynum_line.split(',')])
    den_labels = ','.join([line_names_dict[item] for item in args.yden_line.split(',')])
    ratio_name = f'{num_labels}/{den_labels}'

    if args.slice_at_quantity2 is not None:
        print(f'Slicing the model at {args.quantity2} == {args.slice_at_quantity2}')
        df_ratios = df_ratios[df_ratios[args.quantity2].isin(args.slice_at_quantity2)] # to plot grid for only one value of quantity2, to reduce clutter
    if args.slice_at_quantity3 is not None:
        print(f'Slicing the model at {args.quantity3} == {args.slice_at_quantity3}')
        df_ratios = df_ratios[df_ratios[args.quantity3].isin(args.slice_at_quantity3)] # to plot grid for only one value of quantity3, to reduce clutter

    for i, quant2 in enumerate(np.unique(df_ratios[args.quantity2])):
        df_sub = df_ratios[df_ratios[args.quantity2] == quant2]

        for j, quant3 in enumerate(np.unique(df_sub[args.quantity3])):
            df_sub_sub = df_sub[df_sub[args.quantity3] == quant3]
            label_root = f'{args.quantity2} = {quant2}, {args.quantity3} = '
            label = f'{quant3:>{len(label_root)}}' if j else f'{label_root}{quant3}'
            ax.plot(df_sub_sub[args.quantity1], np.log10(df_sub_sub[ratio_name]), color=col_arr[i], linestyle=style_arr[j % len(style_arr)], linewidth=thickness_arr[j], label=label)

    return ax, ratio_name

# --------------------------------------------------------------------------------------------------------------------
def plot_ratio_model_fig(df_ratios, args):
    '''
    Wrapper for plot_ratio_model()
    Returns the figure handle
    '''
    # ------declare the figure-----------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95, wspace=0.2)

    # --------plot the model ratios---------------
    ax, ratio_name = plot_ratio_model(df_ratios, ax, args)

    # ------annotate figure-----------------------
    ax.grid(which='both', color='gray', linestyle='solid', linewidth=1, alpha=0.3)
    ax.legend(fontsize=args.fontsize/2)
    ax.set_xlabel('log(O/H) + 12' if args.quantity1 == 'Z' else args.quantity1, fontsize=args.fontsize)
    ax.set_ylabel(f'Log {ratio_name}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)


    xmin = args.xmin if args.xmin is not None else np.min(np.unique(df_ratios[args.quantity1])) * 0.9
    xmax = args.xmax if args.xmax is not None else  np.max(np.unique(df_ratios[args.quantity1])) * 1.1
    ymin = args.ymin if args.ymin is not None else np.log10(np.min(df_ratios[ratio_name]) * 0.9)
    ymax = args.ymax if args.ymax is not None else np.log10(np.max(df_ratios[ratio_name]) * 1.1)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # ------for fitting to upper y-axis envelope of data------------
    if args.fit_y_envelope:
        nbins = 25

        def func(x, *popt): return np.poly1d(popt)(x)
        p_init_arr = [[1, 1]]
        func_arr = [func]

        if args.slice_at_quantity2 is not None: df_ratios = df_ratios[df_ratios[args.quantity2].isin(args.slice_at_quantity2)]  # to plot grid for only one value of quantity2, to reduce clutter
        if args.slice_at_quantity3 is not None: df_ratios = df_ratios[df_ratios[args.quantity3].isin(args.slice_at_quantity3)]  # to plot grid for only one value of quantity3, to reduce clutter
        df_ratios = df_ratios.dropna()

        xarr = np.linspace(np.min(df_ratios[args.quantity1]),np.max(df_ratios[args.quantity1]), nbins)
        df_ratios['bin'] = pd.cut(df_ratios[args.quantity1], bins=xarr)
        grouped = df_ratios.groupby('bin')
        xbins = grouped[args.quantity1].mean().values
        ybins = np.log10(grouped[ratio_name].max().values)
        good_mask = np.isfinite(xbins) & np.isfinite(ybins)
        xbins = xbins[good_mask]
        ybins = ybins[good_mask]
        ax.plot(xbins, ybins, c='r', lw=0.5)

        col_arr = ['brown', 'darkgreen', 'cornflowerblue']
        for index, (p_init, func) in enumerate(zip(p_init_arr, func_arr)):
            popt, pcov = curve_fit(func, xbins, ybins, p0=p_init)
            ax.plot(xarr, func(xarr, *p_init), c=col_arr[index], lw=1, ls='--')
            ax.plot(xarr, func(xarr, *popt), c=col_arr[index], lw=2)
            ax.text(0.98, 0.01 + index * 0.07, f'fit{index+1} = [{",".join([f"{item:.2f}" for item in popt])}]', c=col_arr[index], ha='right', va='bottom', fontsize=args.fontsize, transform=ax.transAxes)

    # --------for talk plots--------------
    if args.fortalk:
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    # ---------saving figure--------------
    if args.phot_models.lower() in ['mappings', 'map']: phot_model_text = 'mappings'
    elif args.phot_models.lower() in ['nebulabayes', 'nb']: phot_model_text = 'NebulaBayes'
    slice_text = ''
    if args.slice_at_quantity2 is not None: slice_text += f'_{args.quantity2}={",".join([str(item) for item in args.slice_at_quantity2])}'.replace('/', '-')
    if args.slice_at_quantity3 is not None: slice_text += f'_{args.quantity3}={",".join([str(item) for item in args.slice_at_quantity3])}'.replace('/', '-')
    figname = args.mappings_dir / 'plots' / f'{phot_model_text}_model_{ratio_name.replace("/", "-")}_vs_{args.quantity1}{slice_text}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def number(s):
    '''
    Function to extract float numbers from MAPPINGS model files
    Inputs string, returns float
    '''
    if s[-1].isdigit():
        return float(s)
    else:
        return float(s[:-1])

# --------------------------------------------------------------------------------------------------------------------
def create_grid_file(df_ratio_names, outfilename, args):
    '''
    Creates 3D MAPPINGS grid of a given list of line ratios
    Writes the grid as an ASCII file into a specified outfilename
    '''
    print(f'Creating new MAPPINGS grid file..')
    # ----------------hard coded values-----------------------
    logq_arr = np.arange(6.5, 8.5 + 0.25, 0.25)
    metallicity_arr = np.array((0.05, 0.2, 0.4, 1., 2.))
    density_arr = np.arange(0, 5 + 0.5, 0.5)
    pressure_arr = np.arange(4, 9 + 0.5, 0.5)
    logOHSun = 8.93
    geom_path_dict = {'s':['spherical', 'sp'], 'p':['plane_par', 'pp']} # to determine directory structures based on geometry and iso parameters
    iso_label_dict = {'d': 'log(n_e)', 'P': 'log(P/k)'}

    # --------------------------------------------------------
    if args.iso == 'd': iso_arr = density_arr
    elif args.iso == 'P': iso_arr = pressure_arr
    model_dir_geom = args.mappings_dir / f'{args.iso}_{geom_path_dict[args.geometry][0]}'
    df_ratios = pd.DataFrame(columns=np.hstack(([iso_label_dict[args.iso], 'log(q)', 'Z'], df_ratio_names['numerator_label'] +'/'+ df_ratio_names['denominator_label'])))

    # --------------------------------------------------------
    for iso in iso_arr:
        model_dir_iso = model_dir_geom / f'{geom_path_dict[args.geometry][1]}_{args.iso}{iso * 10:.0f}_a05modelfiles'
        for logq in logq_arr:
            model_dir_q = model_dir_iso / f'Q{logq * 100:.0f}'
            for index, Z in enumerate(metallicity_arr):
                model_file_Z = model_dir_q / f'spec{index + 1:04d}.csv'
                print(f'Reading in MAPPINGS file {model_file_Z}..')
                fname = open(model_file_Z, 'r')
                lines = fname.readlines()
                ratios = []
                for this_num, this_den in zip(df_ratio_names['numerator_lambda_air'], df_ratio_names['denominator_lambda_air']):
                    numerator, denominator = 0., 0.
                    for line in lines:
                        if any([num in line for num in this_num.split(',')]):
                            numerator += number(line.split()[2])
                            continue
                        if any([den in line for den in this_den.split(',')]):
                            denominator += number(line.split()[2])
                            continue
                    if denominator == 0: ratios.append(np.inf) # to avoid division by zero error
                    else: ratios.append(numerator / denominator)

                df_ratios.loc[len(df_ratios)] = np.hstack(([iso, logq, np.log10(Z) + logOHSun], ratios))
                if 'fname' in locals(): fname.close()

    df_ratios.to_csv(outfilename, sep='\t', index=None)
    print(f'Saved ratio grid as file {outfilename}')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if args.fontsize == 10: args.fontsize = 15
    geom_path_dict = {'s':['spherical', 'sp'], 'p':['plane_par', 'pp']} # to determine directory structures based on geometry and iso parameters
    grid_dir = args.mappings_dir / 'grids'
    grid_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.mappings_dir / 'plots'
    plot_dir.mkdir(parents=True, exist_ok=True)
    if not args.keep: plt.close('all')

    # ---------making/reading in the grid-----------------------------
    if args.phot_models.lower() in ['mappings', 'map']:
        df_ratio_names = pd.read_table(args.mappings_dir / 'ratio_list.txt', delim_whitespace=True, skiprows=3)
        grid_filename = grid_dir / f'mappings_grid_{geom_path_dict[args.geometry][1]}_iso_{args.iso}.txt'
        if not grid_filename.exists() or args.clobber: create_grid_file(df_ratio_names, grid_filename, args)
        else: print(f'Reading in existing grid from {grid_filename}')
        df_ratios = pd.read_table(grid_filename, delim_whitespace=True)

    elif args.phot_models.lower() in ['nebulabayes', 'nb']:
        import NebulaBayes
        grid_filename = Path(NebulaBayes.__path__[0]) / 'grids' / 'NB_HII_grid.fits.gz'
        print(f'Reading in existing grid from {grid_filename}')
        df_grid = Table(fits.getdata(grid_filename)).to_pandas()
        df_grid['log q'] = np.round(df_grid['log U'] + np.log10(3e10), 1)

        quant_names_dict = {'Z':'12 + log O/H', 'log(q)':'log q', 'log(P/k)':'log P/k', 'log(U)':'log U'}
        line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OIII-4363': 'OIII4363', 'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'SII': 'SII6716_31', 'NeIII': 'NeIII3869'})

        quant_names = [quant_names_dict[quant] for quant in [args.quantity1, args.quantity2, args.quantity3]]
        df_ratios = df_grid[quant_names].sort_values(by=quant_names)
        df_ratios = df_ratios.rename(columns={v: k for k, v in quant_names_dict.items()})
        df_ratios[f'{line_label_dict[args.xnum_line]}/{line_label_dict[args.xden_line]}'] = df_grid[line_label_dict[args.xnum_line]] / df_grid[line_label_dict[args.xden_line]]
        df_ratios[f'{line_label_dict[args.ynum_line]}/{line_label_dict[args.yden_line]}'] = df_grid[line_label_dict[args.ynum_line]] / df_grid[line_label_dict[args.yden_line]]

    # --------plot the ratio grid---------------
    if args.plot_grid: fig_grid = plot_ratio_grid_fig(df_ratios, args)

    # --------plot the model ratios---------------
    if args.plot_model: fig_model = plot_ratio_model_fig(df_ratios, args)


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
