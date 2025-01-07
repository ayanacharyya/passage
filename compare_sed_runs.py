'''
    Filename: compare_sed_runs.py
    Notes: Plots various quantities from 2 different dataframes corresponding to 2 different runs of SED fitting, on the same data (produced by compute_stellar_masses.py)
    Author : Ayan
    Created: 20-12-24
    Example: run compare_sed_runs.py --plot_conditions EW,mass,PA --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z_narrow_mass,only_st_bands
             run compare_sed_runs.py --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z_narrow_mass,only_st_bands
'''

from header import *
from util import *
from make_passage_plots import bounds_dict, label_dict, colormap_dict

start_time = datetime.now()

# -------------------------------------global dictionaries-------------------------------------------------------------------------------
bounds_dict.update({'log_mass_bgp': (6.5, 10.5)})
run_labels_dict = {'narrow_z': 'All good filters', 'narrow_z_narrow_mass': 'All good filters', 'only_st_bands': 'Only ACS+NIRISS', 'increased_flux_err':' 30% more flux err'}

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -----------declaring fig handles--------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    # -------determining SED runs and file names----------------
    runs = args.run.split(',')
    suffix_dict = defaultdict(lambda: '', _x=f': {run_labels_dict[runs[0]]}', _y=f': {run_labels_dict[runs[1]]}')
    figname = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df_{args.xcol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}_vs_{args.ycol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}_colorby_{args.colorcol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}.png'

    # -------reading in and merging dataframe produced by compute_stellar_masses.py----------------
    df_infilename = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df.txt'
    df_infilename_x = Path(str(df_infilename).replace('.txt', f'_withSED_{runs[0]}.csv'))
    df_infilename_y = Path(str(df_infilename).replace('.txt', f'_withSED_{runs[1]}.csv'))

    df_x = pd.read_csv(df_infilename_x)
    df_y = pd.read_csv(df_infilename_y)
    print(f'Read in dfs to compare from {df_infilename_x} and {df_infilename_y}')

    bgp_columns = [item for item in df_x.columns if 'bgp' in item]
    other_columns = list(set(df_x.columns) - set(bgp_columns))
    df = df_x.merge(df_y, on=other_columns)

    # ---------making the plot-----------
    p = ax.scatter(df[args.xcol], df[args.ycol], c=df[args.colorcol], plotnonfinite=True, s=100, lw=1, edgecolor='w' if args.fortalk else 'k', vmin=bounds_dict[args.colorcol.rstrip('_x').rstrip('_y')][0] if args.colorcol.rstrip('_x').rstrip('_y') in bounds_dict else None, vmax=bounds_dict[args.colorcol.rstrip('_x').rstrip('_y')][1] if args.colorcol.rstrip('_x').rstrip('_y') in bounds_dict else None)
    if args.ycol + '_u' in df and not args.foggie_comp:  # if uncertainty column exists
        ax.errorbar(df[args.xcol], df[args.ycol], yerr=df[args.ycol + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
    if args.xcol + '_u' in df and not args.foggie_comp:  # if uncertainty column exists
        ax.errorbar(df[args.xcol], df[args.ycol], xerr=df[args.xcol + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)

    cbar = plt.colorbar(p)
    cbar.set_label(label_dict[args.colorcol.rstrip('_x').rstrip('_y')] + suffix_dict[args.colorcol[-2:]] if args.colorcol.rstrip('_x').rstrip('_y') in label_dict else args.colorcol)

    # ---------annotate axes and save figure-------
    if len(ax.get_legend_handles_labels()[0]) > 0: plt.legend()
    ax.set_xlabel(label_dict[args.xcol.rstrip('_x').rstrip('_y')] + suffix_dict[args.xcol[-2:]] if args.xcol.rstrip('_x').rstrip('_y') in label_dict else args.xcol)
    ax.set_ylabel(label_dict[args.ycol.rstrip('_x').rstrip('_y')] + suffix_dict[args.ycol[-2:]] if args.ycol.rstrip('_x').rstrip('_y') in label_dict else args.ycol)

    if args.xcol.rstrip('_x').rstrip('_y') in bounds_dict: ax.set_xlim(bounds_dict[args.xcol.rstrip('_x').rstrip('_y')][0], bounds_dict[args.xcol.rstrip('_x').rstrip('_y')][1])
    if args.ycol.rstrip('_x').rstrip('_y') in bounds_dict: ax.set_ylim(bounds_dict[args.ycol.rstrip('_x').rstrip('_y')][0], bounds_dict[args.ycol.rstrip('_x').rstrip('_y')][1])

    # ---------comparing SFRs-------
    if args.xcol[:-2] == args.ycol[:-2]:
        line = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
        ax.plot(line, line, ls='dashed', c='k', lw=1)

    # --------for talk plots--------------
    if args.fortalk:
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
