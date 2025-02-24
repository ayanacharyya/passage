'''
    Filename: compare_sed_runs.py
    Notes: Plots various quantities from 2 different dataframes corresponding to 2 different runs of SED fitting, on the same data (produced by compute_stellar_masses.py)
    Author : Ayan
    Created: 20-12-24
    Example: run compare_sed_runs.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z_narrow_mass,only_st_bands
             run compare_sed_runs.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z_narrow_mass,only_st_bands

             run compare_sed_runs.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z,all_st_bands
             run compare_sed_runs.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp_x --ycol log_mass_bgp_y --colorcol redshift --run narrow_z,all_st_bands

             run compare_sed_runs.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --drv 0.1 --run all_ground_based,all_space_based --id 3139
             run compare_sed_runs.py --line_list OIII,Ha --plot_conditions EW,mass,PA --drv 0.1 --run all_ground_based,all_space_based --do_all_obj
'''

from header import *
from util import *
import h5py
from make_passage_plots import bounds_dict, label_dict

start_time = datetime.now()

# -------------------------------------global dictionaries-------------------------------------------------------------------------------
bounds_dict.update({'log_mass_bgp': (6.5, 10.5)})
run_labels_dict = {'Par028_including_nircam':'COSMOS+PASSAGE+COSMOS-Webb', 'Par028_all_st_bands':'ACS+NIRISS+NIRCAM+MIRI', 'Par052_including_nircam':'COSMOS+PASSAGE+COSMOS-Webb', 'Par052_all_st_bands':'ACS+NIRISS+NIRCAM+MIRI', 'all_st_bands':'ACS+NIRISS+NIRCAM+MIRI', 'all_ground_based':'CFHT+UVISTA+HSC+SC', 'all_space_based':'GALEX+IRAC+ACS+NIRISS+NIRCAM+MIRI', 'including_nircam':'COSMOS+PASSAGE+COSMOS Webb', 'narrow_z': 'COSMOS+PASSAGE', 'narrow_z_narrow_mass': 'All good filters', 'only_st_bands': 'ACS+NIRISS', 'increased_flux_err':' 30% more flux err'}

# ----------------------------------------------------------------------------
def display_file_in_ax(filename, ax):
    '''
    Loads the given filename and displays it in the given axis handle
    Returns axis handle
    '''
    pdf_file = fitz.open(filename)
    rgb = pdf_file[0].get_pixmap()
    pil_image = Image.open(io.BytesIO(rgb.tobytes()))
    ax.imshow(pil_image.convert('RGB'), origin='upper')

    # sfh = mpimg.imread(filename)
    # ax.imshow(sfh, origin='upper')

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -------determining SED runs and file names----------------
    plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
    plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
    runs = args.run.split(',')

    # --------------compare individual object properties---------------------
    if args.id is not None or args.do_all_obj:
        all_ids = np.sort([int(os.path.split(item)[-1][:-3]) for item in glob.glob(str(args.output_dir / f'pipes/posterior/{runs[0]}/*.h5'))])
        if args.do_all_obj: args.id_arr = all_ids
        else: args.id_arr = args.id

        output_dir = args.output_dir / f'pipes/plots/comparisons_{",".join(runs)}/'
        output_dir.mkdir(parents=True, exist_ok=True)

        for index, args.id in enumerate(args.id_arr):
            print(f'\nCommencing ID {args.id} which is {index+1} of {len(args.id_arr)}..')
            # -------compare individual SED fits------------------
            fig, axes = plt.subplots(len(runs), 2, figsize=(12, 1 + 2 * len(runs)))
            fig.subplots_adjust(left=0.04, right=0.99, bottom=0.01, top=0.99, hspace=0.02, wspace=0.02)
            figname = output_dir / f'allpar_{args.drv}_{args.id}_SED_SFH_comp_{",".join(runs)}.png'

            plot_paths = [args.output_dir / 'pipes/plots' / item for item in runs]
            posterior_paths = [args.output_dir / 'pipes/posterior' / item for item in runs]

            for index in range(len(runs)):
                axes[index][0] = display_file_in_ax(plot_paths[index] / f'{args.id}_fit_log_x.pdf', axes[index][0])
                axes[index][1] = display_file_in_ax(plot_paths[index] / f'{args.id}_sfh.pdf', axes[index][1])

                posterior_file = h5py.File(posterior_paths[index] / f'{args.id}.h5', 'r')
                axes[index][1].text(0.9, 0.8, f'log M* = {posterior_file["median"][4]: .02f}', c='k', fontsize=args.fontsize, ha='right', va='top', transform=axes[index][1].transAxes)
                fig.text(0.02, 0.5 / len(runs) + (1 - (index + 1) / len(runs)), f'{run_labels_dict[runs[index]] if runs[index] in run_labels_dict else runs[index]}', c='k', fontsize=args.fontsize/1.3, rotation=90, ha='center', va='center')

            if len(args.id_arr) > 20: plt.close()
            # --------for talk plots--------------
            if args.fortalk:
                try:
                    mplcyberpunk.make_scatter_glow()
                except:
                    pass

            fig.savefig(figname, transparent=args.fortalk)
            print(f'\nSaved figure as {figname}')
            plt.show(block=False)

    # --------------compare population properties---------------------
    else:
        # -----------declaring fig handles--------
        fig, ax = plt.subplots(1, figsize=(8, 6))
        fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

        # -------determining SED runs and file names----------------
        try: suffix_dict = defaultdict(lambda: '', _x=f': {run_labels_dict[runs[0]]}', _y=f': {run_labels_dict[runs[1]]}')
        except KeyError: suffix_dict = defaultdict(lambda: '', _x=f': {runs[0]}', _y=f': {runs[1]}')
        figname = args.output_dir / 'plots' / f'allpar_venn_{plot_conditions_text}_SEDcomp_{args.xcol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}_vs_{args.ycol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}_colorby_{args.colorcol.replace("_x", "_" + runs[0]).replace("_y", "_" + runs[1])}.png'

        # -------reading in and merging dataframe produced by compute_stellar_masses.py----------------
        if args.do_field is None:
            plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
            plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
            args.field_set_plot_conditions_text = f'allpar_{args.drv}_venn_{plot_conditions_text}'
            df_infilename = args.output_dir / 'catalogs' / f'{args.field_set_plot_conditions_text}_df.txt'
        else:
            args.field = f'Par{int(args.do_field.split("Par")[1]):03d}'
            args.field_set_plot_conditions_text = f'{args.field}_{args.drv}_allmatch'
            df_infilename = args.output_dir / args.field / f'{args.field}_all_diag_results.txt'

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
            if 'mass_bgp' in args.xcol and 'mass_bgp' in args.ycol: ax.plot(line, line - np.median(df[args.xcol].values - df[args.ycol].values), ls='dashed', c='r', lw=1)

        # --------for talk plots--------------
        if args.fortalk:
            try: mplcyberpunk.make_scatter_glow()
            except: pass

        fig.savefig(figname, transparent=args.fortalk)
        print(f'\nSaved figure as {figname}')
        plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
