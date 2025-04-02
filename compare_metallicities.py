'''
    Filename: compare_metallicities.py
    Notes: Compares spatially resolved metallicity map across different diagnostics for a given set of galaxies
    Author : Ayan
    Created: 02-04-25
    Example: run compare_metallicities.py --field Par028 --id 300,1303,1634,1849,2171,2727,2867 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,P25,R3,NB --use_original_NB_grid
             run compare_metallicities.py --field Par028 --id 300 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,R3,NB
'''
from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -----figuring out the array fo diagnostics to use------------
    Zdiag_arr = args.Zdiag.split(',')
    if 'NB' in args.Zdiag and args.use_original_NB_grid: Zdiag_arr += ['NB_orig_grid']
    Zdiag_combinations = list(itertools.combinations(Zdiag_arr, 2))

    # ---------setting up filenames-----------------------------
    plots_dir = args.output_dir / f'{args.field}' / f'Zdiag_comparisons'
    plots_dir.mkdir(parents=True, exist_ok=True)

    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'

    # --------setting up full figure----------------------------------
    nrow, ncol = len(Zdiag_arr) - 1, len(Zdiag_arr) - 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(8, 7), layout='constrained')  # layout = 'tight' or 'constrained'

    col_arr = ['salmon', 'sienna', 'cornflowerblue', 'darkolivegreen', 'darkgoldenrod', 'darkorchid', 'darkcyan', 'hotpink']
    Zdiag_label_dict = smart_dict({'NB_orig_grid': 'NB (original grid)'})
    Z_limits = [7.5, 9.5]

    # --------looping over objects---------------------
    args.id_arr = np.atleast_1d(args.id)
    for index, args.id in enumerate(args.id_arr):
        print(f'\nDoing object {args.id} which is {index + 1} out of {len(args.id_arr)} objects..')

        df = pd.DataFrame()
        log_int_array = []
        # --------looping over diagnostics---------------------
        for index2, Zdiag in enumerate(Zdiag_arr):
            if Zdiag == 'NB_orig_grid':
                fitsname = args.output_dir / 'catalogs' / f'{args.field}_{args.id:05d}_logOH_map{snr_text}{only_seg_text}{vorbin_text}_Zdiag_NB_AGNdiag_{args.AGN_diag}_orig_grid.fits'
            else:
                fitsname = args.output_dir / 'catalogs' / f'{args.field}_{args.id:05d}_logOH_map{snr_text}{only_seg_text}{vorbin_text}_Zdiag_{Zdiag}_AGNdiag_{args.AGN_diag}.fits'
            if os.path.exists(fitsname):
                hdul = fits.open(fitsname)
                this_df = Table(hdul['tab'].data).to_pandas()
                if index2 == 0:
                    df = this_df
                else:
                    df = pd.merge(df, this_df, on = ['radius', 'bin_ID', 'agn_dist'])
                df = df.rename(columns={'log_OH': f'log_OH_{Zdiag}', 'log_OH_u': f'log_OH_u_{Zdiag}'})

                # ----getting the integrated values-------
                header = hdul['log_OH'].header
                this_logOH_int = ufloat(header['LOG_OH_INT'], header['LOG_OH_INT_ERR']) if header['LOG_OH_INT'] is not None else ufloat(np.nan, np.nan)
                log_int_array.append(this_logOH_int)

            else:
                log_int_array.append(this_logOH_int)
                print(f'{fitsname} not found, so skipping this diagnostic {Zdiag}')

        # ------now plotting for every diag combination--------------
        for col_index in range(ncol):
            for row_index in range(nrow):
                ax = axes[row_index][col_index]
                if row_index < col_index:
                    try: fig.delaxes(ax)
                    except: pass
                else:
                    Z1_index, Z2_index = col_index, row_index + 1
                    Zdiag1 = Zdiag_arr[Z1_index]
                    Zdiag2 = Zdiag_arr[Z2_index]
                    color = df['bin_ID'] # col_arr[index] # col_arr[index OR df['bin_ID'] OR df['radius]

                    if f'log_OH_{Zdiag1}' in df and f'log_OH_{Zdiag2}' in df:
                        # ---plotting integrated--
                        ax.scatter(log_int_array[Z1_index].n, log_int_array[Z2_index].n, s=20, c=col_arr[index], lw=0.5, edgecolor='k')
                        ax.errorbar(log_int_array[Z1_index].n, log_int_array[Z2_index].n, xerr=log_int_array[Z1_index].s, yerr=log_int_array[Z2_index].s, c=col_arr[index], fmt='none', lw=1, alpha=0.5, zorder=-10)

                        # ---plotting spatially resolved--
                        ax.scatter(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], s=5, c=color, lw=0, cmap='viridis')
                        ax.errorbar(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], xerr=df[f'log_OH_u_{Zdiag1}'], yerr=df[f'log_OH_u_{Zdiag2}'], c='grey', fmt='none', lw=0.5, alpha=0.1, zorder=-10)

                        # ----annotate axis----------
                        if index == len(args.id_arr) - 1:
                            ax.plot(Z_limits, Z_limits, ls='--', lw=0.1, c='k')

                            ax.set_xlim(Z_limits)
                            ax.set_ylim(Z_limits)

                            if Z1_index == 0:
                                ax.set_ylabel(Zdiag_label_dict[Zdiag2], fontsize=args.fontsize)
                                ax.tick_params(axis='y', which='major', labelsize=args.fontsize)
                            else:
                                ax.set_yticklabels([])

                            if Z2_index == nrow:
                                ax.set_xlabel(Zdiag_label_dict[Zdiag1], fontsize=args.fontsize)
                                ax.tick_params(axis='x', which='major', labelsize=args.fontsize)
                            else:
                                ax.set_xticklabels([])
                    else:
                        print(f'Cannot plot {Zdiag1} vs {Zdiag2} comparison for object {args.id}')

        fig.text(0.98, 0.98 - index * 0.03, f'ID #{args.id}', c=col_arr[index], ha='right', va='top', fontsize=args.fontsize)

    # ------------saving the full figure--------------------------
    figname = plots_dir / f'{",".join(args.id_arr.astype(str))}Zdiag_comparison_{",".join(Zdiag_arr)}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
