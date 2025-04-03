'''
    Filename: compare_metallicities.py
    Notes: Compares spatially resolved metallicity map across different diagnostics for a given set of galaxies
    Author : Ayan
    Created: 02-04-25
    Example: run compare_metallicities.py --field Par028 --id 300,1303,1634,1849,2171,2727,2867 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,P25,R3,NB --use_original_NB_grid
             run compare_metallicities.py --field Par028 --id 300 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,R3,NB --colorcol radius --debug_Zdiag
'''
from header import *
from util import *
from make_diagnostic_maps import plot_2D_map, plot_radial_profile

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.colorcol == 'ez_z_phot': args.colorcol = 'color'

    # -----figuring out the array fo diagnostics to use------------
    Zdiag_arr = args.Zdiag.split(',')
    if 'NB' in args.Zdiag and args.use_original_NB_grid: Zdiag_arr += ['NB_orig_grid']
    args.id_arr = np.atleast_1d(args.id)
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)

    # ---------setting up filenames-----------------------------
    plots_dir = args.output_dir / f'{args.field}' / f'Zdiag_comparisons'
    plots_dir.mkdir(parents=True, exist_ok=True)

    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'

    # --------setting up full figure----------------------------------
    nrow, ncol = len(Zdiag_arr) - 1, len(Zdiag_arr) - 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(8, 7), layout='tight')  # layout = 'tight' or 'constrained'
    axes = np.atleast_2d(axes)

    col_arr = ['salmon', 'sienna', 'cornflowerblue', 'darkolivegreen', 'darkgoldenrod', 'darkorchid', 'darkcyan', 'hotpink']
    Zdiag_label_dict = smart_dict({'NB_orig_grid': 'NB (original grid)'})
    Z_limits = [7.1, 9.1]
    color_lim_dict = {'color':[None, None, '', ''], 'bin_ID':[None, None, 'Voronoi bin ID', 'rainbow'], 'radius':[0, 5, 'Galactocentric distance (kpc)', 'cividis'], 'agn_dist':[-1, 1, f'Distance from {args.AGN_diag} SF line', args.diverging_cmap]}

    args.radius_max = color_lim_dict['radius'][1]
    offset = 2 if args.plot_radial_profiles else 1

    # -------setting up diagnostic figure----------
    if args.debug_Zdiag:
        nrow2, ncol2 = len(args.id_arr), len(Zdiag_arr) + 1
        if args.plot_radial_profiles: nrow2 *= 2
        fig2, axes2 = plt.subplots(nrow2, ncol2, figsize=(8, 7), layout='tight')
        axes2 = np.atleast_2d(axes2)

    # --------looping over objects---------------------
    for index, args.id in enumerate(args.id_arr):
        print(f'\nDoing object {args.id} which is {index + 1} out of {len(args.id_arr)} objects..')
        markersize = (400 - index * 300 / len(args.id_arr)) / len(Zdiag_arr)
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
                df = df.rename(columns={'log_OH': f'log_OH_{Zdiag}', 'log_OH_u': f'log_OH_{Zdiag}_err'})
                df['color'] = col_arr[index]

                # ----getting the integrated values-------
                header = hdul['log_OH'].header
                if header['LOG_OH_INT'] is None:
                    print(f'Integrated logOH not found for ID {args.id} diag {Zdiag}')
                    this_logOH_int = ufloat(np.nan, np.nan)
                else:
                    this_logOH_int = ufloat(header['LOG_OH_INT'], header['LOG_OH_INT_ERR'])
                log_int_array.append(this_logOH_int)

                # -------plotting metallicity map----------
                if args.debug_Zdiag:
                    logOH_map = hdul['log_OH'].data
                    axes2[index * offset][index2], _ = plot_2D_map(logOH_map, axes2[index * offset][index2], args, takelog=False, label=f'{Zdiag}: Z = {this_logOH_int.n:.1f}', cmap='viridis', hide_yaxis=index2 > 0, hide_xaxis=index < len(args.id_arr) - 1, vmin=Z_limits[0], vmax=Z_limits[1], hide_cbar=index2 < len(Zdiag_arr) - 1)
                    if args.plot_radial_profiles: axes2[index * offset + 1][index2], _ = plot_radial_profile(df, axes2[index * offset + 1][index2], args, label=f'log(O/H) + 12', ymin=Z_limits[0], ymax=Z_limits[1], hide_xaxis=index < len(args.id_arr) - 1, hide_yaxis=index2 > 0, xcol='radius', ycol=f'log_OH_{Zdiag}')

            else:
                log_int_array.append(this_logOH_int)
                print(f'{fitsname} not found, so skipping this diagnostic {Zdiag}')

        # -------plotting voronoi bin map----------
        if args.debug_Zdiag:
            bin_ID_map = hdul['bin_ID'].data
            axes2[index * offset][index2 + 1], _ = plot_2D_map(bin_ID_map, axes2[index * offset][index2 + 1], args, takelog=False, label=f'ID {args.id}: Vorbin ID', cmap='rainbow', hide_yaxis=True, hide_xaxis=index < len(args.id_arr) - 1, hide_cbar=False)
            if args.plot_radial_profiles: fig2.delaxes(axes2[index * offset + 1][index2 + 1])

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

                    if f'log_OH_{Zdiag1}' in df and f'log_OH_{Zdiag2}' in df:
                        # ---plotting integrated--
                        ax.scatter(log_int_array[Z1_index].n, log_int_array[Z2_index].n, s=markersize, c=col_arr[index], lw=0.5, edgecolor='k', marker='s')
                        ax.errorbar(log_int_array[Z1_index].n, log_int_array[Z2_index].n, xerr=log_int_array[Z1_index].s, yerr=log_int_array[Z2_index].s, c=col_arr[index], fmt='none', lw=0.5, alpha=0.5)

                        # ---plotting spatially resolved--
                        p = ax.scatter(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], s=markersize/2, c=df[args.colorcol], lw=0, cmap=color_lim_dict[args.colorcol][3], vmin=color_lim_dict[args.colorcol][0], vmax=color_lim_dict[args.colorcol][1])
                        ax.errorbar(df[f'log_OH_{Zdiag1}'], df[f'log_OH_{Zdiag2}'], xerr=df[f'log_OH_{Zdiag1}_err'], yerr=df[f'log_OH_{Zdiag2}_err'], c='grey', fmt='none', lw=0.1, alpha=0.5)

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

        # ------annotating figure----------------------
        fig.text(0.88, 0.96 - index * 0.03, f'ID #{args.id}', c=col_arr[index], ha='left', va='top', fontsize=args.fontsize)

        # ------making colorbar----------------------
        if args.colorcol != 'color' and index == len(args.id_arr) - 1:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.89, 0.3, 0.02, 0.4]) # left, bottom, width, height
            cbar = fig.colorbar(p, cax=cbar_ax, orientation='vertical')
            cbar.set_label(color_lim_dict[args.colorcol][2], fontsize=args.fontsize)

    # -------save diagnostic figure----------
    if args.debug_Zdiag:
        figname2 = plots_dir / f'{",".join(args.id_arr.astype(str))}Zdiag_maps_{",".join(Zdiag_arr)}.png'
        fig2.savefig(figname2, transparent=args.fortalk)
        print(f'Saved figure as {figname2}')

    # ------------saving the full figure--------------------------
    colorby_text = f'_colorby_{args.colorcol}' if args.colorcol != 'color' else ''
    figname = plots_dir / f'{",".join(args.id_arr.astype(str))}Zdiag_comparison_{",".join(Zdiag_arr)}{colorby_text}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
