'''
    Filename: compare_stacked_line_ratios.py
    Notes: Compared integraed emission line flux ratios between all available 1D vs 2D stacks 
    Author : Ayan
    Created: 16-11-26
    Example: run compare_stacked_line_ratios.py --system ssd --do_all_fields --cut_z_flag 4 --adaptive_bins --bin_by_distance_mass --skip_deproject --include_par682 --include_hst --inc_cut 30 --bin_spatial --nb_mass_prior --phot nb --slice_at_quantity2 8.5
'''

from header import *
from util import *
from make_sfms_bins import get_binned_df, required_lines, sfms
from plot_stacked_maps import read_stacked_maps, get_emission_line_maps
from plot_mappings_grid import plot_ratio_model

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.5

    # ------------reading and binning dataframe-------------
    df_binned, bin_list, args = get_binned_df(args, required_lines=required_lines, sfms=sfms)

    df_2dratios_filename = args.fits_dir / f'stacked_{args.binby_text}{args.binned_text}{args.fold_text}_ratios_all_2dstack_{args.deproject_text}{args.rescale_text}.csv'
    df_1dratios_filename = args.output_dir / 'catalogs' / f'ratios_all_1dstack.csv' # args.output_dir / 'catalogs' / 'stacked_spec_ms' / 

    # -------making 1D stacked ratios list-----------
    if not os.path.exists(df_1dratios_filename) or args.clobber:
        print(f'Making list of 1D stacked ratios..')

        line_list = ['OII', 'OIIIboth', 'Hb']

        # -----dict for converting line label names to those acceptable to NB---------
        mylabel_to_stacklabel_dict = {'OII':'OII', 'OIII':'OIIIboth', 'Hb':'Hb', 'Ha':'HaNII', 'SII':'SII'}

        # --------finding files---------------
        files = glob.glob(str(args.output_dir / 'catalogs/stacked_spec_ms') + '/meas_logm*.dat')
        df_1dratios = pd.DataFrame(columns=['log_mass_min', 'log_mass_max', 'log_mass_median', 'R23_int', 'R23_int_u', 'R3_int', 'R3_int_u'])

        # ------looping over stacks-----------
        for index, file in enumerate(files):
            df = pd.read_csv(file, sep=r'\s*\|\s*', engine='python').iloc[:, 1:-1]
            log_mass_min, log_mass_max = map(float, Path(file).stem.split('_')[-1].split('-'))
            log_mass_mean = np.mean(np.array([log_mass_min, log_mass_max]))
            print(f'\nStarting ({index + 1}/{len(files)}) stack {log_mass_min} <= log_mass <= {log_mass_max}..')

            hb = df[df['line'] == 'Hb']
            hb_flux = ufloat(hb['flux_norm'].values[0], hb['flux_norm_err'].values[0])
            oii = df[df['line'] == 'OII']
            oii_flux = ufloat(oii['flux_norm'].values[0], oii['flux_norm_err'].values[0])
            oiii = df[df['line'] == 'OIIIboth']
            oiii_both_flux = ufloat(oiii['flux_norm'].values[0], oiii['flux_norm_err'].values[0])

            # ----------deblending flux--------------------
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            print(f'\tCorrecting OIII for 4959 component, by factor of {factor:.3f}')
            oiii_flux = oiii_both_flux * factor

            # --------determine ratios---------------
            R3_int = unp.log10(oiii_flux / hb_flux)
            R23_int = unp.log10((oii_flux + oiii_both_flux) / hb_flux)

            # ------appending to master df----------
            thisrow = {'log_mass_min': log_mass_min, 'log_mass_max': log_mass_max, 'log_mass_median': log_mass_mean, \
                       'R3_int': R3_int.tolist().n, \
                       'R3_int_u': R3_int.tolist().s, \
                       'R23_int': R23_int.tolist().n, \
                       'R23_int_u': R23_int.tolist().s, \
                        }
            df_1dratios.loc[len(df_1dratios)] = thisrow

        # ----------saving master dataframe-----------
        df_1dratios.to_csv(df_1dratios_filename, index=None)
        print(f'Saved results as {df_1dratios_filename}')
    else:
        print(f'Reading 1D ratios from existing file {df_1dratios_filename}')
        df_1dratios = pd.read_csv(df_1dratios_filename)

    # -------making 2D stacked ratios list-----------
    if not os.path.exists(df_2dratios_filename) or args.clobber:
        print(f'Making list of 2D stacked ratios..')

        # ------------setting up the dataframe----------------------
        if args.bin_by_distance_mass:
            df2 = df_binned.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'bin_intervals':'delta_sfms_bin'})
            df2 = df2.groupby(['delta_sfms_bin', 'log_mass_bin']).agg(log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
            df2['delta_sfms_min'] = pd.IntervalIndex(df2['delta_sfms_bin']).left
            df2['delta_sfms_max'] = pd.IntervalIndex(df2['delta_sfms_bin']).right
            df2.drop(['delta_sfms_bin'], axis=1, inplace=True)
        elif args.bin_by_distance:
            df2 = df_binned.rename(columns={'bin_intervals':'delta_sfms_bin'})
            df2 = df2.groupby(['delta_sfms_bin']).agg( log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
            df2['delta_sfms_min'] = pd.IntervalIndex(df2['delta_sfms_bin']).left
            df2['delta_sfms_max'] = pd.IntervalIndex(df2['delta_sfms_bin']).right
            df2.drop(['delta_sfms_bin'], axis=1, inplace=True)
        elif args.bin_by_sfh_mass:
            df2 = df_binned.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'bin_intervals':'tform_ratio_bin'})
            df2 = df2.groupby(['tform_ratio_bin', 'log_mass_bin']).agg(log_mass_median=('log_mass', 'median'), tform_ratio_median=('delta_tform_ratio', 'median')).reset_index()
            df2['tform_ratio_min'] = pd.IntervalIndex(df2['tform_ratio_bin']).left
            df2['tform_ratio_max'] = pd.IntervalIndex(df2['tform_ratio_bin']).right
            df2.drop(['tform_ratio_bin'], axis=1, inplace=True)
        else:
            df2 = df_binned.rename(columns={'mass_interval':'log_mass_bin', 'mass_intervals':'log_mass_bin', 'sfr_intervals':'log_sfr_bin', 'sfr_interval':'log_sfr_bin'})
            df2 = df2.groupby(['log_mass_bin', 'log_sfr_bin']).agg(log_mass_median=('log_mass', 'median'), delta_sfms_median=('delta_sfms', 'median')).reset_index()
            df2['log_sfr_min'] = pd.IntervalIndex(df2['log_sfr_bin']).left
            df2['log_sfr_max'] = pd.IntervalIndex(df2['log_sfr_bin']).right
            df2.drop(['log_sfr_bin'], axis=1, inplace=True)

        df2['log_mass_min'] = pd.IntervalIndex(df2['log_mass_bin']).left
        df2['log_mass_max'] = pd.IntervalIndex(df2['log_mass_bin']).right
        df_2dratios = df2.drop(['log_mass_bin'], axis=1)
        common_cols = ['nobj', 'R23_int', 'R23_int_u', 'R3_int', 'R3_int_u']
        for col in common_cols:
            df_2dratios[col] = np.nan
        
        # ------------looping over each bin-----------------------
        for index, obj in df_2dratios.iterrows():
            if args.bin_by_distance:
                args.bin_text = f'delta_sfms_bin_{obj["delta_sfms_min"]}-{obj["delta_sfms_max"]}'
            elif args.bin_by_distance_mass:
                args.bin_text = f'delta_sfms_bin_{obj["delta_sfms_min"]}-{obj["delta_sfms_max"]}_mass_bin_{obj["log_mass_min"]}-{obj["log_mass_max"]}'
            elif args.bin_by_sfh_mass:
                args.bin_text = f'tform_ratio_bin_{obj["tform_ratio_min"]}-{obj["tform_ratio_max"]}_mass_bin_{obj["log_mass_min"]}-{obj["log_mass_max"]}'
            else:
                args.bin_text = f'logmassbin_{obj["log_mass_min"]}-{obj["log_mass_max"]}_logsfrbin_{obj["log_sfr_min"]}-{obj["log_sfr_max"]}'
            print(f'\n\tStarting ({index + 1}/{len(bin_list)}) {args.bin_text}..', end=' ')
            
            # -------reading previously saved stacked fits file------------
            stacked_filename = args.fits_dir / f'stacked_maps{args.deproject_text}{args.rescale_text}_{args.bin_text}.fits'
            if not stacked_filename.exists():
                print(f'No stacked fits file found for {args.bin_text}, so skipping this bin.')
                continue
            line_dict = read_stacked_maps(stacked_filename, args)
            if not line_dict:
                print(f'No lines found for {args.bin_text}. So Skipping.')
                continue

            # -----extract integrated line ratios------------
            _, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'OII', 'Hb'], args, silent=args.only_integrated)
            R3_int = unp.log10(line_int_arr[0] / line_int_arr[2])
            R23_int = unp.log10((line_int_arr[0] + line_int_arr[1]) / line_int_arr[2])

            # --------appending to df-----------            
            thisrow = {'nobj': nobj, 'R3_int': R3_int.tolist().n, 'R3_int_u': R3_int.tolist().s, 'R23_int': R23_int.tolist().n, 'R23_int_u': R23_int.tolist().s}
            df_2dratios.loc[index, list(thisrow.keys())] = list(thisrow.values())

        # ----------saving master dataframe-----------
        df_2dratios.to_csv(df_2dratios_filename, index=None)
        print(f'Saved results as {df_2dratios_filename}')
    else:
        print(f'Reading 2D ratios from existing file {df_2dratios_filename}')
        df_2dratios = pd.read_csv(df_2dratios_filename)

    # ---------set up the figure-------------
    xcol, colorcol = 'log_mass_median', 'delta_sfms_median'
    ycol_arr = ['R23_int']#, 'R3_int']

    nrow = len(ycol_arr)
    fig, axes = plt.subplots(nrow, 2, figsize=(9, 2 + 2 * nrow))
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(left=0.08, right=0.88, top=0.98, bottom=0.14 if nrow == 1 else 0.1, wspace=0.1, hspace=0.1)
    cmin, cmax, cmap = -0.6, 0.6, 'RdBu'

    label_dict = {
                'delta_sfms_median': r'$\delta_{SFMS}$ [$\log$ M$_\odot$yr$^{-1}$]',\
                'log_mass_median': r'Median $\log$ (M/M$_{\odot}$) of stack',\
                #'R3_int': r'R3 = $\log$ OIII 5007/H$\beta$',\
                #'R23_int': r'R23 = $\log$ (OII 3727 + OIII 4959,5007)/H$\beta$',\
                'R3_int': r'R3',\
                'R23_int': r'R23',\
                'Z': r'$log$ (O/H) + 12', \
                }
    
    # --------plotting observations-----------
    color = df_2dratios[colorcol] if colorcol is not None else 'salmon'
    for index, ycol in enumerate(ycol_arr):
        ax = axes[index][1]
        ax.scatter(df_1dratios[xcol], df_1dratios[ycol], s=100, c='yellow', lw=1, edgecolors='brown', marker='D', label='1D stacked')
        ax.errorbar(df_1dratios[xcol], df_1dratios[ycol], yerr=df_1dratios[f'{ycol}_u'], c='brown', lw=0.7, fmt='none', alpha=1)

        p = ax.scatter(df_2dratios[xcol], df_2dratios[ycol], s=100, c=color, lw=1, vmin=cmin, vmax=cmax, edgecolors='k', cmap=cmap, marker='o', label='2D stacked')
        ax.errorbar(df_2dratios[xcol], df_2dratios[ycol], yerr=df_2dratios[f'{ycol}_u'], c='grey', lw=0.7, fmt='none', alpha=1)

        ax = annotate_axes(ax, label_dict[xcol], label_dict[ycol], ylim=[0, 1.1], args=args, hide_yaxis=True,
                                    hide_xaxis=index < len(ycol_arr) - 1, clabel=label_dict[colorcol] if colorcol is not None else '', 
                                    hide_cbar=colorcol is None, cbar_width=2, p=p)
    ax.legend(fontsize = args.fontsize / args.fontfactor, loc='best')

    # ---------making/reading in the grid-----------------------------
    if args.phot_models.lower() in ['mappings', 'map']:
        geom_path_dict = {'s':['spherical', 'sp'], 'p':['plane_par', 'pp']} # to determine directory structures based on geometry and iso parameters
        grid_dir = args.mappings_dir / 'grids'
        df_ratio_names = pd.read_table(args.mappings_dir / 'ratio_list.txt', delim_whitespace=True, skiprows=3)
        grid_filename = grid_dir / f'mappings_grid_{geom_path_dict[args.geometry][1]}_iso_{args.iso}.txt'
        df_ratios = pd.read_table(grid_filename, delim_whitespace=True)

    elif args.phot_models.lower() in ['nebulabayes', 'nb', 'neb']:
        import NebulaBayes
        grid_filename = Path(NebulaBayes.__path__[0]) / 'grids' / 'NB_HII_grid.fits.gz'
        print(f'Reading in existing grid from {grid_filename}')
        df_grid = Table(fits.getdata(grid_filename)).to_pandas()
        df_grid['log q'] = np.round(df_grid['log U'] + np.log10(3e10), 1)

        quant_names_dict = {'Z':'12 + log O/H', 'log(q)':'log q', 'log(P/k)':'log P/k', 'log(U)':'log U'}
        line_label_dict = smart_dict({'OII': 'OII3726_29', 'Hb': 'Hbeta', 'OIII': 'OIII5007', 'OII,OIII':'OII3726_29,OIII5007', 'OIII-4363': 'OIII4363', 
                                      'OI-6302': 'OI6300', 'Ha': 'Halpha', 'NII':'NII6583', 'Ha,NII': 'Halpha,NII6583', 
                                      'SII': 'SII6716_31', 'NeIII': 'NeIII3869', 'NeIII,HeI': 'NeIII3869_HeI3889',
                                      'NeIII155000':'[NeIII]155000', 'NeII128000':'[NeII]128000'})

        quant_names = [quant_names_dict[quant] for quant in [args.quantity1, args.quantity2, args.quantity3]]
        df_ratios = df_grid[quant_names].sort_values(by=quant_names)
        df_ratios = df_ratios.rename(columns={v: k for k, v in quant_names_dict.items()})
        df_ratios[f'OIII5007/Hbeta'] = df_grid['OIII5007'] / df_grid['Hbeta']
        df_ratios[f'OII3726_29,OIII5007/Hbeta'] = (df_grid['OII3726_29'] + df_grid['OIII5007']) / df_grid['Hbeta']

    # ------plotting grid-----------
    for index, ycol in enumerate(ycol_arr):
        ax = axes[index][0]
        if 'R3' in ycol:
            args.ynum_line, args.yden_line = 'OIII', 'Hb'
        elif 'R23' in ycol:
            args.ynum_line, args.yden_line = 'OII,OIII', 'Hb'
        ax, ratio_name = plot_ratio_model(df_ratios, ax, args)
        ax = annotate_axes(ax, label_dict[args.quantity1], label_dict[ycol], ylim=[0, 1.1], args=args, hide_xaxis=index < len(ycol_arr) - 1)
    ax.legend(fontsize = args.fontsize / args.fontfactor, loc='upper left')

    # ------saving---------
    figname = f'ratio_2d_vs_1d_comparison_{args.binby_text}{args.deproject_text}{args.rescale_text}.png'
    save_fig(fig, args.fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')