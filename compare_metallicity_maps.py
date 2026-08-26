'''
    Filename: compare_metallicity_maps.py
    Notes: Plots the stacked metallicity maps from different ways of stacking, and plots the OIII scaling factor used in each case
    Author : Ayan
    Created: 25-08-26
    Example: run compare_metallicity_maps.py --system ssd
'''
from header import *
from util import *
setup_plot_style()
from plot_stacked_maps import plot_2D_map, read_metallicity_map

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def parse_bin_intervals(file_path):
    # Convert Path object to string
    filename = str(file_path)
    
    # Matches float numbers (including negative ones like -0.500)
    pattern = r'delta_sfms_bin_([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)_mass_bin_([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)'
    match = re.search(pattern, filename)
    
    if match:
        sfms_low, sfms_high = float(match.group(1)), float(match.group(2))
        mass_low, mass_high = float(match.group(3)), float(match.group(4))
        return (sfms_low, sfms_high, mass_low, mass_high)
    
    return (float('inf'), float('inf'), float('inf'), float('inf'))

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.5
    args.extent = (0, args.re_limit, 0, args.re_limit) # because these are folded maps
        
    # -----------define colorbar properties-----------
    Zmin, Zmax, cmap = 7, 9, 'plasma'
    limits = [-18, -15]
    
    # ------------determining directory and files-------------
    df_scaling1 = pd.read_csv(args.output_dir / f'stacking_lines_OIII,Ha_zflag_cut_4/scaling_data.csv')
    folder1 = args.output_dir / 'stacking_with_halpha_correction_by_constant' / f'stacking_lines_OIII,Ha_zflag_cut_4' / 'adap_binby_distance_8_sfms_PASSAGE_mass_4' / 'maps_nodeproject'
    files1 = glob.glob(str(folder1) + f'/stacked_folded_metallicity_map_NB_nodeproject_delta_sfms_bin_*_mass_bin_*.fits')
    all_sorted_files1 = sorted(files1, key=parse_bin_intervals)

    df_scaling2 = pd.read_csv(args.output_dir / f'stacking_lines_OIII,Ha_zflag_cut_4/scaling_data_oldlinecut.csv')
    folder2 = args.output_dir / 'stacking_with_halpha_correction_by_constant' / f'stacking_lines_OIII,Ha_zflag_cut_4_oldlinecut' / 'adap_binby_distance_8_sfms_PASSAGE_mass_4' / 'maps_nodeproject'
    #df_scaling2 = pd.read_csv(args.output_dir / f'stacking_lines_OIII,Ha_zflag_cut_4/scaling_data.csv')
    #folder2 = args.output_dir / f'stacking_lines_OIII,Ha_zflag_cut_4' / 'adap_binby_distance_8_sfms_PASSAGE_mass_4' / 'maps_nodeproject'

    files2 = glob.glob(str(folder2) + f'/stacked_folded_metallicity_map_NB_nodeproject_delta_sfms_bin_*_mass_bin_*.fits')
    all_sorted_files2 = sorted(files2, key=parse_bin_intervals)

    df_scaling_arr = [df_scaling1, df_scaling2]

    nfiles = max(len(all_sorted_files1), len(all_sorted_files2))
    ncols = 8
    npages = int(np.ceil(nfiles / ncols))

    # ------------looping over several images------------------
    for img_index in range(npages):
        print(f'\n\nDoing page {img_index + 1} of {npages}..')
        sorted_files1 = all_sorted_files1[img_index * ncols : (img_index + 1) * ncols]
        sorted_files2 = all_sorted_files2[img_index * ncols : (img_index + 1) * ncols]

        # ---------setting up figure--------------
        nrows = 4
        ncols = ncols
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 1.5, nrows * 1.5), sharex='row', sharey='row')
        fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.1)

        # ----------looping over bins-----------
        for index, (logOH_file1, logOH_file2) in enumerate(zip(sorted_files1, sorted_files2)):
            print(f'\n\nCommencing pair {index} of {ncols}..')
            for index2, logOH_file in enumerate([logOH_file1, logOH_file2]):
                sfms_low, sfms_high, mass_low, mass_high = parse_bin_intervals(logOH_file)
                print(f'\nCommencing bins: delta_sfms=[{sfms_low}-{sfms_high}], log_mass=[{mass_low}-{mass_high}]')

                logOH_map, _, _ = read_metallicity_map(logOH_file, args)

                line_map_filename = logOH_file.replace('folded_metallicity_map_NB', 'maps')
                line_hdu = fits.open(line_map_filename)

                df_id = Table(line_hdu['OIII_IDS'].data).to_pandas()
                df_id[['field', 'id']] = df_id['OBJID'].str.split('-', expand=True)
                df_id['id'] = df_id['id'].astype(int)
                df_id.drop(columns=['OBJID'], inplace=True)
                df_id = df_id.merge(df_scaling_arr[index2], on=['field', 'id'], how='left')

                # -------making metallicity plots--------------
                ax = axes[index2][index]
                ax = plot_2D_map(logOH_map, ax, '', args, cmap=cmap, clabel=r'$\log$(O/H) + 12', takelog=False, vmin=Zmin, vmax=Zmax, hide_cbar=True, hide_yaxis=index > 0, hide_xaxis=index2 == 0)
                ax.set_box_aspect(1)
                label = f'{sfms_low:.2f}<' + r'$\delta_{SFMS}$' + f'<{sfms_high:.2f}\n{mass_low:.2f}<' + r'$\log M_*$' + f'<{mass_high:.2f}'
                ax.text(0.9, 0.95, label, ha='right', va='top', fontsize = args.fontsize / (args.fontfactor**2), transform=ax.transAxes)

                # ------making the scaling factor plots------------
                ax = axes[index2 + 2][index]
                ax.scatter(np.log10(df_id['scaling_flux_sum']), np.log10(df_id['scaling_flux_speccat']), s=5, lw=0.5, c='cornflowerblue')
                ax.plot(limits, limits, ls='dashed', c='k', lw=0.5)
                ax.set_box_aspect(1)
                ax = annotate_axes(ax, '2D Sum', 'Grizli', xlim=limits, ylim=limits, args=args, label=f'{len(df_id)}', labelx=0.05, labely=0.9, hide_xaxis=index2 == 0, hide_yaxis=index > 0)

        save_fig(fig, args.output_dir, f'logOH_map_comparison_{img_index}.png', args, dpi=300)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
