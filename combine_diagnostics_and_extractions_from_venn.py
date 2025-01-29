'''
    Filename: combine_diagnostics_and_extractions_from_venn.py
    Notes: Combine extracted 1D and 2D spectra and diagnostic maps of same object based on a given list of objects (IDs and fields) from the intersection of Venn diagram (produced by get_field_stats.py
    Author : Ayan
    Created: 26-07-24
    Example: run combine_diagnostics_and_extractions_from_venn.py --plot_radial_profiles --only_seg --do_all_fields --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --SNR_thresh 10
             run combine_diagnostics_and_extractions_from_venn.py --plot_radial_profiles --only_seg --do_all_fields --line_list OIII,Ha --plot_conditions EW,mass,PA --EW_thresh 300
             run combine_diagnostics_and_extractions_from_venn.py --plot_radial_profiles --only_seg --do_all_fields --line_list OIII,Ha --plot_conditions SNR,mass,PA,a_image --SNR_thresh 10 --a_thresh 2.4
             run combine_diagnostics_and_extractions_from_venn.py --plot_radial_profiles --only_seg --drv 0.5 --field Par028 --vorbin --voronoi_line Ha --voronoi_snr 5 --do_not_correct_pixel --use_O3S2 --plot_conditions SNR --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --clobber
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------determining directories and global variables---------
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    snr_text = f'_snr{args.snr_cut:.1f}' if args.snr_cut is not None else ''
    pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
    vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'
    description_text = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'
    fields_text = 'allpar' if args.do_all_fields else ','.join(args.field_arr)

    plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
    plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
    description_text2 = f'diagnostics_and_extractions_{plot_conditions_text}'

    df_int_filename = args.output_dir / 'catalogs' / f'{fields_text}_{args.drv}_venn_{plot_conditions_text}_df.txt'
    df_int = pd.read_csv(df_int_filename)
    df_int = df_int.sort_values(by=['field', 'objid']).reset_index(drop=True)
    if args.use_only_good and args.drv == 'v0.5' and set(args.plot_conditions) == set(['SNR']) and set(args.line_list) == set(['OIII', 'Ha', 'OII', 'Hb', 'SII']):
        df_int = df_int[df_int['objid'].isin([1303, 1934, 2734, 2867, 300, 2903])].reset_index(drop=True)  # only choosing the pre-determined good galaxies
        print(f'Using only the pre-determined good galaxies, and there are {len(df_int)} of them..')

    output_dir = args.output_dir / f'{description_text2}'
    output_dir.mkdir(parents=True, exist_ok=True)

    quant_arr = ['line', 'stack', 'full']
    if len(df_int) > 20: args.hide = True

    # -----------looping over all objects----------------------------------
    for index, row in df_int.iterrows():
        print(f'\nCombining extractions and diagnostic maps for object {index + 1} of {len(df_int)}...')
        args.field = row['field']
        args.id = row['objid']

        extraction_path = args.input_dir / args.drv / args.field / 'Extractions'
        products_path = args.input_dir / args.drv / args.field / 'Products'
        diag_img_dir = args.output_dir / args.field / f'{description_text}'

        # ------------read in grizli images and plot comparison--------------------------------
        alternate_path = args.output_dir / args.field / f'{args.id:05d}'
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(len(quant_arr), 2)
        axes_ex = [fig.add_subplot(gs[item, 0]) for item in range(len(quant_arr))]
        axes_diag = fig.add_subplot(gs[:, 1])

        found_individual_plots = False
        for ind, quant in enumerate(quant_arr):
            print(f'Reading in {quant} png files..')
            try:
                ex = mpimg.imread(alternate_path / f'{args.field}_{args.id:05d}.{quant}.png')
                axes_ex[ind].imshow(ex, origin='upper')
                found_individual_plots = True
            except FileNotFoundError:
                try:
                    ex = mpimg.imread(extraction_path / f'{args.field}_{args.id:05d}.{quant}.png')
                    axes_ex[ind].imshow(ex, origin='upper')
                    found_individual_plots = True
                except FileNotFoundError:
                    print('Could not find file, so trying to plot pre-made plot-group')
                    continue

        # ------------------opening grizli-made group of plots if individual plots not found-----------------------
        if not found_individual_plots:
            try:
                print(f'Reading in group-plot png file..')
                ex = mpimg.imread(products_path / 'plots' / f'{args.field}_{args.id:05d}.png')
                ex = ex[int(np.shape(ex)[0] * 0.23):, :, :] # to crop out the additional direct image that this group of plot includes, and we do not need that

                # ------to delete existing fig with 3 separate axis and make a new one with only axis because this is already a group of plots------
                plt.close(fig)
                fig = plt.figure(figsize=(12, 6))
                gs = fig.add_gridspec(1, 2)
                axes_merged, axes_diag = [fig.add_subplot(gs[:, item]) for item in range(2)]
                axes_merged.imshow(ex, origin='upper') # use origin ='lower' for cases when PASSAGEPipe originally stitched these plots in a flipped way

            except FileNotFoundError:
                print('Could not find grizli-made plot-group, so skipping')
                continue

        # ------------------opening diagnostic plots-----------------------
        print(f'Reading in diagnostic png files..')

        diag_figname = f'{args.field}_{args.id:05d}_{description_text}{vorbin_text}.png'
        if args.clobber or (not os.path.exists(diag_img_dir / diag_figname) and not os.path.exists(alternate_path / diag_figname)):
            print('Could not diagnostic maps image, so..')

            command = ['python', 'make_diagnostic_maps.py', '--field', f'{args.field}', '--id', f'{args.id}', '--drv', f'{args.drv}', '--plot_AGN_frac']
            if args.plot_radial_profiles: command += ['--plot_radial_profiles']
            if args.only_seg: command += ['--only_seg']
            if args.snr_cut is not None: command += ['--snr_cut', f'{args.snr_cut}']
            if args.do_not_correct_pixel: command += ['--do_not_correct_pixel']
            if args.do_not_correct_flux: command += ['--do_not_correct_flux']
            if args.use_O3S2: command += ['--use_O3S2']
            if args.vorbin:
                command += ['--vorbin']
                command += ['--voronoi_line', f'{args.voronoi_line}']
                command += ['--voronoi_snr', f'{args.voronoi_snr}']

            print(f'Trying to run following command:\n{" ".join(command)}')
            dummy = subprocess.run(command)

        try:
            diag = mpimg.imread(diag_img_dir / diag_figname)
            axes_diag.imshow(diag, origin='upper')
        except FileNotFoundError:
            try:
                diag = mpimg.imread(alternate_path / diag_figname)
                axes_diag.imshow(diag, origin='upper')
            except FileNotFoundError:
                print('Could not find diagnostic plots, so skipping')
                continue

        for ax in fig.axes:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
        fig.subplots_adjust(left=0.01, bottom=0.01, top=0.99, right=0.99, hspace=0.0, wspace=0.0)

        figname =  output_dir / f'{args.field}_{args.id:05d}_diagnostics_and_extractions_{description_text}{vorbin_text}.png'
        fig.savefig(figname, dpi=400)
        print(f'Saved figure at {figname}')
        if args.hide: plt.close('all')
        else: plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
