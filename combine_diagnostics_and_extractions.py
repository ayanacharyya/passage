'''
    Filename: combine_diagnostics_and_extractions.py
    Notes: Combine extracted 1D and 2D spectra and diagnostic maps of same object args.ids, for a given object/s in a given field
    Author : Ayan
    Created: 26-07-24
    Example: run combine_diagnostics_and_extractions.py --field Par61 --plot_radial_profiles --only_seg --snr_cut 3 --keep --id 38
             run combine_diagnostics_and_extractions.py --field Par61 --plot_radial_profiles --only_seg --snr_cut 3 --keep --do_all_obj
             run combine_diagnostics_and_extractions.py --field Par28 --drv 0.5 --plot_radial_profiles --only_seg --do_all_obj
    Afterwards, to make the animation: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par061/diagnostics_and_extractions/ --rootname Par061_*_diagnostics_and_extractions.png --delay 0.1
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
    description_text2 = f'diagnostics_and_extractions'

    extraction_path = args.input_dir / args.drv / args.field / 'Extractions'
    products_path = args.input_dir / args.drv / args.field / 'Products'
    input_dir = args.output_dir / args.field / f'{description_text}'
    quant_arr = ['line', 'stack', 'full']
    id_arr = args.id

    if args.do_all_obj:
        try:
            catalog_file = extraction_path / f'{args.field}-ir.cat.fits'
            catalog = GTable.read(catalog_file)
        except:
            catalog_file = args.input_dir / args.drv / args.field / 'Products' / f'{args.field}_photcat.fits'
            catalog = GTable.read(catalog_file)

        args.id_arr = catalog['NUMBER'] if 'NUMBER' in catalog.columns else catalog['id']
    else:
        args.id_arr = args.id

    if args.start_id: args.id_arr = args.id_arr[args.start_id - 1:]
    if len(args.id_arr) > 10: args.hide = True # if too many plots, do not display them, just save them

    output_dir = args.output_dir / args.field / f'{description_text2}'
    output_dir.mkdir(parents=True, exist_ok=True)


    # ------------read in grizli images and plot comparison--------------------------------
    for index, args.id in enumerate(args.id_arr):
        print(f'\nCombining extractions and diagnostic maps for id {args.id} which is {index+1} out of {len(args.id_arr)}..')
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
        if not os.path.exists(alternate_path / diag_figname) and not os.path.exists(input_dir / diag_figname):
            print('Could not diagnostic maps image, so..')

            command = ['python', 'make_diagnostic_maps.py', '--field', f'{args.field}', '--id', f'{args.id}']
            if args.plot_radial_profiles: command += ['--plot_radial_profiles']
            if args.only_seg: command += ['--only_seg']
            if args.snr_cut is not None: command += ['--snr_cut', f'{args.snr_cut}']
            if args.vorbin:
                command += ['--vorbin']
                command += ['--voronoi_line', f'{args.voronoi_line}']
                command += ['--voronoi_snr', f'{args.voronoi_snr}']

            print(f'Trying to run following command:\n{" ".join(command)}..')
            dummy = subprocess.run(command)

        try:
            diag = mpimg.imread(alternate_path / diag_figname)
            axes_diag.imshow(diag, origin='upper')
        except FileNotFoundError:
            try:
                diag = mpimg.imread(input_dir / diag_figname)
                axes_diag.imshow(diag, origin='upper')
            except FileNotFoundError:
                print('Could not find diagnostic plots, so skipping')
                continue

        for ax in fig.axes:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
        fig.subplots_adjust(left=0.01, bottom=0.01, top=0.99, right=0.99, hspace=0.0, wspace=0.0)

        figname =  output_dir / f'{args.field}_{args.id:05d}_{description_text2}{vorbin_text}.png'
        fig.savefig(figname, dpi=400)
        print(f'Saved figure at {figname}')
        if args.hide: plt.close('all')
        else: plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
