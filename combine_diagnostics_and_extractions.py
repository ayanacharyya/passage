'''
    Filename: combine_diagnostics_and_extractions.py
    Notes: Combine extracted 1D and 2D spectra and diagnostic maps of same object args.ids, for a given object/s in a given field
    Author : Ayan
    Created: 26-07-24
    Example: run combine_diagnostics_and_extractions.py --field Par61 --keep --id 38
             run combine_diagnostics_and_extractions.py --field Par61 --keep --do_all_obj
    Afterwards, to make the animation: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par061/ --rootname Par061_*_diagnostics_and_extractions.png --delay 0.1
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------determining directories and global variables---------
    extraction_path = args.input_dir / args.field / 'Extractions'
    output_dir = args.output_dir / args.field
    quant_arr = ['line', 'stack', 'full']
    id_arr = args.id

    args.plot_radial_profiles, args.only_seg, args.snr_cut = True, True, 3 #
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    snr_text = f'_snr{args.snr_cut:.1f}' if args.snr_cut is not None else ''

    if args.do_all_obj:
        try:
            catalog_file = extraction_path / f'{args.field}-ir.cat.fits'
            catalog = GTable.read(catalog_file)
        except:
            catalog_file = args.input_dir / args.field / 'Products' / f'{args.field}_photcat.fits'
            catalog = GTable.read(catalog_file)

        args.id_arr = catalog['NUMBER'] if 'NUMBER' in catalog.columns else catalog['id']
    else:
        args.id_arr = args.id

    if args.start_id: args.id_arr = args.id_arr[args.start_id - 1:]
    if len(args.id_arr) > 10: args.hide = True # if too many plots, do not display them, just save them

    # ------------read in grizli images and plot comparison--------------------------------
    for index, args.id in enumerate(args.id_arr):
        print(f'\nCombining extractions and diagnostic maps for id {args.id} which is {index+1} out of {len(args.id_arr)}..')

        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(len(quant_arr), 2)
        axes_ex = [fig.add_subplot(gs[item, 0]) for item in range(len(quant_arr))]
        axes_diag = fig.add_subplot(gs[:, 1])

        for ind, quant in enumerate(quant_arr):
            print(f'Reading in {quant} png files..')
            try:
                ex = mpimg.imread(extraction_path / f'{args.field}_{args.id:05d}.{quant}.png')
                axes_ex[ind].imshow(ex, origin='upper')
            except FileNotFoundError:
                print('Could not find file, so skipping')
                continue

        print(f'Reading in diagnostic png files..')

        try:
            diag = mpimg.imread(output_dir / f'{args.field}_{args.id:05d}_all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}.png')
            axes_diag.imshow(diag, origin='upper')
        except FileNotFoundError:
            print('Could not find file, so skipping')
            continue

        for ax in fig.axes:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
        fig.subplots_adjust(left=0.01, bottom=0.01, top=0.99, right=0.99, hspace=0.0, wspace=0.0)

        figname = output_dir / f'{args.field}_{args.id:05d}_diagnostics_and_extractions.png'
        fig.savefig(figname, dpi=400)
        print(f'Saved figure at {figname}')
        if args.hide: plt.close('all')
        else: plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
