'''
    Filename: compare_extractions.py
    Notes: Compare extracted 1D and 2D spectra of same object args.ids, for a given object/s in a given field
    Author : Ayan
    Created: 26-07-24
    Example: run compare_extractions.py --field Par51 --keep --id 420
             run compare_extractions.py --field Par51 --keep --re_extract --do_all_obj
'''

from header import *
from util import *
import matplotlib.image as mpimg

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------
def plot_redshift_comparison(field, directory='/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/'):
    '''
    Plot redshift distribution of newer and older extractions
    Saves plot as png file, and returns figure handle
    '''
    path = directory / field / 're_extracted' / 'comparisons'
    df = pd.read_table(path / f'{field}_all_diag_results_old_vs_new_comp.txt')

    fig, ax = plt.subplots(111)

    ax.hist(df['redshift_old'], bins=20, range=(0, 8), histtype='step', fill=False, lw=2, label='0.5<z<3.8')
    ax.hist(df['redshift_new'], bins=20, range=(0, 8), histtype='step', fill=False, lw=2, label='0.1<z<8')

    plt.legend()
    ax.set_xlabel('Redshift')
    ax.set_ylabel('No. of objects')

    ax.text(8, ax.get_ylim()[1] * 0.7, f'Par051; total {len(df)} objects', c='k', ha='right', va='top')

    figname = path / f'{field}_redshift_old_vs_new_comp.png'
    fig.savefig(figname)
    print(f'Saved figure at {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------determining directories and global variables---------
    orig_path = args.input_dir / args.field / 'Extractions'
    output_dir = args.output_dir / args.field
    re_extraction_path = output_dir / 're_extracted'
    quant_arr = ['line', 'stack', 'full']
    id_arr = args.id

    args.plot_radial_profiles, args.only_seg, args.snr_cut = True, True, 3 #
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    snr_text = f'_snr{args.snr_cut:.1f}' if args.snr_cut is not None else ''

    # fig_zcomp = plot_redshift_comparison(args.field, directory=args.output_dir)

    if args.do_all_obj:
        if args.re_extract:
            args.id_arr = ids_to_re_extract_dict[args.field]
        else:
            try:
                catalog_file = orig_path / f'{args.field}-ir.cat.fits'
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
        print(f'\nComparing extractions for id {args.id} which is {index+1} out of {len(args.id_arr)}..')
        new_path = re_extraction_path / f'{args.id:05d}'

        fig = plt.figure(figsize=(8, 8))
        gs = fig.add_gridspec(2 * len(quant_arr), 2)
        axes_ex = [fig.add_subplot(gs[item, 0]) for item in range(2 * len(quant_arr))]
        axes_diag = [fig.add_subplot(gs[item * len(quant_arr) : len(quant_arr) * (item + 1), 1]) for item in range(2)]

        for ind, quant in enumerate(quant_arr):
            print(f'Reading in {quant} png files..')

            old_ex = mpimg.imread(orig_path / f'Par051_{args.id:05d}.{quant}.png')
            axes_ex[ind].imshow(old_ex, origin='upper')

            new_ex = mpimg.imread(new_path / f'Par051_{args.id:05d}.{quant}.png')
            axes_ex[ind + len(quant_arr)].imshow(new_ex, origin='upper')

        print(f'Reading in diagnostic png files..')

        old_diag = mpimg.imread(output_dir / f'{args.field}_{args.id:05d}_all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}.png')
        axes_diag[0].imshow(old_diag, origin='upper')

        new_diag = mpimg.imread(re_extraction_path / f'{args.field}_{args.id:05d}_all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}.png')
        axes_diag[1].imshow(new_diag, origin='upper')

        fig.subplots_adjust(left=0.01, bottom=0.01, top=0.99, right=0.99, hspace=0.05, wspace=0.0)

        figname = re_extraction_path / 'comparisons' / f'{args.field}_{args.id:05d}_extraction_old_vs_new_comp.png'
        fig.savefig(figname, dpi=800)
        print(f'Saved figure at {figname}')
        if args.hide: plt.close('all')
        else: plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
