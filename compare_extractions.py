'''
    Filename: compare_extractions.py
    Notes: Compare extracted 1D and 2D spectra of same object args.ids, for a given object/s in a given field
    Author : Ayan
    Created: 26-07-24
    Example: run compare_extractions.py --field Par51 --keep --id 420
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
    re_extraction_path = args.output_dir / args.field / 're_extracted'
    quant_arr = ['line', 'stack', 'full']
    id_arr = args.id

    # fig_zcomp = plot_redshift_comparison(args.field, directory=args.output_dir)

    # ------------read in grizli images and plot comparison--------------------------------
    for index, args.id in enumerate(id_arr):
        print(f'\nFitting spectra for id {args.id} which is {index+1} out of {len(id_arr)}..')

        fig, axes = plt.subplots(len(quant_arr), 2, figsize=(10, 8))
        for ind, quant in enumerate(quant_arr):
            print(f'Reading in {quant} png file')
            old = mpimg.imread(orig_path / f'Par051_{args.id:05d}.{quant}.png')
            axes[ind, 0].imshow(old, origin='upper')
            new = mpimg.imread(re_extraction_path / f'{args.id:05d}' / f'Par051_{args.id:05d}.{quant}.png')
            axes[ind, 1].imshow(new, origin='upper')

        fig.subplots_adjust(left=0.05, bottom=0.01, top=0.99, right=0.98, hspace=0.01, wspace=0.1)

        figname = re_extraction_path / 'comparisons' / f'{args.field}_{args.id:05d}_extraction_old_vs_new_comp.png'
        fig.savefig(figname)
        print(f'Saved figure at {figname}')
        plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
