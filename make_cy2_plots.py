'''
    Filename: make_cy2_plots.py
    Notes: Makes some test plots of a given Cy2 NIRISS WFSS field, based on PJW's catalog
    Author : Ayan
    Created: 02-04-26
    Example: run make_cy2_plots.py --drv pjw --field Par682
             run make_cy2_plots.py
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def read_catalog(filename, remove_cols=['cdf_z'], remove_ids=[7161, 7169, 7171]):
    '''
    Reads the grizli redshift catalog
    Returns pandas dataframe
    '''
    tab = Table.read(filename)
    tab.remove_columns(remove_cols)
    df = tab.to_pandas()

    df = df[~(df['id'].isin(remove_ids))]

    return df

# --------------------------------------------------------------------------------------------------------------------
def read_image(filename, remove_ids=[]):
    '''
    Reads the grizli segmentation map fits file
    Returns 2D array and header
    '''
    hdul = fits.open(filename)
    data = hdul[0].data
    header = hdul[0].header

    if 'seg' in str(filename): data = np.where(np.isin(data, remove_ids), np.nan, data)

    return data, header

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2

    # ----------determining directory paths---------------
    if args.drv == 'v0.5': args.drv = 'vpjw'
    if args.field == 'Par003': args.field = 'Par682'
    args.input_dir = args.root_dir / 'passage_data' / args.drv / args.field

    # ------------reading in files-------------------
    remove_ids = [7159, 7161, 7162, 7163, 7164, 7169, 7171, 7172, 7172, 7173, 7175, 7176, 7179, 7180, 7181, 7185, 7190] # these are found by eye, looking at the seg map and deciding which areas are affected by the stellar spikes
    df = read_catalog(args.input_dir / 'compiled_grizli_v1.0.0.fits', remove_ids=remove_ids)
    seg_map, seg_header = read_image(args.input_dir / 'passage-par682-ir_seg.fits', remove_ids=remove_ids)

    # ----------setup plot parameters------------------
    quant = 'redshift'
    delta_quant1 = 0.1
    delta_quant2 = 0.01
    target_quant_arr = [1.552, 1.656, 2.637, 4.876]

    # ----------setup figure------------------
    fig = plt.figure(figsize=(13, 6), layout='constrained')
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection=pywcs.WCS(seg_header))

    # ----------plot histogram------------------
    bins = np.arange(np.min(df[quant]), np.max(df[quant]), delta_quant1)
    ax1.hist(df[quant], bins=bins, color='grey', histtype='step', lw=1.)
    ax1 = annotate_axes(ax1, quant, 'Number of objects', args=args)

    # ----------plot seg map------------------
    ax2.imshow(seg_map, cmap='Greys', alpha=1, origin='lower', interpolation='nearest')
    ax2 = annotate_axes(ax2, 'RA', 'Dec', args=args)
    
    # -----------loop over overplots--------------------
    cmap_arr = ['Reds', 'Greens', 'Blues', 'Purples', 'Oranges']
    for index, target_quant in enumerate(target_quant_arr):
        color = mplcolormaps[cmap_arr[index]](0.5)
        df_sub = df[df[quant].between(target_quant - delta_quant2/2, target_quant + delta_quant2/2)]
        ax1.hist(df_sub[quant], bins=bins, color=color, label=f'{quant}={target_quant}' + r'$\pm$' + f'{delta_quant2/2}')

        highlight_map = np.where(np.isin(seg_map, df_sub['id'].values), seg_map, np.nan)
        ax2.imshow(highlight_map, cmap=cmap_arr[index], alpha=1, origin='lower', interpolation='nearest')
        ax2.scatter(df_sub['ra'], df_sub['dec'], marker='x', color=color, transform=ax2.get_transform('world'))

    ax1.legend(loc='upper right', fontsize=args.fontsize / args.fontfactor)

    # --------------save figure--------------
    save_fig(fig, args.input_dir / 'figs', f'{args.field}_{quant}_overdensity.png', args, dpi=800)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

