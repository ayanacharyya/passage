'''
    Filename: compare_pjw_sed_fits.py
    Notes: Compares various runs of Peter's SED fit results with one another
    Author : Ayan
    Created: 01-09-26
    Example: run compare_pjw_sed_fits.py --system ssd
'''

from header import *
from util import *
from make_sfms_bins import read_passage_sed_catalog

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.5
    args.fontsize = 10

    # ---------reading in the SED catalogs----------------
    passage_catalog1 = 'SED_fits_v1.0.2_cosmosweb.fits'
    passage_catalog2 = 'SED_fits_v1.3.2-photonly_cosmosweb.fits'

    df1 = read_passage_sed_catalog(args.output_dir / 'catalogs' / passage_catalog1, include_cosmos2020=False)
    df2 = read_passage_sed_catalog(args.output_dir / 'catalogs' / passage_catalog2, include_cosmos2020=False)

    # -------obtaining comparison columns-----------------
    cols_to_match = ['field', 'id', 'cosmoswebid']
    df_merged = pd.merge(df1, df2, on=cols_to_match, how='inner')
    df_merged = df_merged.dropna(axis=0)
    print(f'\n{len(df_merged)} common galaxies found..')

    # ----------making comparison plots-----------------
    remaining_cols = list(set(df1.columns) - set(cols_to_match))
    cols_to_compare = np.sort([item for item in remaining_cols if not item.endswith('_u') and 'redshift' not in item])
    ncomp = len(cols_to_compare)

    ncols = 5
    nrows = int(np.ceil(ncomp / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(2. * ncols, 2. * nrows), layout='constrained')

    for index, thiscol in enumerate(cols_to_compare):
        print(f'Plotting comparison ({index + 1}/{ncomp}) {thiscol}..')
        row = index // ncols
        col = index % ncols
        ax = axes[row][col]

        ax.scatter(df_merged[f'{thiscol}_x'], df_merged[f'{thiscol}_y'], s=5, lw=0.5, c='cornflowerblue')
        if f'{thiscol}_u_x' in df_merged:
            ax.errorbar(df_merged[f'{thiscol}_x'], df_merged[f'{thiscol}_y'],xerr=df_merged[f'{thiscol}_u_x'], yerr=df_merged[f'{thiscol}_u_y'], fmt='none', lw=0.5, c='gray', zorder=-5)

        limits = [min(df_merged[f'{thiscol}_x'].min(), df_merged[f'{thiscol}_y'].min()), max(df_merged[f'{thiscol}_x'].max(), df_merged[f'{thiscol}_y'].max())]
        ax.plot(limits, limits, lw=0.5, ls='--', c='k')

        mad = (df_merged[f'{thiscol}_x'] - df_merged[f'{thiscol}_y']).abs().median()
        ax = annotate_axes(ax, f'{thiscol}_x', f'{thiscol}_y', xlim=limits, ylim=limits, label=f'MAD = {mad:.2f}', labelx=0.5, labely=0.95, args=args)

    fig.supxlabel(f'{passage_catalog1[:-5]}', fontsize=args.fontsize)
    fig.supylabel(f'{passage_catalog2[:-5]}', fontsize=args.fontsize)
    save_fig(fig, args.output_dir / 'plots', f'comparison_{passage_catalog1[:-5]}_{passage_catalog2[:-5]}.png', args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
