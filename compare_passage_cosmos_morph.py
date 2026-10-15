'''
    Filename: compare_passage_cosmos_morph.py
    Notes: Compare morphological parameters between PASSAGE and COSMSOWeb
    Author : Ayan
    Created: 15-08-26
    Example: run compare_passage_cosmos_morph.py --system ssd --do_all_fields --do_all_obj --cut_z_flag 4
'''

from header import *
from util import *
setup_plot_style()

from make_sfms_bins import get_binned_df, required_lines, sfms
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    if args.re_limit is None: args.re_limit = 2.
    colorby_col = 'redshift'
    
    # ---------reading in the master SED catalog----------------
    df, _, args = get_binned_df(args, required_lines=required_lines, sfms=sfms)

    # ------reading COSMOS-Web (subset) catalog-----------------
    df_cw = Table.read(args.root_dir / 'passage_data/v0.5/COSMOS/COSMOSWeb_mastercatalog_v1_photom_primary_only_objects_in_Huberty26.fits').to_pandas()

    # -----------merging--------------
    cols_to_extract = ['a_image', 'b_image', 'theta_image']
    df2 = df.merge(df_cw[np.hstack((['id'], cols_to_extract))], left_on='cosmoswebid', right_on='id', suffixes=['_passage', '_cweb'])
    df2['theta_image_passage'] = np.degrees(df2['theta_image_passage']) # originally PASSAGE PAs are in radians, but COSMOSWeb PAs are in degrees

    # --------setting up figure---------------
    ncols = len(cols_to_extract)
    fig, axes=plt.subplots(1, ncols, figsize=(12, 4))
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.15, top=0.99, wspace=0.3)

    for index, col in enumerate(cols_to_extract):
        ax = axes[index]
        mad = MAD(df2[f'{col}_passage'] - df2[f'{col}_cweb'])
        p = ax.scatter(df2[f'{col}_passage'], df2[f'{col}_cweb'], c='cornflowerblue' if colorby_col is None else df2[colorby_col], s=10, ec='k', lw=0.5, vmin=0.8, vmax=2.8, cmap='Greens')
        llim = min(ax.get_xlim()[0], ax.get_ylim()[0])
        ulim = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([llim, ulim], [llim, ulim], ls='dashed', c='k', lw=1)
        ax = annotate_axes(ax, 'PASSAGE', 'COSMOS-Web', label=f'{col}; MAD={mad:.1f}', args=args, p=p, 
                           hide_cbar=colorby_col is None or index < ncols - 1, clabel=colorby_col)

    # --------------saving figure----------------
    args.required_lines_text = '_lines_' + ','.join(required_lines) if len(required_lines) > 0 else ''
    args.zcut_text = f'_zcut_{args.zcut}' if args.zcut is not None else ''
    args.cut_z_flag_text = f'_zflag_cut_{args.cut_z_flag}{args.zcut_text}' if args.cut_z_flag is not None else f'{args.zcut_text}'

    if args.do_all_fields: root_dir = args.output_dir
    else: root_dir = args.output_dir / args.field
    args.stacking_dir = root_dir / f'stacking{args.required_lines_text}{args.cut_z_flag_text}'

    save_fig(fig, args.stacking_dir, 'passage_cosmosweb_morph_comparison.png', args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
