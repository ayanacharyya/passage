'''
    Filename: make_cosmos_plots.py
    Notes: Plots various quantities from the COSMOS2020 catalog
    Author : Ayan
    Created: 27-08-24
    Example: run make_cosmos_plots.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_cosmos_plots.py --xcol lp_mass_best --ycol lp_SFR_best --colorcol ez_z_phot
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------
def make_datashader_plot_mpl(df, args):
    '''
    Function to make data shader plot of y_field vs x_field, colored in bins of color_field
    This function is based on foggie.render.shade_maps.render_image()
    :return figure
    '''
    # --------to make the main datashader plot--------------------------
    fig, ax = plt.subplots()
    color_key = [mplcolors.to_hex(item) for item in mplcolormaps[args.colormap].colors]
    artist = dsshow(df, dsh.Point(args.xcol, args.ycol), dsh.mean(args.colorcol), norm='linear', cmap=color_key, x_range=(args.xmin, args.xmax), y_range=(args.ymin, args.ymax), vmin=args.cmin, vmax=args.cmax, aspect = 'auto', ax=ax)

    # ------to make the axes-------------
    ax.xaxis.set_label_text(label_dict[args.xcol])
    ax.yaxis.set_label_text(label_dict[args.ycol])

    # ------to make the colorbar axis-------------
    cax_xpos, cax_ypos, cax_width, cax_height = 0.2, 0.835, 0.25, 0.035
    cax = fig.add_axes([cax_xpos, cax_ypos, cax_width, cax_height])
    plt.colorbar(artist, cax=cax, orientation='horizontal')

    cax.set_xticklabels(['%.0F' % index for index in cax.get_xticks()])
    fig.text(cax_xpos + cax_width / 2, cax_ypos + cax_height + 0.005, label_dict[args.colorcol], ha='center', va='bottom')

    # ---------to annotate and save the figure----------------------
    figname = args.output_dir / 'plots' / f'cosmos2020_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'
    plt.savefig(figname, transparent=False)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
label_dict = {'lp_mass_best': r'log M$_*$/M$_{\odot}$', 'lp_SFR_best': r'log SFR (M$_{\odot}$/yr)', 'ez_z_phot': 'Redshift'}
bounds_dict = {'lp_mass_best': (6, 12), 'lp_SFR_best': (-3, 3), 'ez_z_phot': (0, 3)}
colormap_dict = defaultdict(lambda: 'viridis', ez_z_phot='plasma')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ----------to initialise axes limits and colormaps-----------------------
    args.xmin, args.xmax = bounds_dict[args.xcol]
    args.ymin, args.ymax = bounds_dict[args.ycol]
    args.cmin, args.cmax = bounds_dict[args.colorcol]
    args.colormap = colormap_dict[args.colorcol]

    # -----------------reading and trimming the dataframe--------------------------------------
    df = read_COSMOS2020_catalog(args=args)
    df = df[(df[args.xcol].between(args.xmin, args.xmax)) & (df[args.ycol].between(args.ymin, args.ymax)) & (df[args.colorcol].between(args.cmin, args.cmax))]
    df = df[(np.isfinite(df[args.xcol])) & (np.isfinite(df[args.ycol])) & (np.isfinite(df[args.colorcol]))]

    # -----------------making the plots--------------------------------------
    fig = make_datashader_plot_mpl(df, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
