'''
    Filename: read_line_catalog.py
    Notes: Makes 1D and 2D spectra from beam files, for a given object in a given field
    Author : Ayan
    Created: 18-06-24
    Last modified: 18-06-24
    Example: run read_line_catalog.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par21
             run read_line_catalog.py --zmin 0 --zmax 0.7 --line_list OII,OIII,Ha --keep --field Par21 --mag_lim 25
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def plot_zhist(df, args):
    '''
    Plots the redshift distribution function for a given input redshift catalog
    '''
    if not args.keep: plt.close('all')
    z_min, z_max = 0, 6

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(right=0.95, top=0.95, bottom=0.1, left=0.1)
    hist, bin_edges = np.histogram(df['redshift'], bins=args.nbins, range=(z_min, z_max))
    ax.hist(df['redshift'], bins=bin_edges, color='k', histtype='step')
    ax.set_xlim(z_min, z_max) #

    line_color_dict = defaultdict(lambda:'grey', OII='mediumpurple', Hb='lightskyblue', OIII='palegreen', Ha='salmon', SII='maroon', SIII='sienna', PaD='peru', PaG='tan', PaB='orange', PaA='firebrick')
    df_all_lines = df

    if args.line_list != 'all':
        for index, thisline in enumerate(args.line_list):
            z_limits, filters = get_zranges_for_filters(thisline)
            print(f'{thisline} is captured between ' + ', '.join(['z=[%.2f, %.2f]' % (z[0], z[1]) for z in z_limits]))

            df_line = pd.DataFrame()
            for i,z in enumerate(z_limits):
                df_line = pd.concat([df_line, df[df['redshift'].between(z[0], z[1])]])
                ax.axvspan(z[0], z[1], ymin=0.05, ymax=0.9 - 0.1 * index, alpha=0.3, color=line_color_dict[thisline])
                if index == 0: ax.text(z[1], ax.get_ylim()[0] + (0.9 - 0.1 * index) * np.diff(ax.get_ylim())[0], filters[i], c=line_color_dict[thisline], ha='right', va='bottom', fontsize=args.fontsize / 1.5)

            df_all_lines = df_all_lines[df_all_lines['objid'].isin(df_line['objid'])]
            ax.text(ax.get_xlim()[1] * 0.99, ax.get_ylim()[1] * 0.3 - 0.05 * np.diff(ax.get_ylim())[0] * index, '%s: %d' %(thisline, len(df_line)), c=line_color_dict[thisline], ha='right', va='bottom', fontsize=args.fontsize)

    ax.hist(df_all_lines['redshift'], bins=bin_edges, color='grey')

    ax.text(ax.get_xlim()[0] * 1.1 + 0.1, ax.get_ylim()[1] * 0.98, f'Total: {len(df)}', c='k', ha='left', va='top', fontsize=args.fontsize)
    ax.text(ax.get_xlim()[0] * 1.1 + 0.1, ax.get_ylim()[1] * 0.93, f'All {len(args.line_list)} lines: {len(df_all_lines)}', c='grey', ha='left', va='top', fontsize=args.fontsize)

    if args.mag_lim is not None:
        mag_col = 'f140w_mag'
        selection_color = 'saddlebrown'
        df_all_lines = df_all_lines[df_all_lines[mag_col] <= args.mag_lim]

        ax2 = ax.inset_axes([0.68, 0.45, 0.3, 0.3]) # [lower left corner coordinates, width, height] # make inset plot for magnitude distribution
        hist2, bin_edges2 = np.histogram(df[mag_col], bins=args.nbins)
        ax2.hist(df[mag_col], bins=bin_edges2, color='k', histtype='step')
        ax2.hist(df[df[mag_col] <= args.mag_lim][mag_col], bins=bin_edges2, color='grey')
        ax2.hist(df_all_lines[mag_col], bins=bin_edges2, color=selection_color)
        ax2.set_xlabel(mag_col, fontsize=args.fontsize/2)
        ax2.set_ylabel('# of objects', fontsize=args.fontsize/2)

        ax.hist(df_all_lines['redshift'], bins=bin_edges, color=selection_color)

        ax.text(ax.get_xlim()[1] * 0.99, ax.get_ylim()[1] * 0.88, f'After mag cut at {args.mag_lim}: {len(df[df[mag_col] <= args.mag_lim])}', c='grey', ha='right', va='top', fontsize=args.fontsize)
        ax.text(ax.get_xlim()[1] * 0.99, ax.get_ylim()[1] * 0.83, f'All {len(args.line_list)} lines: {len(df_all_lines)}', c=selection_color, ha='right', va='top', fontsize=args.fontsize)

    ax.set_xlabel('Redshift', fontsize=args.fontsize)
    ax.set_ylabel('# of objects', fontsize=args.fontsize)
    ax.set_title(args.field, fontsize=args.fontsize)

    figname = args.output_dir / 'zdist_histogram.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')

    plt.show(block=False)

    if args.line_list != 'all': return fig, df_all_lines
    else: return fig

# -------------------------------------------------------------------------------------------------------
def get_catalog(filename, nheader=217):
    '''
    Reads and returns the redshift catalog for a given input filename
    Returns pandas dataframe
    '''
    df_col = pd.read_table(filename, index_col=False, header=None, delim_whitespace=True, nrows=nheader)
    df = pd.read_table(filename, index_col=False, header=None, delim_whitespace=True, skiprows=nheader)
    df.columns = df_col[1].tolist()
    df = df[['objid', 'redshift', 'redshift_error', 'ra_obj', 'dec_obj', 'f140w_mag', 'snr_tot_others', 'chisq']]

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    args.input_dir = args.input_dir / args.field
    args.output_dir = args.output_dir / args.field
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # -----reading the file---------------
    filenames = list(args.input_dir.glob('*catalog*'))
    if len(filenames) > 0: df = get_catalog(filenames[0], nheader=217)
    else: sys.exit('No catalog file in %s' %args.input_dir)

    # -----computing stats---------------
    df_subz = df[df['redshift'].between(args.zmin, args.zmax)]
    #print(df)
    print(f'Field {args.field} has {len(df)} objects..out of which {len(df_subz)} objects, i.e., {len(df_subz) * 100 / len(df)}%, are within z=[{args.zmin:.2f}, {args.zmax:.2f}]')
    if args.mag_lim is not None: print(f'..out of which, {len(df[df["f140w_mag"] <= args.mag_lim])} are below magnitude limit of {args.mag_lim}')

    # -----making histogram---------------
    if args.line_list != 'all': fig, df_all_lines = plot_zhist(df, args)
    else: fig = plot_zhist(df, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
