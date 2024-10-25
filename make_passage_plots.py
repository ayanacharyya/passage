'''
    Filename: make_passage_plots.py
    Notes: Plots various quantities from a given passage dataframe (produced by get_field_stats.py by intersecting Venn diagrams of various properties)
    Author : Ayan
    Created: 10-09-24
    Example: run make_passage_plots.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_passage_plots.py --plot_conditions EW,mass,PA --xcol lp_mass --ycol lp_SFR --colorcol OIII_EW
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def convert_redshift_to_time(redshift):
    '''
    Computes time (age of the Universe) for a given redshift, and returns the time (in Gyr)
    Currently not complete
    '''
    if redshift == 1.8: time = 10 # Gyr
    elif redshift == 1.2: time = 11 # Gyr

    return time

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso22(ax, redshift, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+22 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a0, a1, b0, b1, b2 = 0.2, -0.034, -26.134, 4.722, -0.1925  # Table 2, Eq 10

    log_mass = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)
    time = convert_redshift_to_time(redshift) # Gyr
    log_SFR = (a1 * time + b1) * log_mass + b2 * (log_mass) ** 2 + b0 + a0 * time

    ax.plot(log_mass, log_SFR, c=color, lw=2, label=f'Popesso+22: z~{redshift}')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_MZR_literature(ax):
    '''
    Computes empirical MZRs based on several studies
    Then overplots them on a given existing axis handle
    Returns axis handle
    '''
    log_mass = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)

    # ---------Andrews+13-----------
    log_mass_TO, logOH_asym, gamma = 8.901, 8.798, 0.640 # First row of Table 3 in https://iopscience.iop.org/article/10.1088/0004-637X/765/2/140
    log_OH = logOH_asym - np.log10(1 + (log_mass_TO / log_mass)**gamma) # Eq 5 of A13
    ax.plot(log_mass, log_OH, c='cornflowerblue', lw=2, label=f'Andrews+13')

    # ---------Zahid+14-----------
    log_mass_TO, logOH_asym, gamma = 9.08, 10.06, 0.61 # COSMOS row of Table 2 in https://iopscience.iop.org/article/10.1088/0004-637X/791/2/130#apj498704t2
    log_OH = logOH_asym + np.log10(1 - np.exp(-(log_mass/log_mass_TO)**gamma)) # Eq 5 of Z14
    ax.plot(log_mass, log_OH, c='blue', lw=2, label=f'Zahid+14:COSMOS')

    return ax

# --------------------------------------------------------------------------------------------------------------------
label_dict = {'lp_mass': r'log M$_*$/M$_{\odot}$', 'lp_SFR': r'log SFR (M$_{\odot}$/yr)', 'ez_z_phot': 'Redshift', 'redshift': 'Redshift'}
bounds_dict = {'lp_mass': (6, 12), 'lp_SFR': (-3, 3), 'ez_z_phot': (0, 3), 'redshift': (0, 3)}
colormap_dict = defaultdict(lambda: 'viridis', ez_z_phot='plasma')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ----------to initialise axes limits and colormaps-----------------------
    # args.xmin, args.xmax = bounds_dict[args.xcol]
    # args.ymin, args.ymax = bounds_dict[args.ycol]
    # args.cmin, args.cmax = bounds_dict[args.colorcol]
    # args.colormap = colormap_dict[args.colorcol]

    # -------reading in dataframe produced by get_field_stats.py----------------
    df_infilename = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df.txt'
    df = pd.read_csv(df_infilename)

    # -----------plotting stuff with the resultant intersecting dataframe--------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.95)
    figname = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'

    # ---------SFMS from df-------
    p = ax.scatter(df[args.xcol], df[args.ycol], c=df[args.colorcol], marker='s', s=100, lw=1, edgecolor='k')
    cbar = plt.colorbar(p)
    cbar.set_label(label_dict[args.colorcol] if args.colorcol in label_dict else args.colorcol)

    # ---------SFMS from literature-------
    if args.xcol == 'lp_mass' and args.ycol == 'lp_SFR':
        ax = plot_SFMS_Popesso22(ax, 1.8, color='cornflowerblue')
        ax = plot_SFMS_Popesso22(ax, 1.2, color='darkblue')

    # ---------MZR from literature-------
    if args.xcol == 'lp_mass' and 'logOH' in args.ycol:
        ax = plot_MZR_literature(ax)

    # ---------annotate axes and save figure-------
    plt.legend()
    ax.set_xlabel(label_dict[args.xcol] if args.xcol in label_dict else args.xcol)
    ax.set_ylabel(label_dict[args.ycol] if args.ycol in label_dict else args.ycol)

    # --------for talk plots--------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
