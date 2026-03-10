'''
    Filename: analyse_jades_spec.py
    Notes: Reads in JADES NIRSpec spectrsocopic catalog, and plots various quantitites including BPT diagram
    Author : Ayan
    Created: 06-03-26
    Example: run analyse_jades_spec.py --plot_BPT --plot_DIG
             run analyse_jades_spec.py
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def AGN_func(x, theoretical_line):
    '''
    Equation for AGN demarcation line on R3-S2 BPT, from different literature sources
    '''
    if theoretical_line == 'K01': # Eq 5 of Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)
        y = np.where(x < 0.47, 1.19 + 0.61 / (x - 0.47), np.nan)
    elif theoretical_line == 'K13': # Eq 1 of Kauffmann+2013
        y = np.where(x < 0.05, 1.3 + 0.61 / (x - 0.05), np.nan)
    else:
        sys.exit(f'Requested theoreitcal line {theoretical_line} should be one of K01,K13')
    return y

# --------------------------------------------------------------------------------------------------------------
def get_sum(df, line1, line2):
    '''
    Computes the flux sum for given numerator and denominator line labels, appropriately propagating uncertainties
    Returns dataframe with two new (ratio and uncertainty) columns as well as the name of the line ratio
    '''
    sum_name = line1 + '-' + line2
    sum_with_unc = unp.uarray(df[f'{line1}_flux'], df[f'{line1}_err']) + unp.uarray(df[f'{line2}_flux'], df[f'{line2}_err'])
    df[f'{sum_name}_flux'] = unp.nominal_values(sum_with_unc)
    df[f'{sum_name}_err'] = unp.std_devs(sum_with_unc)

    df = df[(~np.isnan(df[f'{sum_name}_flux'])) & (~np.isnan(df[f'{sum_name}_err']))]

    return df, sum_name

# --------------------------------------------------------------------------------------------------------------
def get_log_ratio(df, num_line, den_line):
    '''
    Computes the log flux ratio for given numerator and denominator line labels, appropriately propagating uncertainties
    Returns dataframe with two new (ratio and uncertainty) columns as well as the name of the line ratio
    '''
    ratio_name = num_line.split('_')[0] + den_line.split('_')[0]
    ratio_with_unc = unp.log10(unp.uarray(df[f'{num_line}_flux'], df[f'{num_line}_err']) / unp.uarray(df[f'{den_line}_flux'], df[f'{den_line}_err']))
    df[f'{ratio_name}'] = unp.nominal_values(ratio_with_unc)
    df[f'{ratio_name}_err'] = unp.std_devs(ratio_with_unc)

    df = df[(~np.isnan(df[f'{ratio_name}'])) & (~np.isnan(df[f'{ratio_name}_err']))]

    return df, ratio_name

# --------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------------user defined quantities-------------
    snr_thresh = 2
    #colorcol = 'C3_1907_snr'
    colorcol = 'z_Spec'
    cmap = 'plasma'

    # -----------get filenames----------
    input_dir = Path('/Users/acharyya/Work/astro/jades')
    fig_dir = input_dir / 'plots'
    cat_dir = input_dir / 'catalogs'

    # -------read in catalogs-------------
    files = glob.glob(str(cat_dir) + '/*catalog.fits')
    
    df = pd.DataFrame()
    for file in files:
        this_df = Table.read(file).to_pandas()
        df = pd.concat([df, this_df])

    # ----------modify dataframe--------------
    lines = [item[:-5] for item in df.columns if '_flux' in item]
    for line in lines: df[f'{line}_snr'] = df[f'{line}_flux'] / df[f'{line}_err']

    # ----------------for BPT diagram-------------
    if args.plot_BPT:
        # ----------curtail dataframe-------------
        lines_to_filter_by = ['N2_6584', 'HB_4861']
        dfsub = df.copy()
        for line in lines_to_filter_by: dfsub = dfsub[(dfsub[f'{line}_snr'] > snr_thresh)]

        dfsub, O3Hb = get_log_ratio(dfsub, 'O3_5007', 'HB_4861')
        dfsub, N2Ha = get_log_ratio(dfsub, 'N2_6584', 'HA_6563')
        print(f'{len(dfsub)} of {len(df)} objects survive.')

        # -------plot BPT diagram---------------
        fig, ax = plt.subplots(1, layout='constrained')
        
        p = ax.scatter(dfsub[N2Ha], dfsub[O3Hb], c=dfsub[colorcol], cmap=cmap, s=20, lw=0.1, edgecolor='w' if args.fortalk else 'k', zorder=10)
        ax.errorbar(dfsub[N2Ha], dfsub[O3Hb], xerr=dfsub[f'{N2Ha}_err'], yerr=dfsub[f'{O3Hb}_err'], c='gray', fmt='none', lw=1, alpha=0.8)

        # -----------annotate axes-------------------
        cbar = plt.colorbar(p)
        cbar.set_label(colorcol, fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)

        ax.set_xlim(-2.5, 0.5)
        ax.set_ylim(-1, 1.5)
        ax.set_xlabel(f'Log {N2Ha}', fontsize=args.fontsize)
        ax.set_ylabel(f'Log {O3Hb}', fontsize=args.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

        ax.text(0.95, 0.95, f'{len(dfsub)} objects', c='k', ha='right', va='top', fontsize=args.fontsize, transform=ax.transAxes)

        # ---------adding literature AGN demarcation lines----------
        label_dict = {'K01': 'Kewley+2001', 'K13':'Kauffmann+2013'}
        theoretical_lines = ['K01', 'K13']
        ls_arr = ['solid', 'dashed']
        x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)

        for index, theoretical_line in enumerate(theoretical_lines):
            y = AGN_func(x, theoretical_line)
            ax.plot(x, y, c='k', ls=ls_arr[index], lw=1, label=label_dict[theoretical_line])
        ax.legend(loc='lower left', fontsize=args.fontsize)

        # --------saving the fig-------------
        figname = f'JADES_{N2Ha}_vs_{O3Hb}.png'
        save_fig(fig, fig_dir, figname, args)


    # ---------------for DIG----------------
    if args.plot_DIG:
        # ----------curtail dataframe-------------
        lines_to_filter_by = ['S2_6718', 'S2_6733', 'HA_6563', 'N2_6584']
        dfsub = df.copy()
        for line in lines_to_filter_by: dfsub = dfsub[(dfsub[f'{line}_snr'] > snr_thresh)]

        dfsub, S2sum = get_sum(dfsub, 'S2_6733', 'S2_6718')
        dfsub, O3Hb = get_log_ratio(dfsub, S2sum, 'HA_6563')
        dfsub, N2Ha = get_log_ratio(dfsub, 'N2_6584', 'HA_6563')
        print(f'{len(dfsub)} of {len(df)} objects survive.')

        # -------plot BPT diagram---------------
        fig, ax = plt.subplots(1, layout='constrained')
        
        p = ax.scatter(dfsub[N2Ha], dfsub[O3Hb], c=dfsub[colorcol], cmap=cmap, s=20, lw=0.1, edgecolor='w' if args.fortalk else 'k', zorder=10)
        ax.errorbar(dfsub[N2Ha], dfsub[O3Hb], xerr=dfsub[f'{N2Ha}_err'], yerr=dfsub[f'{O3Hb}_err'], c='gray', fmt='none', lw=1, alpha=0.8)

        # -----------annotate axes-------------------
        cbar = plt.colorbar(p)
        cbar.set_label(colorcol, fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)

        ax.set_xlim(-2.5, 0.5)
        ax.set_ylim(-1.5, 0.5)
        ax.set_xlabel(f'Log {N2Ha}', fontsize=args.fontsize)
        ax.set_ylabel(f'Log {O3Hb}', fontsize=args.fontsize)
        ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

        ax.text(0.95, 0.95, f'{len(dfsub)} objects', c='k', ha='right', va='top', fontsize=args.fontsize, transform=ax.transAxes)

        # --------saving the fig-------------
        figname = f'JADES_{N2Ha}_vs_{O3Hb}.png'
        save_fig(fig, fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
   

