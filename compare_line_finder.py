'''
    Filename: compare_line_finder.py
    Notes: Compares fluxes and redshifts measurments from various sources, given a spreadsheet
    Author : Ayan
    Created: 22-11-24
    Example: run compare_line_finder.py
'''
from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ---------reading in dataframe from google sheets----------------
    link = 'https://docs.google.com/spreadsheets/d/15fDDGHmLugFD-fo1pI2SYzbuRjt0BEPDsw8zSJRBijw/edit?gid=0#gid=0'
    ghseet_id = link[link.find('/d/') + 3: link.find('/edit')]
    url = f'https://docs.google.com/spreadsheet/ccc?key={ghseet_id}&output=xlsx'
    df = pd.read_excel(url, sheet_name='Sheet2')
    df = df.replace('-', np.nan).replace('x', np.nan)
    print(f'Read in spreadsheet from url {url}')

    # -------determining how many panels to plot-------------------
    example = 'Ha'
    sources = [item[3:] for item in df.columns if example in item]
    if 'G1D' in sources and 'z_G' in df: df.rename(columns={'z_G':'z_G1D'}, inplace=True)
    if 'G2D' in sources: df['z_G2D'] = np.nan
    quantities = [item[:-len(sources[0])-1] for item in df.columns if sources[0] in item]
    df['index'] = df.index

    # --------plotting figure---------------
    fig, axes = plt.subplots(len(quantities), 1, figsize=(8, 7.5), sharex=True)
    fig.subplots_adjust(right=0.98, top=0.98, bottom=0.07, left=0.1, hspace=0.1)
    color_dict = {'G1D':'brown', 'G2D':'salmon', 'AA':'lightblue', 'KN':'cornflowerblue', 'FH':'darkblue'}

    for index, quantity in enumerate(quantities):
        print(f'Doing quantity {index + 1} of {len(quantities)}..')
        quant_cols = [item for item in df.columns if quantity + '_' in item]
        df_sub = df[np.hstack(['index', quant_cols])]
        df_sub[quantity + '_min'] = df_sub[quant_cols].min(axis=1)
        for col in quant_cols: df_sub[col] = df_sub[col] - df_sub[quantity + '_min'] # normalise by the minimum measurement found by any method
        for source in sources: axes[index].plot(df_sub['index'], df_sub[quantity + '_' + source], marker='o', lw=1, c=color_dict[source], label=source)
        if index == 0: axes[index].legend(fontsize=args.fontsize/1.5, loc='best')
        axes[index].set_ylabel(r'$\Delta$ ' + quantity, fontsize=args.fontsize)
        #axes[index].set_ylim(-0.1, 1)

    # ---------annotate axes and save figure-------
    figname = args.output_dir / 'plots' / f'line_finder_comparison.png'
    axes[-1].set_xlabel(r'Galaxies', fontsize=args.fontsize)
    axes[-1].set_xticklabels([f'#{df.iloc[int(item)]["objid"]}' if item >=0 and item < len(df) else '' for item in axes[-1].get_xticks()], fontsize=args.fontsize)

    # --------for talk plots--------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try:
            mplcyberpunk.make_lines_glow()
        except:
            pass
        try:
            mplcyberpunk.make_scatter_glow()
        except:
            pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
