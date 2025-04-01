'''
    Filename: make_glass_speccat.py
    Notes: Reads in all existing full.fits file from GLASS-NIRISS (from PJW) and makes a catalog of the integrated line fluxes from the individual file headers, and plot the redshift distribution
    Author : Ayan
    Created: 27-03-25
    Example: run make_glass_speccat.py --field glass-a2744 --line_list OIII,OII,Hb,NeIII-3867 --SNR_thresh 2
'''
from header import *
from util import *
from make_diagnostic_maps import get_emission_line_int

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.drv = 'v0.5'

    # ---------setting up directory paths-------------------
    which_files = 'full'
    full_dir = args.input_dir / args.field / 'Products' / which_files
    if not os.path.exists(full_dir):
        which_files = 'maps'
        full_dir = args.input_dir / args.field / 'Products' / which_files
    plots_dir = args.output_dir / 'plots'
    catalog_dir = args.output_dir / 'catalogs'
    for this_dir in [plots_dir, catalog_dir]: this_dir.mkdir(exist_ok=True, parents=True)
    catfilename = catalog_dir / f'{args.field}_line_fluxes.csv'

    # ----------setting up dataframe---------------------
    basic_cols = ['field', 'objid', 'ra', 'dec', 'redshift']
    lines_to_save = rest_wave_dict.keys()
    line_cols = np.hstack([[item + '_int', item + '_int_u', item + '_snr'] for item in lines_to_save])
    cols_in_df = np.hstack([basic_cols, line_cols])
    df = pd.DataFrame(columns=cols_in_df)

    # -----------looping over all files---------------------------
    if not os.path.exists(catfilename) or args.clobber:
        print(f'Preparing line flux catalog file..\n')
        full_files = glob.glob(str(full_dir) + f'/*{which_files}.fits')
        for index, file in enumerate(full_files):
            objid = int(Path(file).stem.split('.')[0][-5:])
            print(f'Doing file {objid}, which is {index + 1} of {len(full_files)}')

            full_hdu = fits.open(file)
            header = full_hdu[0].header
            args.available_lines = np.array(header['HASLINES'].split(' '))
            args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines])  # replace 'OIII-5007' with 'OIII'

            basic_cols_quant = [args.field, objid, header['RA'], header['DEC'], header['REDSHIFT']]
            lines_cols_quant = []
            for line in lines_to_save:
                try:
                    line_int = get_emission_line_int(line, full_hdu, args, dered=False)
                    snr = line_int.n / line_int.s
                except:
                    #print(f'Could not find {line} for object {objid}')
                    line_int, snr = ufloat(0, 0), 0
                lines_cols_quant.append([line_int.n, line_int.s, snr])

            this_row = np.hstack(([basic_cols_quant, np.hstack(lines_cols_quant)]))
            df.loc[len(df)] = this_row

        # ------writing out the catalog----------------
        for col in line_cols: df[col] = df[col].astype(np.float64)
        df.to_csv(catfilename, index=None, header='column_names')
        print(f'Saved catalog in {catfilename}')
    else:
        print(f'Reading in existing flux catalog from {catfilename}')
        df = pd.read_csv(catfilename)

    # -------setting up the figure-----------------
    histcol = 'redshift'
    fig, ax = plt.subplots(figsize=(8, 6))

    line_list_array = [args.line_list] # [['OIII', 'OII', 'Hb', 'NeIII-3867'], ['OIII', 'OII', 'Hb', 'Ha', 'SII'], ['OIII', 'OII', 'Hb', 'Ha', 'SII', 'NeIII-3867']]
    col_arr = ['salmon', 'cornflowerblue', 'darkgreen']

    # ----------filtering dataframe-------------------
    for index, args.line_list in enumerate(line_list_array):
        df_sub = df.copy()
        for line in args.line_list:
            print(f'Applying SNR >= {args.SNR_thresh} cut on {line}..')
            df_sub = df_sub[df_sub[f'{line}_snr'] >= args.SNR_thresh]

        line_cols = np.hstack([[item + '_int', item + '_int_u', item + '_snr'] for item in args.line_list])
        df_sub = df_sub[np.hstack([basic_cols, line_cols])].drop(labels=['field'], axis=1)
        print(f'{len(df_sub)} objects remain out of original {len(df)}')


        ax.hist(df_sub[histcol], bins=50, range=[1.5,3.5], color=col_arr[index], histtype='step')
        ax.text(0.99, 0.99 - index * 0.1, f'Total {len(df_sub)} {args.field} objects\nwith {",".join(args.line_list)} SNR > {args.SNR_thresh}', color=col_arr[index], ha='right', va='top', transform=ax.transAxes)

    # -------annotating the figure-----------------
    ax.set_xlabel(histcol, fontsize=args.fontsize)
    ax.set_ylabel('# of objects', fontsize=args.fontsize)

    ax.set_ylim(0, 6)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ------------saving the full figure--------------------------
    figname = plots_dir / f'{args.field}_{histcol}_histogram.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
