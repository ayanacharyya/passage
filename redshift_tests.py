'''
    Filename: redshift_tests.py
    Notes: Compare redshifts given by grizli automatically and the wavelet decomposition, and make plots
    Author : Ayan
    Created: 19-02-26
    Example: run redshift_tests.py --field Par028
'''

from header import *
from util import *

start_time = datetime.now()
# --------------------------------------------------------------------------------------------------------------
def setup_ax(ax, ratios=None):
    xmin, xmax = list(ax.get_xlim())
    ymin, ymax = list(ax.get_ylim())

    if ratios is None: ratios = [1]
    for ratio in ratios:
        xarr = np.array([0, 5])
        yarr = (1 + xarr) * ratio - 1
        ax.plot(xarr, yarr, ls='dashed', color='cornflowerblue')
        ax.fill_between(xarr, yarr - delta_z, yarr + delta_z, alpha=0.2, color='cornflowerblue', zorder=-5)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    return ax

# -------------------------------------------------------------------------------------------------------------
def curtail_df_along_trace(df, xcol, ycol):
    for this_line in list(lines_dict.keys()):
        lambda_ratio = lines_dict[this_line] / (rest_wave_dict[line] * 10)
        expected_z_this_line = (1 + df[xcol]) * lambda_ratio - 1
        df[f'deltaz_{this_line}'] = np.abs(df[ycol] - expected_z_this_line)
    
    target_cols = [f'deltaz_{this_line}' for this_line in list(lines_dict.keys())]
    df['good_flag'] = (df[target_cols] < delta_z).any(axis=1).astype(int)
    df = df[df['good_flag'] == 1]

    return df

# --------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    line = 'Ha'
    delta_z = 0.05
    lines_dict = {'Ha':6563, 'OIII':5007, 'OII':3727, 'MgII':2802, 'CIII':1907}
   
    # -----------define filenames----------------
    orig_speccat = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5') / args.field / 'Products' / f'{args.field}_speccat.fits'

    new_photcat = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.5') / args.field / 'Products' / f'{args.field}_photcat.fits'
    new_linelist = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.5') / args.field / f'{args.field}lines.dat'
    #new_linelist = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.5') / args.field / f'{args.field}lines_modified.dat'

    new_speccat = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_output/v0.5/catalogs/passage_cosmos_redshift_catalog_v2.dat')
    
    # ------------read in speccat and photcat files---------------------
    df_orig_speccat = Table.read(orig_speccat)
    df_orig_speccat.remove_column('cdf_z')
    df_orig_speccat = df_orig_speccat.to_pandas()
    df_orig_speccat = df_orig_speccat.rename(columns={'id': 'id_orig'})

    df_new_photcat = Table.read(new_photcat).to_pandas()
    df_new_photcat = df_new_photcat.rename(columns={'id': 'id_new'})

    df_new_speccat = pd.read_csv(new_speccat, header=0, sep='\t')
    df_new_speccat = df_new_speccat.rename(columns={'field_id': 'id_new'})
    df_new_speccat = df_new_speccat[df_new_speccat['field'] == args.field]

    # ------------read in linelist file---------------------
    if 'modified' in str(new_linelist):
        df_new_linelist = pd.read_csv(new_linelist, delim_whitespace=True)
    else:
        df_new_linelist = pd.read_csv(new_linelist, delim_whitespace=True, names=['field_dummy', 'filter', 'id_new', 'obs_wave', 'line_id', 'snr'])
        # -----------determing redshift from only strongest line in new linelist---------------
        idx = df_new_linelist.groupby('id_new')['snr'].idxmax()
        df_new_linelist = df_new_linelist.loc[idx].reset_index(drop=True)
        df_new_linelist = pd.merge(df_new_linelist, df_new_photcat[['id_new', 'ra', 'dec']], on='id_new', how='left')
        df_new_linelist[f'z_linelist_default_{line}'] = (df_new_linelist['obs_wave'] / (rest_wave_dict[line] * 10)) - 1    

    # -----------cross matching to new linelist---------------
    df_match1 = get_crossmatch(df_orig_speccat, df_new_linelist, sep_threshold=0.1, df1_idcol='id_orig', df2_idcol='id_new')
    df_match1 = df_match1.rename(columns={'df1_id':'id_orig', 'df2_id':'id_new'})
    df_match1 = pd.merge(df_match1, df_orig_speccat, on='id_orig', how='inner')
    df_match1 = pd.merge(df_match1, df_new_linelist, on='id_new', how='inner')

    # -----------cross matching to new speccat---------------
    # df_match2 = get_crossmatch(df_orig_speccat, df_new_speccat, sep_threshold=0.1, df1_idcol='id_orig', df2_idcol='id_new')
    # df_match2 = df_match2.rename(columns={'df1_id':'id_orig', 'df2_id':'id_new'})
    # df_match2 = pd.merge(df_match2, df_orig_speccat, on='id_orig', how='inner')
    # df_match2 = pd.merge(df_match2, df_new_speccat, on='id_new', how='inner')

    # ------------computing new columns------------------
    redshift_col_x = 'redshift'
    redshift_col_y1 = f'z_linelist_default_{line}'
    redshift_col_y2 = 'zpassagelinefinder'
    delta_z_col = 'delta_z'
    zwidth_col = 'zwidth2'

    colorcol1 = 'zwidth2'
    colorcol2 = f'snr_{line}'

    df_match1[delta_z_col] = (df_match1[redshift_col_x] - df_match1[redshift_col_y1]) / (df_match1[redshift_col_x] + 1)
    # df_match2[delta_z_col] = (df_match2[redshift_col_x] - df_match2[redshift_col_y2]) / (df_match2[redshift_col_x] + 1)

    # -------------getting objects along the trace------------
    possible_wavelengths = [lines_dict[item] for item in list(lines_dict.keys())]
    possible_wavelength_ratios = [item / (rest_wave_dict[line] * 10) for item in possible_wavelengths]

    # -----------determing redshift from all lines in new linelist---------------
    df_match1_good = curtail_df_along_trace(df_match1, redshift_col_x, redshift_col_y1)
    # df_match2_good = curtail_df_along_trace(df_match2, redshift_col_x, redshift_col_y2)

    # ----------------matching the good galaxies in orig to new speccat-----------------
    df_match1_good_speccat = pd.merge(df_match1_good, df_new_speccat, on='id_new', how='inner')
    df_match1_good_speccat[delta_z_col] = (df_match1_good_speccat[redshift_col_x] - df_match1_good_speccat[redshift_col_y2]) / (df_match1_good_speccat[redshift_col_x] + 1)

    # ---------making plots--------------
    ax0 = df_match1.plot.scatter(redshift_col_x, redshift_col_y1, c=colorcol1)
    ax0 = setup_ax(ax0, ratios=possible_wavelength_ratios)

    ax1 = df_match1_good.plot.scatter(redshift_col_x, redshift_col_y1, c=colorcol1)
    ax1 = setup_ax(ax1, ratios=possible_wavelength_ratios)

    ax2 = df_match1_good.plot.scatter(redshift_col_x, redshift_col_y1, c=colorcol1)
    ax2 = setup_ax(ax2, ratios=possible_wavelength_ratios)
    ax2.scatter(df_match1_good_speccat[redshift_col_x], df_match1_good_speccat[redshift_col_y1], facecolor='none', edgecolor='green', lw=2, s=50, zorder=10)

    # ax3 = df_match2_good.plot.scatter(redshift_col_x, redshift_col_y2, c=colorcol2)
    # ax3 = setup_ax(ax3, ratios=possible_wavelength_ratios)
    
    ax4 = df_match1_good_speccat.plot.scatter(redshift_col_x, redshift_col_y2, c=colorcol2)
    ax4 = setup_ax(ax4, ratios=possible_wavelength_ratios)

    # ---------------------testing with SNR cuts---------------------
    lines = ['Ha', 'OIII', 'OII']
    col_arr = ['salmon', 'cornflowerblue', 'olive']
    snr_values = np.arange(0, 20, 1)
    n_survive = len(df_match1_good_speccat)

    print(f'{len(df_match1)} galaxies in initial linelist, out of which {len(df_match1_good)} are along a locus,\n'\
          #f'{len(df_match2)} galaxies in inspected linelist, out of which {len(df_match2_good)} are along a locus,\n'\
          f'{n_survive} out of {len(df_match1)} galaxies made it through inspection without any cut')

    # ---------------------making SNR cut plots---------------------
    fig, ax = plt.subplots()
    ax.scatter(0, n_survive / len(df_match1_good), marker='*', s=70, c='gray')
    
    for index, line in enumerate(lines):
        n_fracs = []
        for snr in snr_values:
            df_match1_good_sub = df_match1_good[df_match1_good[f'sn_{line}'] > snr]
            df_match1_good_speccat_sub = pd.merge(df_match1_good_sub, df_new_speccat, on='id_new', how='inner')
            n_fracs.append(len(df_match1_good_speccat_sub) / len(df_match1_good_sub))
 
        ax.scatter(snr_values, n_fracs, c=col_arr[index], edgecolor='k', lw=1, s=50, label=f'{line} SNR')

    ax.legend()
    ax.set_xlabel(f'SNR threshold')
    ax.set_ylabel('Fraction of galaxies above given SNR')



    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
