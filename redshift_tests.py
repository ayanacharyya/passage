'''
    Filename: redshift_tests.py
    Notes: Compare redshifts given by grizli automatically and the wavelet decomposition, and make plots
    Author : Ayan
    Created: 19-02-26
    Example: run redshift_tests.py --field Par028
             run redshift_tests.py --field Par028 --include_all_lines_in_linelist
             run redshift_tests.py --do_all_fields
             run redshift_tests.py --do_all_fields --include_all_lines_in_linelist
'''

from header import *
from util import *

start_time = datetime.now()
# --------------------------------------------------------------------------------------------------------------
def setup_ax(ax, ratios=None, len_df=None, color='cornflowerblue'):
    # xmin, xmax = list(ax.get_xlim())
    # ymin, ymax = list(ax.get_ylim())

    xmin, xmax = 0, 5.2
    ymin, ymax = 0, 5.2

    if ratios is None: ratios = [1]
    for ratio in ratios:
        xarr = np.array([0, 5])
        yarr1 = (1 + xarr) * ratio - 1
        yarr2 = (1 + xarr) * (1/ratio) - 1
        
        ax.plot(xarr, yarr1, ls='dashed', color=color)
        ax.fill_between(xarr, yarr1 - delta_z, yarr1 + delta_z, alpha=0.2, color=color, zorder=-5)
       
        ax.plot(xarr, yarr2, ls='dashed', color=color)
        ax.fill_between(xarr, yarr2 - delta_z, yarr2 + delta_z, alpha=0.2, color=color, zorder=-5)

        ax.plot([xmin, xmax], [xmin, xmax], c='k', lw=0.5)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.text(0.05, 0.95, f'{args.field}', ha='left', va='top', transform=ax.transAxes)
    if len_df is not None:
        ax.text(0.05, 0.9, f'#obj = {len_df}', ha='left', va='top', transform=ax.transAxes)
    
    return ax

# -------------------------------------------------------------------------------------------------------------
def curtail_df_along_trace(df, xcol, ycol, lines, inside=True):
    for this_line in lines:
        lambda_ratio = lines_dict[this_line] / (rest_wave_dict[strongest_line] * 10)
        expected_z_this_line = (1 + df[xcol]) * lambda_ratio - 1
        df[f'deltaz_{this_line}'] = np.abs(df[ycol] - expected_z_this_line)
    
    target_cols = [f'deltaz_{this_line}' for this_line in lines]
    df['good_flag'] = (df[target_cols] < delta_z).any(axis=1).astype(int)
    df = df[df['good_flag'].astype(bool) == inside]

    return df

# -------------------------------------------------------------------------------------------------------------
def identify_line_match(row, labels, rest_waves, res):
    '''
    Computes expected observed wavelengths based on grizli redshifts
    '''
    obs_waves = np.fromstring(row['obs_wave_arr'], sep=',')
    snrs = np.fromstring(row['snr_arr'], sep=',')
    expected_waves = rest_waves * (1 + row[redshift_col_x])

    diff_matrix = np.abs(expected_waves[:, np.newaxis] - obs_waves)
    allowed_deltas = obs_waves / res
    norm_diff = diff_matrix / allowed_deltas[np.newaxis, :]

    # min_idx[0] = index of the line label
    # min_idx[1] = index of the observed wavelength array
    min_idx = np.unravel_index(np.argmin(norm_diff), norm_diff.shape)
    best_norm_val = norm_diff[min_idx]

    if best_norm_val < 1.0:
        matched_label = labels[min_idx[0]]
        matched_obs_wave = obs_waves[min_idx[1]] # The specific wavelength found
        matched_snr = snrs[min_idx[1]] # The specific snr found
        return pd.Series([True, matched_label, matched_obs_wave, matched_snr])
    else:
        return pd.Series([False, None, np.nan, np.nan])
    
# --------------------------------------------------------------------------------------------------------------
def get_redshift_from_strongest_line(row, strongest_line):
    '''
    Computes redshift assuming the given line label is the strongest line
    '''
    waves = np.fromstring(row['obs_wave_arr'], sep=',')
    snrs = np.fromstring(row['snr_arr'], sep=',')
    
    best_idx = np.argmax(snrs)
    
    strongest_wave = waves[best_idx]
    
    return (strongest_wave / lines_dict[strongest_line]) - 1
    
# --------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.field_arr

    strongest_line = 'Ha'
    delta_z = 0.05
    res = 200
    lines_dict = {'Ha':6563, 'OIII':5007, 'OII':3727, 'MgII':2802, 'CIII':1907}
   
    redshift_col_x = 'grizli_redshift'
    redshift_col_y0 = f'z_linelist_{strongest_line}'
    redshift_col_y1 = f'z_linelist_multi_line'
    redshift_col_y2 = 'zbest'
    redshift_col_y3 = 'zpassagelinefinder'
    delta_z_col = 'delta_z'
    zwidth_col = 'zwidth2'

    colorcol1 = 'zwidth2'
    colorcol2 = f'emline_flag'
    
    lines_for_snr_cut = [strongest_line]# ['Ha', 'OIII', 'OII']
    col_arr = ['salmon', 'cornflowerblue', 'olive', 'teal', 'crimson', 'k']
    snr_values = np.arange(0, 20, 1)
    marker_dict = {'Par028': 'o', 'Par052': 'd', 'Par024': 'D', 'Par025': 's', 'Par026': '^'}
    lw_dict = {'Par028': 2, 'Par052': 2, 'Par024': 1, 'Par025': 1, 'Par026': 1}
    ls_arr = ['solid', 'dashed', 'dotted']

    # ------------finding all available original Par speccats-------------------
    orig_speccat_dir = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5') / 'orig_speccats'
    new_speccat = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_output/v0.5/catalogs/passage_cosmos_redshift_catalog_v2.dat')
    
    fig_dir = Path('/Users/acharyya/Work/astro/passage/passage_output/v0.5') / 'plots' / 'grizli_redshift_comp'
    fig_dir.mkdir(exist_ok=True, parents=True)
    
    if args.do_all_fields:
        files = glob.glob(str(orig_speccat_dir) + '/*speccat.fits')
        args.field_arr = [Path(item).stem[:6] for item in files]
        args.field_arr = [item for item in args.field_arr if item != 'Par052']

    # ---------------------setting up multi-field figure--------------------
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), layout='constrained')
    ax_twinx = ax.twinx()
    incl_text = '_all_ll' if args.include_all_lines_in_linelist else '_strongest_ll'

    # -----------------getting restframe wavelengths-----------------
    lines = list(lines_dict.keys())
    possible_wavelengths = np.array([lines_dict[item] for item in lines])
    possible_wavelength_ratios = [item / (rest_wave_dict[strongest_line] * 10) for item in possible_wavelengths]

    # -----------define filenames----------------
    for ind, args.field in enumerate(args.field_arr):
        print(f'\nDoing ({ind + 1}/{len(args.field_arr)}) field {args.field}..')
        orig_speccat = orig_speccat_dir / f'{args.field}_speccat.fits'

        new_photcat = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.5') / args.field / 'Products' / f'{args.field}_photcat.fits'
        new_linelist = Path('/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.5') / args.field / f'{args.field}lines.dat'

        # ------------read in speccat and photcat files---------------------
        df_orig_speccat = Table.read(orig_speccat)
        df_orig_speccat.remove_column('cdf_z')
        df_orig_speccat = df_orig_speccat.to_pandas()
        df_orig_speccat = df_orig_speccat.rename(columns={'id': 'id_orig', 'redshift': 'grizli_redshift'})

        df_new_photcat = Table.read(new_photcat).to_pandas()
        df_new_photcat = df_new_photcat.rename(columns={'id': 'id_new'})

        df_new_speccat = pd.read_csv(new_speccat, header=0, sep='\t')
        df_new_speccat = df_new_speccat.rename(columns={'field_id': 'id_new'})
        df_new_speccat = df_new_speccat[df_new_speccat['field'] == args.field]

        # ------------read in linelist file---------------------
        df_new_linelist = pd.read_csv(new_linelist, delim_whitespace=True, names=['field_dummy', 'filter', 'id_new', 'obs_wave', 'line_id', 'snr'])
        #df_new_linelist['snr']= df_new_linelist['snr'] / np.sqrt(3)
        
        if args.include_all_lines_in_linelist:
            print(f'Making modified line list, this can take several seconds..')
            df_dummy = pd.DataFrame(columns=['field', 'id_new', 'filter_arr', 'obs_wave_arr', 'line_id_arr', 'snr_arr'])
            for obj in pd.unique(df_new_linelist['id_new']):
                df_sub = df_new_linelist[df_new_linelist['id_new'] == obj]
                df_dummy.loc[len(df_dummy)] = [args.field, obj, 
                                               ','.join(map(str, df_sub['filter'])), 
                                               ','.join(map(str, df_sub['obs_wave'])),
                                               ','.join(map(str, df_sub['line_id'])),
                                               ','.join(map(str, df_sub['snr'])),
                                               ]
            df_new_linelist = df_dummy.copy()
        else:
            # -----------determing redshift from only strongest line in new linelist---------------
            idx = df_new_linelist.groupby('id_new')['snr'].idxmax()
            df_new_linelist = df_new_linelist.loc[idx].reset_index(drop=True)
            df_new_linelist[f'z_linelist_{strongest_line}'] = (df_new_linelist['obs_wave'] / (rest_wave_dict[strongest_line] * 10)) - 1
        
        df_new_linelist = pd.merge(df_new_linelist, df_new_photcat[['id_new', 'ra', 'dec']], on='id_new', how='left')

        # -----------cross matching to new linelist---------------
        df_match = get_crossmatch(df_orig_speccat, df_new_linelist, sep_threshold=0.1, df1_idcol='id_orig', df2_idcol='id_new')
        df_match = df_match.rename(columns={'df1_id':'id_orig', 'df2_id':'id_new'})
        df_match = pd.merge(df_match, df_orig_speccat, on='id_orig', how='inner')
        df_match = pd.merge(df_match, df_new_linelist, on='id_new', how='inner')

        # ------------getting potential match to ANY line in the line list------------------
        if args.include_all_lines_in_linelist:
            df_match[f'z_linelist_{strongest_line}'] = df_match.apply(get_redshift_from_strongest_line, axis=1, args=(strongest_line,))

            df_match[['line_match_flag', 'line_match_id', 'obs_wave_match', 'snr_match']] = df_match.apply(identify_line_match, axis=1, args=(lines, possible_wavelengths, res))
            df_match_good = df_match[df_match['line_match_flag'] == True]
            df_match_good['snr'] = df_match_good['snr_match']
            these_rest_waves = df_match_good['line_match_id'].map(lines_dict)
            df_match_good[redshift_col_y1] =  (df_match_good['obs_wave_match'] / these_rest_waves) - 1
        else:
            # -------------getting objects along the trace------------
            redshift_col_y1 = redshift_col_y0
            df_match_good = curtail_df_along_trace(df_match, redshift_col_x, redshift_col_y1, lines)
            #df_match_good = df_match

        # ----------------matching the good galaxies in orig to new speccat-----------------
        df_match_good_speccat = pd.merge(df_match_good, df_new_speccat, on='id_new', how='inner')

        df_match_good_speccat_disagree = curtail_df_along_trace(df_match_good_speccat, redshift_col_x, redshift_col_y2, [strongest_line], inside=False)

        # ---------making plots--------------
        ax0 = df_match.plot.scatter(redshift_col_x, redshift_col_y0, c=colorcol1)
        ax0 = setup_ax(ax0, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match))
        save_fig(ax0.figure, fig_dir, f'{args.field}{incl_text}_linelist_{redshift_col_x}_vs_{redshift_col_y0}.png', args)

        ax1 = df_match_good.plot.scatter(redshift_col_x, redshift_col_y1, c=colorcol1)
        ax1 = setup_ax(ax1, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match_good))
        save_fig(ax1.figure, fig_dir, f'{args.field}{incl_text}_linelist_agree_{redshift_col_x}_vs_{redshift_col_y1}.png', args)

        ax2 = df_match_good.plot.scatter(redshift_col_x, redshift_col_y1, c=colorcol1)
        ax2 = setup_ax(ax2, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match_good_speccat))
        ax2.scatter(df_match_good_speccat[redshift_col_x], df_match_good_speccat[redshift_col_y1], facecolor='none', edgecolor='green', lw=2, s=50, zorder=10)
        save_fig(ax2.figure, fig_dir, f'{args.field}{incl_text}_linelist_agree_overplotted_{redshift_col_x}_vs_{redshift_col_y1}.png', args)
        
        df_match_good_speccat = df_match_good_speccat[df_match_good_speccat[colorcol2] > 0]
        ax4 = df_match_good_speccat.plot.scatter(redshift_col_x, redshift_col_y2, c=colorcol2, cmap='rainbow')
        ax4 = setup_ax(ax4, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match_good_speccat))
        save_fig(ax4.figure, fig_dir, f'{args.field}{incl_text}_linefinding_{redshift_col_x}_vs_{redshift_col_y2}.png', args)

        ax5 = df_match_good_speccat.plot.scatter(redshift_col_x, redshift_col_y3, c=colorcol2, cmap='rainbow')
        ax5 = setup_ax(ax5, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match_good_speccat))
        save_fig(ax5.figure, fig_dir, f'{args.field}{incl_text}_linefinding_{redshift_col_x}_vs_{redshift_col_y3}.png', args)

        ax6 = df_match_good_speccat_disagree.plot.scatter(redshift_col_x, redshift_col_y3, c=colorcol2, cmap='rainbow')
        ax6 = setup_ax(ax6, ratios=possible_wavelength_ratios, color=col_arr[ind], len_df=len(df_match_good_speccat_disagree))
        save_fig(ax6.figure, fig_dir, f'{args.field}{incl_text}_linefinding_grizli_disagree_{redshift_col_x}_vs_{redshift_col_y3}.png', args)

        # ---------------------making SNR cut plot---------------------        
        print(f'{len(df_match)} galaxies in initial linelist, out of which {len(df_match_good)} are along a locus,\n'\
            f'{len(df_match_good_speccat)} out of {len(df_match)} galaxies made it through inspection without any cut')

        n_fracs_purity, n_fracs_completeness, n_fracs_disagree = [], [], []
        index = 0
        for snr in snr_values:
            df_match_good_sub = df_match_good[df_match_good['snr'] > snr]
            df_sub = pd.merge(df_match_good_sub, df_new_speccat, on='id_new', how='inner')
            n_fracs_purity.append(len(df_sub) / len(df_match_good_sub))
            
            df_sub = df_match_good_speccat[df_match_good_speccat['snr'] > snr]
            n_fracs_completeness.append(len(df_sub) / len(df_new_speccat))

            df_sub = df_match_good_speccat_disagree[df_match_good_speccat_disagree[f'snr'] > snr]
            n_fracs_disagree.append(len(df_sub) / len(df_match_good_speccat))

        ax.plot(snr_values, n_fracs_purity, c=col_arr[ind], ls='solid', lw=lw_dict[args.field], label=f'{args.field}: Purity' if ind == 0 else f'{args.field}' if index == 0 else None)
        ax.plot(snr_values, n_fracs_completeness, c=col_arr[ind], ls='dashed', lw=lw_dict[args.field], label=f'{args.field}: Completeness' if ind == 0 else f'{args.field}' if index == 0 else None)
        
        ax_twinx.plot(snr_values, n_fracs_disagree, c=col_arr[ind], ls='dotted', lw=lw_dict[args.field], label=f'{args.field}: Disagreeing' if ind == 0 and index == 0 else None)

    # ---------annotate and save SNR cut plot-------------------
    ax.legend(loc='upper right')
    ax.set_xlabel(f'Linelist SNR (threshold)')
    ax.set_ylabel('Fraction of objects')
    ax.set_ylim(0, 1)

    ax_twinx.legend(loc='lower center')
    ax_twinx.set_ylabel('Fraction of disagreeing redshifts')
    ax_twinx.set_ylim(0, 0.3)

    save_fig(fig, fig_dir, f'{",".join(args.field_arr)}{incl_text}_grizli_linelist_zcomp_survived_frac_vs_snr.png', args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
