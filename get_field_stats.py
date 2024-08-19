'''
    Filename: get_field_stats.py
    Notes: Computes and plots various statistics (redshift distribution, detected line distribution, magnitude distribution etc.) for a given PASSAGE field which has been already run through make_diagnostic_maps.py
    Author : Ayan
    Created: 19-08-24
    Example: run get_field_stats.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --re_extract
             run get_field_stats.py --field Par61
'''
from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def get_detection_fraction(df, line, args, EW_thresh=300):
    '''
    To compute fraction of objects in a given field that have a given line detected (beyond a given EW threshold)
    out of all objects where the line was accessible in the given filter wavelength regime
    Prints out results, and returns subset dataframe
    '''
    zrange_arr = get_zranges_for_filters(line, filters=args.filters)[0]

    dfline = pd.DataFrame()
    for zrange in zrange_arr:
        dfsub = df[(df['redshift'].between(zrange[0], zrange[1]))]
        print(f'{len(dfsub)} objects within z range {zrange}')
        dfline = pd.concat([dfline, dfsub])

    dfline_detected = dfline[dfline[f'{line}_EW'] > EW_thresh]

    print(f'{len(dfline_detected)} out of {len(dfline)} (out of total {len(df)}) has strong {line} detections!')

    # ------plot histogram of line EWs--------
    fig, ax = plt.subplots()
    ax.hist(np.log10(dfsub[f'{line}_EW']), bins=50)
    ax.axvline(np.log10(EW_thresh), c='k', lw=2, ls='dashed')

    ax.set_ylabel('#objects')
    ax.set_xlabel(f'log {line} EW (A)')
    ax.text()

    figname = args.output_dir / f'{args.field}' / f'{args.field}_{line}_EW_histogram.png'
    fig.savefig(figname)
    print(f'Saved figure as {figname}')
    plt.show(block=False)

    return dfline_detected


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    EW_thresh = 300

    # ---------determining filename suffixes-------------------------------
    output_dir = args.output_dir / args.field
    if args.re_extract: output_dir = output_dir / 're_extracted'
    df_filename = output_dir / f'{args.field}_all_diag_results.txt'
    df = pd.read_table(df_filename, delim_whitespace=True)

    # ----------------------initiliasing dataframe-------------------------------
    all_lines_to_consider = ['OII', 'OIII'] # OR args.line_list #

    # -------------------------------------------------------
    for line in all_lines_to_consider: dfline_detected = get_detection_fraction(df, line, args, EW_thresh=EW_thresh)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
