'''
    Filename: get_targets_cy5.py
    Notes: Plots location and best-fit SED models of galaxies selected in a certain way, for JWST Cycle 5 proposal submission
    Author : Ayan
    Created: 29-09-25
    Example: run get_targets_cy5.py
'''

from header import *
from util import *
from plot_footprints import plot_background, plot_footprints, plot_skycoord_from_df

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------directory structures and file names-------------
    reg_files_dir = args.input_dir / 'footprints/region_files'

    if args.bg_image_dir is None: args.bg_image_dir = args.input_dir / 'footprints/survey_images'
    else: args.bg_image_dir = Path(args.bg_image_dir)

    if args.fg_image_dir is None: args.fg_image_dir = args.input_dir / 'COSMOS/imaging'
    else: args.fg_image_dir = Path(args.fg_image_dir)

    if args.bg_file is None:
        if args.bg == 'COSMOS': bg_filename = args.bg_image_dir / 'COSMOS-HST-ACS_mosaic_Shrink100.fits' # if field is COSMOS, by default use this background image
        else: bg_filename = list(args.bg_image_dir.glob('*%s*.fits' %args.bg))[0] # for other fields, look for available background image files
    else:
        bg_filename = args.bg_image_dir / args.bg_file

    reg_filenames = []
    reg_filenames += list(reg_files_dir.glob('*%s*NIRCam*.reg' %args.bg))
    reg_filenames += list(reg_files_dir.glob('*%s*WFC3+ACS*.reg' %args.bg))

    candels_regfile = args.input_dir / 'footprints/region_files/CANDELS_COSMOS_WFC3+ACS.reg'
    merged_filename =  args.input_dir / 'COSMOS' / 'zCOSMOS_WFC3+ACS_Web_2020.csv'

    zcosmos_idcol, cweb_idcol, c2020_idcol = 'ID_zcosmos', 'ID_cweb', 'ID_c2020'

    # ------plotting the background------------
    fig, bg_img_hdu = plot_background(bg_filename, args)
    wcs_header = pywcs.WCS(bg_img_hdu[0].header)
    fig, sky_regions = plot_footprints(reg_filenames, bg_img_hdu, fig, args)

    # -----reading the selected objects------
    if merged_filename.exists() and not args.clobber:
        print(f'Reading merged catalog from existing {merged_filename}')    
    else:
        print(f'Did not find {merged_filename}; so making new..')
        # -----reading zcosmos catalog and filtering---------
        df = read_zCOSMOS_catalog(args=args)
        candels_regions = Regions.read(candels_regfile, format='ds9')
        print(f'Cross-matching zCOSMOS catalog and WFC3+ACS region...')
        df_candels = df[candels_regions[1].contains(SkyCoord(df['ra'], df['dec'], unit='deg'), wcs_header)]
        print(f'...found {len(df_candels)} objects in this region, now applying z-cut off..')
        df_candels_zsub = df_candels[df_candels['REDSHIFT'].between(1.0,1.7)]
        df_candels_zsub.rename(columns={'REDSHIFT':'redshift', 'id':zcosmos_idcol}, inplace=True)
        print(f'...found {len(df_candels_zsub)} after the z-cut\n')

        # --------load COSMOS-Web catalog and crossmatch---------
        df_cweb = read_COSMOSWebb_catalog(args=args)
        df_cweb.rename(columns={'id':cweb_idcol}, inplace=True)
        print(f'Cross-matching zCOSMOS-sub and COSMOS-Web catalogs; can take a while...')
        df_candels_cweb = get_crossmatch(df_candels_zsub, df_cweb, sep_threshold=0.1, df1_idcol=zcosmos_idcol, df2_idcol=cweb_idcol)
        print(f'...found {len(df_candels_cweb)} cross-matched objects\n')

        df_candels_cweb.rename(columns={'df1_id':zcosmos_idcol, 'df2_id':cweb_idcol}, inplace=True)
        df_candels_cweb = pd.merge(df_candels_cweb, df_candels_zsub, on=zcosmos_idcol, how='inner')
        
        df_cweb.rename(columns={'ra':'ra_cweb', 'dec':'dec_cweb'}, inplace=True)
        df_candels_cweb = pd.merge(df_candels_cweb, df_cweb, on=cweb_idcol, how='inner')

        # ------load COSMOS 2020 catalog and crossmatch---------
        df_c2020 = read_COSMOS2020_catalog(args=args)
        df_c2020.rename(columns={'id':c2020_idcol}, inplace=True)
        print(f'Cross-matching zOOSMOS-sub+COSMOS-Web and COSMOS-2020 catalogs; can take a while...')
        df_candels_cweb_c2020 = get_crossmatch(df_candels_cweb, df_c2020, sep_threshold=0.1, df1_idcol=cweb_idcol, df2_idcol=c2020_idcol)
        print(f'...found {len(df_candels_cweb_c2020)} cross-matched objects\n')

        df_candels_cweb_c2020.rename(columns={'df1_id':cweb_idcol, 'df2_id':c2020_idcol}, inplace=True)
        df_candels_cweb_c2020 = pd.merge(df_candels_cweb_c2020, df_candels_cweb, on=cweb_idcol, how='inner')

        df_c2020.rename(columns={'ra':'ra_c2020', 'dec':'dec_c2020'}, inplace=True)
        df_candels_cweb_c2020 = pd.merge(df_candels_cweb_c2020, df_c2020, on=c2020_idcol, how='inner')
  
        # --------modifying columns for SED fitting-------------
        df_candels_cweb_c2020.drop('ID', axis=1, inplace=True)
        df_candels_cweb_c2020.drop(columns=df_candels_cweb_c2020.columns[df_candels_cweb_c2020.columns.str.contains('FLAG')], inplace=True)
        df_candels_cweb_c2020.drop(columns=df_candels_cweb_c2020.columns[df_candels_cweb_c2020.columns.str.contains('SPLASH')], inplace=True)
        
        for thiscol in df_candels_cweb_c2020.columns:
            if 'FLUX_APER_' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[10:] + '_sci'}, inplace=True)
            elif 'FLUX_ERR_APER_' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[14:] + '_err'}, inplace=True)
            elif '_FLUX_AUTO' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[:-10] + '_sci'}, inplace=True)
            elif '_FLUXERR_AUTO' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[:-13] + '_err'}, inplace=True)
            elif '_FLUXERR' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[:-8] + '_err'}, inplace=True)
            elif '_FLUX' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol[:-5] + '_sci'}, inplace=True)

            if 'F115W' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('F115W', 'NIRCAM_F115W')}, inplace=True)
            elif 'F150W' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('F150W', 'NIRCAM_F150W')}, inplace=True)
            elif 'F277W' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('F277W', 'NIRCAM_F277W')}, inplace=True)
            elif 'F444W' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('F444W', 'NIRCAM_F444W')}, inplace=True)
            elif 'F770W' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('F770W', 'MIRI_F770W')}, inplace=True)
            elif 'IA' in thiscol: df_candels_cweb_c2020.rename(columns={thiscol: thiscol.replace('IA', 'IB')}, inplace=True)

        # -----saving the final catalog--------
        df_candels_cweb_c2020.to_csv(merged_filename, index=None)
        print(f'Saved merged, z-filtered catalog as {merged_filename}')

    # ------plotting the targets on sky--------
    df_targets = pd.read_csv(merged_filename)
    ax = fig.gca()
    for index, row in df_targets.iterrows():
        circle = plt.Circle(xy=(row['ra'], row['dec']), radius=50/3600, color='g', alpha=1, fill=True, transform=ax.get_transform('fk5'))
        ax.add_patch(circle)

    # --------saving the plot----------
    figname = args.input_dir / 'footprints' / f'{args.bg}_with_footprints_with_zCOSMOS_sub.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
