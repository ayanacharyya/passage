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
    args.candels_cut = 0
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

    candels_cut_text = '_candelscut' if args.candels_cut else ''
    candels_regfile = args.input_dir / 'footprints/region_files/CANDELS_COSMOS_WFC3+ACS.reg'
    #merged_filename = args.input_dir / 'COSMOS' / f'zCOSMOS_Web_2020{candels_cut_text}.fits'
    merged_filename = args.input_dir / 'COSMOS' / f'zCOSMOS_2025{candels_cut_text}.fits'
    photcat_filename_sed = args.output_dir / 'catalogs' / Path(merged_filename.stem + '_for_bagpipe.csv')
    
    zcosmos_idcol, cweb_idcol, c2020_idcol, c2025_idcol = 'ID_zcosmos', 'ID_cweb', 'ID_c2020', 'ID_c2025'

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
        if args.candels_cut:
            candels_regions = Regions.read(candels_regfile, format='ds9')
            print(f'Cross-matching zCOSMOS catalog and WFC3+ACS region...')
            df = df[candels_regions[1].contains(SkyCoord(df['ra'], df['dec'], unit='deg'), wcs_header)]
            print(f'...found {len(df)} objects in this region')
        print(f'Now applying z-cut off..')
        df_zsub = df[df['REDSHIFT'].between(1., 1.7)]
        df_zsub.rename(columns={'REDSHIFT':'redshift_zcosmos', 'id':zcosmos_idcol}, inplace=True)
        print(f'...found {len(df_zsub)} after the z-cut\n')
        '''
        # --------load COSMOS-Web catalog and crossmatch---------
        df_cweb = read_COSMOSWebb_catalog(args=args)
        df_cweb.rename(columns={'id':cweb_idcol}, inplace=True)

        ergs_s_cm2_hz_to_ujy_factor = 1e29 # 1 ergs/s/cm^2/Hz = 10^29 uJy
        web_fluxcols = [item for item in df_cweb.columns if 'FLUX_APER' in item]
        web_errcols = [item for item in df_cweb.columns if 'FLUX_ERR_APER' in item]
        for thiscol in web_fluxcols:
            df_cweb[thiscol] = df_cweb[thiscol] * ergs_s_cm2_hz_to_ujy_factor # converting from ergs/s/cm^2/Hz to micro Jansky
            df_cweb[thiscol.replace('FLUX', 'FLUX_ERR')] = df_cweb[thiscol.replace('FLUX', 'FLUX_ERR')] * ergs_s_cm2_hz_to_ujy_factor  # converting from ergs/s/cm^2/Hz to micro Jansky

        print(f'Cross-matching zCOSMOS-sub and COSMOS-Web catalogs; can take a while...')
        df_candels_cweb = get_crossmatch(df_zsub, df_cweb, sep_threshold=0.1, df1_idcol=zcosmos_idcol, df2_idcol=cweb_idcol)
        print(f'...found {len(df_candels_cweb)} cross-matched objects\n')

        df_candels_cweb.rename(columns={'df1_id':zcosmos_idcol, 'df2_id':cweb_idcol}, inplace=True)
        df_candels_cweb = pd.merge(df_candels_cweb, df_zsub, on=zcosmos_idcol, how='inner')
        
        df_cweb.rename(columns={'ra':'ra_cweb', 'dec':'dec_cweb'}, inplace=True)
        df_candels_cweb = pd.merge(df_candels_cweb, df_cweb, on=cweb_idcol, how='inner')

        # ------load COSMOS 2020 catalog and crossmatch---------
        df_c2020 = read_COSMOS2020_catalog(args=args)
        df_c2020.rename(columns={'id':c2020_idcol}, inplace=True)
        print(f'Cross-matching zOOSMOS-sub+COSMOS-Web and COSMOS-2020 catalogs; can take a while...')
        df_merged = get_crossmatch(df_candels_cweb, df_c2020, sep_threshold=0.1, df1_idcol=cweb_idcol, df2_idcol=c2020_idcol)
        print(f'...found {len(df_merged)} cross-matched objects\n')

        df_merged.rename(columns={'df1_id':cweb_idcol, 'df2_id':c2020_idcol}, inplace=True)
        df_merged = pd.merge(df_merged, df_candels_cweb, on=cweb_idcol, how='inner')

        df_c2020.rename(columns={'ra':'ra_c2020', 'dec':'dec_c2020'}, inplace=True)
        df_merged = pd.merge(df_merged, df_c2020, on=c2020_idcol, how='inner')
        '''
        # --------load COSMOS2025 catalog and crossmatch---------
        df_c2025 = read_COSMOS2025_catalog(args=args)
        df_c2025.rename(columns={'id':c2025_idcol}, inplace=True)

        c2025_fluxcols = [item for item in df_c2025.columns if 'flux_model_' in item]
        c2025_errcols = [item.replace('flux', 'flux_err-cal') for item in c2025_fluxcols]
        cols_to_extract = np.hstack([[c2025_idcol, 'ra', 'dec', 'ra_model', 'dec_model', 'sersic', 'sersic_err', 'tile', 'warn_flag', 'a_image', 'b_image', 'mag_auto_f444w'], c2025_fluxcols, c2025_errcols])
        df_c2025 = df_c2025[cols_to_extract]

        # loading LePhare and CIGALE COSMOS 2025 catalogs too
        df_c2025_lephare = Table.read(args.input_dir / 'COSMOS/COSMOSWeb_mastercatalog_v1_lephare.fits').to_pandas()
        df_c2025_cigale = Table.read(args.input_dir / 'COSMOS/COSMOSWeb_mastercatalog_v1_cigale.fits').to_pandas()
        df_c2025 = pd.concat([df_c2025, df_c2025_lephare[['ssfr_minchi2', 'mass_minchi2', 'sfr_minchi2', 'mass_l68', 'mass_med', 'mass_u68', 'sfr_l68', 'sfr_med', 'sfr_u68']], df_c2025_cigale[['mass', 'mass_err', 'sfr_100myr', 'sfr_100myr_err']]], axis=1)

        print(f'Cross-matching zCOSMOS-sub and COSMOS2025 catalogs; can take a while...')
        df_merged = get_crossmatch(df_zsub, df_c2025, sep_threshold=0.1, df1_idcol=zcosmos_idcol, df2_idcol=c2025_idcol)
        print(f'...found {len(df_merged)} cross-matched objects\n')

        df_merged.rename(columns={'df1_id':zcosmos_idcol, 'df2_id':c2025_idcol}, inplace=True)
        df_merged = pd.merge(df_merged, df_zsub, on=zcosmos_idcol, how='inner')
        
        df_merged.rename(columns={'ra':'ra_zcosmos', 'dec':'dec_zcosmos'}, inplace=True)
        df_merged = pd.merge(df_merged, df_c2025, on=c2025_idcol, how='inner')
        df_merged = df_merged[df_merged['warn_flag'] == 0]

        # -----saving the final catalog--------
        Table.from_pandas(df_merged).write(merged_filename, overwrite=True)
        print(f'Saved merged, z-filtered catalog as {merged_filename}')

    # ------plotting the targets on sky--------
    df_targets = Table.read(merged_filename).to_pandas()
    df_targets['log_mass'] = np.log10(df_targets['mass'])
    ax = fig.gca()
    color = 'g'
    for index, row in df_targets.iterrows():
        circle = plt.Circle(xy=(row['ra'], row['dec']), radius=50/3600, color=color, alpha=1, fill=True, transform=ax.get_transform('fk5'))
        #ax.text(row['ra'] + 0.01, row['dec'] + 0.01, row[cweb_idcol], color=color, fontsize=10, transform=ax.get_transform('fk5'), va='bottom', ha='left', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
        ax.add_patch(circle)

    # --------saving the plot----------
    figname = args.input_dir / 'footprints' / f'{args.bg}_with_footprints_with_{merged_filename.stem}.png'
    fig.savefig(figname, transparent=args.fortalk)
    print(f'Saved plot to {figname}')
    plt.show(block=False)
    
    # --------modifying columns for SED fitting-------------
    df_sed = df_targets.copy()
    df_sed.rename(columns={'redshift_zcosmos':'redshift'}, inplace=True)

    cols_to_retain = ['flux', 'id', 'redshift', 'ra', 'dec']
    cols_to_drop1 = ['ID', 'ID_c2020', 'ID_zcosmos', 'ID_COSMOS2015', 'objid', 'ra_cweb', 'dec_cweb', 'ra_c2020', 'dec_c2020', 'ra_c2025', 'dec_c2025', 'ra_zcosmos', 'dec_zcosmos'] + [item for item in df_sed.columns if not np.array([item2 in item.lower() for item2 in cols_to_retain]).any()]
    df_sed = df_sed.drop(columns=cols_to_drop1, errors='ignore')

    c2025_filter_dict = {'hsc-nb1010': 'HSC_nb1010', 'hsc-z': 'HSC_z', 'uvista-nb118': 'UVISTA_NB118', 'hsc-y': 'HSC_Y', 'sc-ib709': 'SC_IB709', 'sc-ia767': 'SC_IB767', 'f115w': 'F115W', 'uvista-y': 'UVISTA_Y', 'sc-ia738': 'SC_IB738', 'irac-ch3': 'IRAC_CH3', 'hsc-g': 'HSC_g', 'hsc-r': 'HSC_r', 'f770w': 'F770W', 'sc-nb711': 'SC_NB711', 'irac-ch2': 'IRAC_CH2', 'sc-nb816': 'SC_NB816', 'sc-ia624': 'SC_IB624', 'hsc-nb0921': 'HSC_nb0921', 'sc-ia527': 'SC_IB527', 'sc-ia484': 'SC_IB484', 'sc-ib827': 'SC_IB827', 'f444w': 'F444W', 'sc-ia679': 'SC_IB679', 'uvista-j': 'UVISTA_J', 'f150w': 'F150W', 'f277w': 'F277W', 'hsc-nb0816': 'HSC_nb0816', 'hsc-i': 'HSC_i', 'hst-f814w': 'ACS_F814W', 'uvista-ks': 'UVISTA_Ks', 'irac-ch4': 'IRAC_CH4', 'sc-ib427': 'SC_IB427', 'sc-ib505': 'SC_IB505', 'cfht-u': 'CFHT_U', 'sc-ib574': 'SC_IB574', 'irac-ch1': 'IRAC_CH1', 'uvista-h': 'UVISTA_H'}

    for thiscol in df_sed.columns:
        if 'FLUX_APER_' in thiscol: df_sed.rename(columns={thiscol: thiscol[10:] + '_sci'}, inplace=True)
        elif 'FLUX_ERR_APER_' in thiscol: df_sed.rename(columns={thiscol: thiscol[14:] + '_err'}, inplace=True)
        elif '_FLUX_AUTO' in thiscol: df_sed.rename(columns={thiscol: thiscol[:-10] + '_sci'}, inplace=True)
        elif '_FLUXERR_AUTO' in thiscol: df_sed.rename(columns={thiscol: thiscol[:-13] + '_err'}, inplace=True)
        elif '_FLUXERR' in thiscol: df_sed.rename(columns={thiscol: thiscol[:-8] + '_err'}, inplace=True)
        elif '_FLUX' in thiscol: df_sed.rename(columns={thiscol: thiscol[:-5] + '_sci'}, inplace=True)
        elif 'flux_model_' in thiscol: df_sed.rename(columns={thiscol: c2025_filter_dict[thiscol[11:]] + '_sci'}, inplace=True)
        elif 'flux_err-cal_model_' in thiscol: df_sed.rename(columns={thiscol: c2025_filter_dict[thiscol[19:]] + '_err'}, inplace=True)


    for thiscol in df_sed.columns:
        if 'F115W' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('F115W', 'NIRCAM_F115W')}, inplace=True)
        elif 'F150W' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('F150W', 'NIRCAM_F150W')}, inplace=True)
        elif 'F277W' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('F277W', 'NIRCAM_F277W')}, inplace=True)
        elif 'F444W' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('F444W', 'NIRCAM_F444W')}, inplace=True)
        elif 'F770W' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('F770W', 'MIRI_F770W')}, inplace=True)
        elif 'IA' in thiscol: df_sed.rename(columns={thiscol: thiscol.replace('IA', 'IB')}, inplace=True)
    
    filters_to_omit = ['SPLASH', 'HSC_nb']
    cols_to_omit = np.hstack([df_sed.columns[df_sed.columns.str.contains(filter_to_omit)].tolist() for filter_to_omit in filters_to_omit]) # to remove SPLASH filters because those transmission files are absent
    cols_with_nan = df_sed.columns[df_sed.isna().any()] # to remove filters with missing fluxes
    cols_with_negative = df_sed.columns[(df_sed < 0).any()] # to remove filters with negative fluxes
    cols_to_drop2 = np.hstack([cols_to_omit, cols_with_nan, cols_with_negative])
    df_sed = df_sed.drop(columns=cols_to_drop2)

    print(f'\nRemoved {len(cols_to_omit)} columns ({cols_to_omit}) with {filters_to_omit}\n\
          {len(cols_with_nan)} columns ({cols_with_nan}) with NaNs\n\
          and {len(cols_with_negative)} columns ({cols_with_negative}) with negative\n\
          total {len(cols_to_drop1) + len(cols_to_drop2)} columns removed, leaving {len(df_sed.columns)} columns\n')

    # -----saving the final catalog--------
    df_sed.to_csv(photcat_filename_sed, index=None)
    print(f'Saved just the flux catalog, for bagpipes, as {photcat_filename_sed}')
    
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
