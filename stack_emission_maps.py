'''
    Filename: stack_emission_maps.py
    Notes: Stacks 2D emission line maps (from grizli) for objects in a given field/s, in bins of mass and/or SFR
    Author : Ayan
    Created: 15-01-26
    Example: run stack_emission_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run stack_emission_maps.py --field Par28 --debug_align
'''

from header import *
from util import *
from make_diagnostic_maps import get_re, get_cutout, get_offsets_from_center, trim_image, get_dereddened_flux, myimshow

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_passage_masses_from_cosmos(df, id_col='objid', field_col='field', cosmos_idcol='id'):
    '''
    Derives stellar masses of PASSAGE galaxies present in the given dataframe from COSMOS-Web catalog
    Returns dataframe
    '''
    passage_fields = [item for item in np.unique(df[field_col]) if 'Par' in item]
    
    df_cosmos = pd.DataFrame()
    for index, thisfield in enumerate(passage_fields):
        cosmosfilename = args.input_dir / 'COSMOS' /  f'cosmoswebb_objects_in_{thisfield}.fits'
        df_cosmos_thisfield = Table.read(cosmosfilename).to_pandas()
        sed_cols_to_extract = ['passage_id', cosmos_idcol, 'ra', 'dec', 'mass_med', 'sfr_med', 'ssfr_med']
        df_cosmos = pd.concat([df_cosmos, df_cosmos_thisfield[sed_cols_to_extract]])
    
    df_cosmos = df_cosmos.rename(columns={cosmos_idcol: 'cosmos_id'})
    df_cosmos['passage_id'] = df_cosmos['passage_id'].astype(str)

    df['passage_id'] = df[field_col].astype(str) + '-' + df[id_col].astype(str)  # making a unique combination of field and object id
    df = pd.merge(df, df_cosmos, on=['passage_id'], how='right')
    df = df.drop('passage_id', axis=1)
    df = df.rename(columns={'mass_med': 'log_mass', 'sfr_med': 'log_sfr', 'ssfr_med': 'log_ssfr'})
    df = df[(df['log_mass'] > 0) & (df['log_sfr'] > 0)]

    return df

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, full_hdu, args, dered=True, silent=True):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line map and the corresponding 2D uncertainty map
    '''
    # ---------getting the spatially resolved line flux map----------------
    try:
        line_hdu = full_hdu['LINE', line]
        line_wht = full_hdu['LINEWHT', line]
    except KeyError:
        if line == 'OIII':
            line_hdu = full_hdu['LINE', 'OIII-5007']
            line_wht = full_hdu['LINEWHT', 'OIII-5007']

    line_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2
    line_map_err = 1e-17 / (line_wht.data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
    line_wave = line_hdu.header['RESTWAVE'] # in Angstrom

    # ----------deblending flux--------------------
    factor = 1.0
    if not args.do_not_correct_flux:
        if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and remove the OIII 4959 component
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            if not silent: print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'Ha': # special treatment for Ha line, in order to account for and remove the NII component
            factor = 0.823 # from James et al. 2023?
            if not silent: print(f'Correcting Ha for NII component, by factor of {factor:.3f}')

    line_map = line_map * factor
    line_map_err = line_map_err * factor

    # -------------pixel offset to true center----------------
    line_map = np.roll(line_map, args.ndelta_xpix, axis=0)
    line_map_err = np.roll(line_map_err, args.ndelta_xpix, axis=0)
    line_map = np.roll(line_map, args.ndelta_ypix, axis=1)
    line_map_err = np.roll(line_map_err, args.ndelta_ypix, axis=1)

    # ----------getting a smaller cutout around the object center-----------
    line_map = trim_image(line_map, args)
    line_map_err = trim_image(line_map_err, args)

    # -------masking outside segmentation map---------
    line_map = np.ma.masked_where(args.segmentation_map != args.id, line_map)
    line_map_err = np.ma.masked_where(args.segmentation_map != args.id, line_map_err)

    # -----------getting the dereddened flux value-----------------
    if dered and args.EB_V != 0:
        line_map_quant = get_dereddened_flux(unp.uarray(line_map, line_map_err), line_wave, args.EB_V)
        line_map = unp.nominal_values(line_map_quant)
        line_map_err = unp.std_devs(line_map_quant)

    # ---------converting from flux to surface brightness units------------
    pixscale_kpc = args.pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(args.z).value # kpc
    line_map /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)
    line_map_err /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)

    return line_map, line_map_err

# --------------------------------------------------------------------------------------------------------------------
def rotate_line_map(line_map, args):
    '''
    Accepts a 2D map, rotates it to align the major axis with horizontal, based on PA values (contained within args.pa)
    Returns rotated 2D map
    '''
    # reshape=False to keep the galaxy centered in the same array size
    rotated_map = ndimage.rotate(line_map, -args.pa, reshape=False, order=1, cval=np.nan) # cval=np.nan ensures pixels rotated "in" from the outside don't mess up the median stack
    
    return rotated_map

# --------------------------------------------------------------------------------------------------------------------
def deproject_line_map(line_map, args):
    '''
    Accepts a 2D map, deprojects it assuming galaxy major axis is horizontal (axis 1), based on a and b values (contained within args.semi_major and args.semi_minor)
    Returns deprojected 2D map
    '''
    stretch_factor = args.semi_major / args.semi_minor
    deprojected_map = ndimage.zoom(line_map, (stretch_factor, 1.0), order=1) # zoom axis 0 (y/minor) by stretch_factor and axis 1 (x/major) by 1.0
    
    return deprojected_map

# --------------------------------------------------------------------------------------------------------------------
def rescale_line_map(line_map, args, re_extent=2, npix_side=20):
    '''
    Accepts a 2D map, rescales it to a standardised grid of +/- re_extent effective radii  (contained within args.re_arcsec and args.re_kpc) from the center plit over npix pixels a side
    Returns rescaled 2D map
    '''
    target_re_px = 2 * re_extent / npix_side
    current_re_px = args.re_arcsec / args.pix_size_arcsec
    
    zoom_factor = target_re_px / current_re_px
    rescaled_map = ndimage.zoom(line_map, zoom_factor, order=1)
    
    center_y, center_x = np.array(rescaled_map.shape) // 2
    rescaled_map = rescaled_map[center_y - int(npix_side/2) : center_y + int(npix_side/2), 
                         center_x - int(npix_side/2) : center_x + int(npix_side/2)]
    
    # if rescaled_map.shape != (npix_side, npix_side):
    #     canvas = np.full((npix_side, npix_side), np.nan)
    #     h, w = rescaled_map.shape
    #     canvas[:h, :w] = rescaled_map # This handles cases where the galaxy is at the edge of the original frame
    #     return canvas

    return rescaled_map

# --------------------------------------------------------------------------------------------------------------------
def stack_line_maps(line_maps_array, line_maps_err_array):
    '''
    Accepts 2 x N array (and its corresponding uncertainties) of N emission line maps from N objects and stacks them, appropriately accounting for the uncertainties
    Returns stacked 2D line map and 2D uncertainty map
    '''

    stacked_line_map = np.median(line_maps_array, axis=2)
    stacked_line_map_err = np.median(line_maps_err_array, axis=2)
    
    return stacked_line_map, stacked_line_map_err

# --------------------------------------------------------------------------------------------------------------------
def quickplot_line_map(objid, field='Par028', line='OIII', arcsec_limit=1, no_plot=False):
    '''
    For quick and dirty displaying of a specific line map extension from the maps/full file of a specific object
    '''
    filename = args.input_dir / field / 'Products'/ 'maps' / f'{field}_{objid:05d}.maps.fits'
    full_hdu = fits.open(filename)
    pix_size_arcsec = full_hdu[5].header['PIXASEC']
    available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
    segmentation_map = full_hdu['SEG'].data

    if line in available_lines:
        line_hdu = full_hdu['LINE', line]
        line_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2
        line_map = trim_image(line_map, arcsec_limit=arcsec_limit, pix_size_arcsec=pix_size_arcsec)
        segmentation_map = trim_image(segmentation_map, arcsec_limit=arcsec_limit, pix_size_arcsec=pix_size_arcsec)

        if no_plot:
            print(f'{objid} should work')
            return objid
        else:
            fig, ax = plt.subplots(1)
            ax = myimshow(line_map, ax, label=f'{field} {objid}: {line}', contour=segmentation_map != objid, cmap='viridis', fontsize=15)
            
            plt.show(block=False)
            return fig
    else:
        if no_plot: print(f'{line} is not available among {available_lines}')
        else: print(f'{objid} will not work')

        return -99

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontsize = 15
    args.line_list = ['OIII', 'OII', 'NeIII-3867', 'Hb', 'OIII-4363', 'Ha', 'SII']

    # ----------declaring mass and SFR bins-------------------
    delta_log_mass, delta_log_sfr = 1, 0.5
    log_mass_bins = np.arange(7.5, 11.5 + delta_log_mass/2, delta_log_mass)
    log_sfr_bins = np.arange(-2, 2 + delta_log_sfr/2, delta_log_sfr)
    
    # ---------determining filename suffixes-------------------------------
    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / 'Par*') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr

    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'\nStarting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / args.field
        fig_dir = output_dir / 'stacked_plots'
        fig_dir.mkdir(parents=True, exist_ok=True)
        if args.re_extract: output_dir = output_dir / 're_extracted'
        output_dir.mkdir(parents=True, exist_ok=True)

        # ---------read the photometric catalog file--------------------
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        df['field'] = args.field

        # ---------crossmatch with cosmos-web to get stellar mass and SFR--------------------
        df = get_passage_masses_from_cosmos(df, id_col='id')

        # ---------merge with effective radius catalog--------------------
        df_re = Table.read(args.output_dir / f'catalogs/{args.field}_re_list.fits').to_pandas()
        if 'redshift' in df_re: df_re.drop(columns=['redshift'], axis=1, inplace=True)
        df_re['id'] = df_re['id'].astype(int)
        df = pd.merge(df, df_re, on='id', how='inner')

        # ----------binning the dataframe by mass and SFR bins-------------------
        df['mass_interval'] = pd.cut(df['log_mass'], bins=log_mass_bins)
        df['sfr_interval'] = pd.cut(df['log_sfr'], bins=log_sfr_bins)
        df = df.dropna(subset=['mass_interval', 'sfr_interval'])
        df['bin_intervals'] = list(zip(df['mass_interval'], df['sfr_interval']))

        all_mass_intervals = df['mass_interval'].cat.categories
        all_sfr_intervals = df['sfr_interval'].cat.categories
        bin_list = list(itertools.product(all_mass_intervals, all_sfr_intervals))
                
        # ------------looping over each bin-----------------------
        for index2, this_mass_sfr_bin in enumerate(bin_list):
            start_time3 = datetime.now()
            print(f'\tStarting bin ({index2 + 1}/{len(bin_list)}) mass={this_mass_sfr_bin[0]}, sfr={this_mass_sfr_bin[1]}..', end=' ')
            # --------determine which objects fall in this bin----------
            mask = df['bin_intervals'] == this_mass_sfr_bin
            id_arr = df.loc[mask, 'id'].tolist()
            id_arr = np.array(id_arr)[np.isin(id_arr, [2698, 2583, 2171, 672])] # subset of IDs of "good-looking" galaxies, just for testing, otherwise comment this line out

            if len(id_arr) > 0:
                print(f'which has {len(id_arr)} objects.')
                line_maps_array, line_maps_err_array = [], [] # each array must be len(args.line_list) long and each element will become a stack of 2D images
                
                # ------------looping over the objects-----------------------
                for index3, this_id in enumerate(id_arr):
                    args.id = this_id
                    print(f'\t\tCommencing ({index3 + 1}/{len(id_arr)}) ID {args.id}..')
                    obj = df.loc[df['id']==this_id].iloc[0]

                    # ------determining directories and filenames---------
                    full_fits_file = product_dir / 'full' / f'{args.field}_{args.id:05d}.full.fits'
                    maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'
                    od_filename = product_dir / 'spec1D' / f'{args.field}_{args.id:05d}.1D.fits'

                    if os.path.exists(maps_fits_file): # if the fits files are in maps/
                        full_filename = maps_fits_file
                    elif os.path.exists(full_fits_file): # if the fits files are in full/
                        full_filename = full_fits_file
                    else:
                        print(f'Could not find {full_fits_file} or {maps_fits_file}, so skipping it.')
                        continue
                    if not os.path.exists(od_filename):
                        od_filename = Path(str(od_filename).replace('.1D.', '.spec1D.'))

                    # ------------read in maps files--------------------------------
                    full_hdu = fits.open(full_filename)
                    od_hdu = fits.open(od_filename)

                    # ----------determining object parameters------------
                    args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
                    args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
                    args.z = full_hdu[0].header['REDSHIFT']
                    args.nlines = full_hdu[0].header['NUMLINES']
                    args.distance = cosmo.luminosity_distance(args.z)
                    args.pix_arcsec = full_hdu[5].header['PIXASEC']
                    args.mag = obj['mag_auto']

                    line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
                    args.pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
                    imsize_arcsec = full_hdu['DSCI'].data.shape[0] * args.pix_size_arcsec
                    offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
                    args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)
                    
                    args.EB_V = 0. # until it gets over-written, if both H alpha and H beta lines are present
                    args.semi_major = obj['a_image']
                    args.semi_minor = obj['b_image']
                    args.pa = obj['theta_image']
                    args.re_kpc, args.re_arcsec = obj['re_kpc'], obj['re_arcsec']

                    # --------determining true center of object---------------------
                    args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W')
                    
                    # ---------------segmentation map---------------
                    segmentation_map = full_hdu['SEG'].data
                    segmentation_map = np.roll(segmentation_map, args.ndelta_xpix, axis=0)
                    segmentation_map = np.roll(segmentation_map, args.ndelta_ypix, axis=1)
                    args.segmentation_map = trim_image(segmentation_map, args)
                    deprojected_segmentation_map = deproject_line_map(args.segmentation_map, args)
                    rotated_segmentation_map = rotate_line_map(deprojected_segmentation_map, args)
                    rescaled_segmentation_map = rescale_line_map(rotated_segmentation_map, args)

                    # ---------looping over all emission lines------------------------------
                    for index4, this_line in enumerate(args.line_list):
                        if this_line not in args.available_lines:
                                print(f'\t\t\tLine {this_line} not available for object {args.id}. So skipping this line..')
                                continue
                        # -------determining true center of object-------------
                        args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W')

                        # -----------extracting the emission line map----------------
                        line_map, line_map_err = get_emission_line_map(this_line, full_hdu, args, dered=True, silent=False)

                        # --------rotating the line map---------------------
                        rotated_line_map = rotate_line_map(line_map, args)
                        rotated_line_map_err = rotate_line_map(line_map_err, args)

                        # --------deprojecting the line map---------------------
                        deprojected_line_map = deproject_line_map(rotated_line_map, args)
                        deprojected_line_map_err = deproject_line_map(rotated_line_map_err, args)

                        # --------scaling the line map by effective radius---------------------
                        rescaled_line_map = rescale_line_map(deprojected_line_map, args)
                        rescaled_line_map_err = rescale_line_map(deprojected_line_map_err, args)

                        # -------plotting the intermediate steps: for debugging--------------
                        if args.debug_align:
                            fig, axes = plt.subplots(1, 4, figsize=(12, 4))
                            fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)
                            cmap = 'viridis'
                            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: {this_line} map: Alignment diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')
                            
                            re_pix = args.re_arcsec / args.pix_size_arcsec
                            axes[0] = myimshow(line_map, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=cmap)
                            axes[1] = myimshow(rotated_line_map, axes[1], contour=rotated_segmentation_map != args.id, re_pix=re_pix, label='Rotated', cmap=cmap)
                            axes[2] = myimshow(deprojected_line_map, axes[2], contour=deprojected_segmentation_map != args.id, re_pix=re_pix, label='Deprojected', cmap=cmap)
                            axes[3] = myimshow(rescaled_line_map, axes[3], contour=rescaled_segmentation_map != args.id, re_pix=1, label='Rescaled', cmap=cmap)

                            plt.show(block=False)
                            sys.exit(f'\t\t\tExiting here because of --debug_re mode; if you want to run the full code as usual then remove the --debug_re option and re-run')

                        # ---------------appending the line map-------------------------
                        line_maps_array[index4].append(rescaled_line_map)
                        line_maps_err_array[index4].append(rescaled_line_map_err)
                
                # -----stacking all line maps for all objects in this bin-----------
                stacked_maps, stacked_maps_err = [], [] # each array must be len(args.line_list) long and each element will be a 2D image
                for index4, this_line in enumerate(args.line_list):
                    stacked_maps[index4], stacked_maps_err[index4] = stack_line_maps(line_maps_array[index4], line_maps_err_array[4])

                # -------writing out stacked line maps as fits files--------------

                # -----------------plot emission line maps of this bin---------------
                
                # -----------------plot emission line ratio maps of this bin---------------
                
            else:
                print(f'which has no object. Skipping this bin.')
                continue
            
            print(f'Completed bin {this_mass_sfr_bin} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')

        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
