'''
    Filename: stack_emission_maps.py
    Notes: Stacks 2D emission line maps (from grizli) for objects in a given field/s, in bins of mass and/or SFR
    Author : Ayan
    Created: 15-01-26
    Example: run stack_emission_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28 --do_all_obj
             run stack_emission_maps.py --field Par28 --debug_align --id 2822,2698,2583,2171,672
             run stack_emission_maps.py --field Par28 --do_all_obj --re_limit 2
             run stack_emission_maps.py --field Par28 --do_all_obj --adaptive_bins --max_gal_per_bin 20
'''

from header import *
from util import *
from make_diagnostic_maps import trim_image, get_dereddened_flux, myimshow, get_offsets_from_center, get_cutout

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_direct_image(full_hdu, filter, args):
    '''
    Loads the direct image for a given filter for a given object
    Returns the image
    '''
    try:
        hdu = full_hdu['DSCI', filter.upper()]
        image = hdu.data
        exptime = 1
    except:
        try:
            hdu = full_hdu['DSCI', f'{filter}-CLEAR']
            image = hdu.data
            exptime = full_hdu[0].header[f'T_{filter.upper()}']
        except:
            try:
                hdu = full_hdu['DSCI', f'{filter.upper()}-{filter.upper()}-CLEAR']
                image = hdu.data
                exptime = full_hdu[0].header[f'T_{filter.upper()}']
            except:
                full_field_filename = args.root_dir / f'passage_data/' / f'{args.drv}' / f'{args.field}' / 'Products' / f'{args.field}_{filter.lower()}-clear_drz_sci.fits'
                print(f'{filter.upper()} not found in full_hdu extension. Therefore trying to get cutout from full field image {full_field_filename}')
            
                exptime = fits.open(full_field_filename)[0].header['EXPTIME']
                pos = SkyCoord(full_hdu[0].header['RA'], full_hdu[0].header['DEC'], unit = 'deg')
                size = 2 * args.arcsec_limit * u.arcsec
                target_header = full_hdu['DSCI', 'F140W'].header
                
                temp1, temp2 = args.only_seg, args.vorbin
                args.only_seg, args.vorbin = False, False
                image = get_cutout(full_field_filename, pos, size, target_header, args, plot_test_axes=None, skip_re_trim=True)
                args.only_seg, args.vorbin = temp1, temp2

    image = ndimage.shift(image, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
    image = trim_image(image, args, skip_re_trim=True)
 
    return image, exptime

# --------------------------------------------------------------------------------------------------------------------
def get_passage_masses_from_cosmos(df, args, id_col='objid', field_col='field', cosmos_idcol='id'):
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
    df = df.dropna(subset=['log_mass', 'log_sfr'], axis=0)

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

    # ----------getting a smaller cutout around the object center-----------
    line_map = trim_image(line_map, args, skip_re_trim=True)
    line_map_err = trim_image(line_map_err, args, skip_re_trim=True)

    # ----------re-center line map-----------
    line_map = ndimage.shift(line_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
    line_map_err = ndimage.shift(line_map_err, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)

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
def rescale_line_map(line_map, args):
    '''
    Accepts a 2D map, rescales it to a standardised grid of +/- re_extent effective radii  (contained within args.re_arcsec and args.re_kpc) from the center plit over npix pixels a side
    Returns rescaled 2D map
    '''
    npix_side = args.npix_side
    stretch_factor = args.semi_major / args.semi_minor
    target_re_px = npix_side / (2 * args.re_limit)
    current_re_px = args.re_arcsec / args.pix_size_arcsec
    zoom_factor = target_re_px / current_re_px
    #if args.debug_align: print(f'Deb133: id {args.id}: re_arcsec={args.re_arcsec}, pix_size_arcsec={args.pix_size_arcsec}, semi_major={args.semi_major}, semi_minor={args.semi_minor}, stretch_factor={stretch_factor}, target_re_px={target_re_px}, current_re_px={current_re_px}, zoom_factor={zoom_factor}') ##

    rescaled_map = ndimage.zoom(line_map, (zoom_factor * stretch_factor, zoom_factor), order=1)
    
    center_y, center_x = np.array(rescaled_map.shape) // 2
    rescaled_map = rescaled_map[center_y - npix_side // 2 : center_y + npix_side // 2, 
                         center_x - npix_side // 2 : center_x + npix_side // 2]
    
    if rescaled_map.shape != (npix_side, npix_side):
        canvas = np.full((npix_side, npix_side), np.nan)
        h, w = rescaled_map.shape
        start_y = (npix_side - h) // 2
        start_x = (npix_side - w) // 2
        canvas[start_y : start_y + h, start_x : start_x + w] = rescaled_map # This handles cases where the galaxy is at the edge of the original frame
        return canvas

    return rescaled_map

# --------------------------------------------------------------------------------------------------------------------
def get_center_offsets(dir_img, args, silent=False):
    '''
    Accepts a 2D map, then computes the offset from the original center of image to the brightest pixel
    Returns two integers (offset in x and y axes)
    '''
    smoothing_kernel = Box2DKernel(5, mode=args.kernel_mode)
    dir_img_smoothed = convolve(dir_img, smoothing_kernel)
    smooth_shape = dir_img_smoothed.shape
    ncells = 5 # only searches within +/- 5 cells of the original center
    dir_img_smoothed_subarea = dir_img_smoothed[smooth_shape[0] // 2 - ncells: smooth_shape[0] // 2 + ncells, smooth_shape[1] // 2 - ncells: smooth_shape[1] // 2 + ncells] # only searching for brightest pixel in the vicinity of the original center, otherwise might pick up neighbouring galaxies
    brightest_coords = np.where(dir_img_smoothed == dir_img_smoothed_subarea.max())
    brightest_x, brightest_y = brightest_coords[0][0], brightest_coords[1][0]
    cen_x, cen_y = int(np.shape(dir_img)[0] / 2), int(np.shape(dir_img)[1] / 2)
    ndelta_xpix = cen_x - brightest_x
    ndelta_ypix = cen_y - brightest_y
    if not silent: print(f'For {args.field}:{args.id}: Determined x and y offsets from direct image = {ndelta_xpix}, {ndelta_ypix}')

    return ndelta_xpix, ndelta_ypix

# --------------------------------------------------------------------------------------------------------------------
def weighted_stack_line_maps(line_maps_array, line_maps_err_array):
    '''
    Accepts 2 x N array (and its corresponding uncertainties) of N emission line maps from N objects and sums them, weighting by the uncertainties
    Returns stacked 2D line map and 2D uncertainty map
    '''
    line_maps_array = np.array(line_maps_array)
    line_maps_err_array = np.array(line_maps_err_array)

    weights_array = 1.0 / (line_maps_err_array ** 2)
    stacked_line_map_err = np.sqrt(1.0 / np.nansum(weights_array, axis=0))    

    weights_array /= 1e37 # dividing out this typical large factor, otherwise upon summing across many objects it becomes inf
    stacked_line_map = np.nansum(line_maps_array * weights_array, axis=0) / np.nansum(weights_array, axis=0)
    
    return stacked_line_map, stacked_line_map_err

# --------------------------------------------------------------------------------------------------------------------
def quickplot_raw_line_map(objid, field='Par028', line='OIII', arcsec_limit=1, no_plot=False, cmin=None, cmax=None, takelog=False, col='w'):
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
            if takelog: line_map = np.log10(line_map)
            ax = myimshow(line_map, ax, label=f'{field} {objid}: {line}', contour=segmentation_map != objid, cmap='viridis', fontsize=15, cmin=cmin, cmax=cmax, col=col)
            
            plt.show(block=False)
            return fig
    else:
        if no_plot: print(f'{line} is not available among {available_lines}')
        else: print(f'{objid} will not work as {line} is not available for this object')

        return -99

# --------------------------------------------------------------------------------------------------------------------
def create_re_wcs(args):
    '''
    Creates a fake WCS for storing the stacked emission line maps in FITS format
    Returns WCS header
    Borrowed from Gemini
    '''
    cen_pix = args.npix_side // 2 + 0.5
    re_per_pix = (2 * args.re_limit) / args.npix_side

    w = pywcs.WCS(naxis=2)
    w.wcs.crpix = [cen_pix, cen_pix] 
    w.wcs.cdelt = [re_per_pix, re_per_pix]
    w.wcs.crval = [0, 0]
    w.wcs.ctype = ["RE-X", "RE-Y"]
    return w.to_header()

# --------------------------------------------------------------------------------------------------------------------
def write_stacked_maps(stacked_maps, stacked_maps_err, ids_arr, output_filename, args):
    '''
    Accept a list of stacked emission line maps (and corresponding uncertainties) and writes them out into one multi-extension fits file
    '''
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['CONTENT'] = 'Stacked Galaxy Emission Line Maps'
    primary_hdu.header['NPIX'] = (args.npix_side, 'Pixels per side')
    pix_per_re = args.npix_side / (2 * args.re_limit)
    primary_hdu.header['RE_PX'] = (pix_per_re, f'Scale: 1 Re = {pix_per_re} pixels')
    
    hdu_list = [primary_hdu]

    for index, this_line in enumerate(args.line_list):
        if type(stacked_maps[index]) != np.ndarray:
            print(f'\t{this_line} extension is completely empty in the stack, hence not writing it in the file.')
            continue
        flux_hdu = fits.ImageHDU(stacked_maps[index], header=create_re_wcs(args), name=this_line.upper())
        err_hdu = fits.ImageHDU(stacked_maps_err[index], header=create_re_wcs(args), name=f'{this_line.upper()}_ERR')
        col = fits.Column(name='OBJID', format='20A', array=np.array(ids_arr[index]))
        table_hdu = fits.BinTableHDU.from_columns([col], name=f'{this_line.upper()}_IDS')
        flux_hdu.header['NOBJ'] = (len(ids_arr[index]), '#galaxies that went into this stack')
        hdu_list.extend([flux_hdu, err_hdu, table_hdu])

    hdul = fits.HDUList(hdu_list)
    hdul.writeto(output_filename, overwrite=True)
    print(f'Saved {len(args.line_list)*2} extensions to {output_filename}')

# --------------------------------------------------------------------------------------------------------------------
def read_stacked_maps(output_filename, args):
    '''
    Reads the stacked emission line maps from one multi-extension fits file
    Returns 3 lists: fluxes, uncertainties, and number of objects, with each element of each list corresponding to a particular emission line
    '''
    hdul = fits.open(output_filename)
    line_names = [hdu.name for hdu in hdul if hdu.name != 'PRIMARY' and '_ERR' not in hdu.name and not '_ID' in hdu.name]
    line_dict = {}
    for this_line in line_names:
        stacked_map = hdul[f'{this_line.upper()}'].data
        stacked_map_err = hdul[f'{this_line.upper()}_ERR'].data
        ids = Table(hdul[f'{this_line.upper()}_IDS'].data)['OBJID'].value
        nobj = hdul[f'{this_line.upper()}'].header['NOBJ']

        line_dict.update({this_line: np.ma.masked_where(False, unp.uarray(stacked_map, stacked_map_err)), f'{this_line}_ids': ids, f'{this_line}_nobj': nobj})
    
    return line_dict

# --------------------------------------------------------------------------------------------------------------------
def get_adaptive_bins(df_subset, m_range, s_range, max_n=20):
    '''
    m_range: (min, max) of log_mass for this specific tile
    s_range: (min, max) of log_sfr for this specific tile
    Courtesy of this function: Gemini
    '''
    # Count how many galaxies are in this specific rectangular area
    mask = (df_subset['log_mass'] >= m_range[0]) & (df_subset['log_mass'] < m_range[1]) & \
           (df_subset['log_sfr'] >= s_range[0]) & (df_subset['log_sfr'] < s_range[1])
    
    subset = df_subset[mask]
    n_count = len(subset)

    # Base case: if count is small OR area is already very tiny, stop splitting
    if n_count <= max_n or (m_range[1] - m_range[0]) < 0.1:
        if n_count == 0: return []
        # Return the coordinates and the mean value for this leaf node
        return [{'m_min': m_range[0], 'm_max': m_range[1], 's_min': s_range[0], 's_max': s_range[1], 'n_count': n_count}]
    
    # Recursive step: Split into 4 quadrants
    m_mid = (m_range[0] + m_range[1]) / 2
    s_mid = (s_range[0] + s_range[1]) / 2
    
    results = []
    results.extend(get_adaptive_bins(subset, (m_range[0], m_mid), (s_range[0], s_mid))) # Bottom-Left
    results.extend(get_adaptive_bins(subset, (m_mid, m_range[1]), (s_range[0], s_mid))) # Bottom-Right
    results.extend(get_adaptive_bins(subset, (m_range[0], m_mid), (s_mid, s_range[1]))) # Top-Left
    results.extend(get_adaptive_bins(subset, (m_mid, m_range[1]), (s_mid, s_range[1]))) # Top-Right
    
    return results

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True, segmentation_map=None, in_re_units=True, seg_col='k'):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''

    if takelog: image =  np.log10(image.data)

    if in_re_units:
        offset = args.pix_size_re / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
        args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)
    else:
        offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
        args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    if segmentation_map is not None:
        ax.contour(segmentation_map != args.id, extent=args.extent, levels=0, colors=seg_col, linewidths=0.5)
        limit = args.re_limit * args.re_arcsec
        ax.add_patch(plt.Rectangle((-limit, -limit), 2 * limit, 2 * limit, lw=0.5, color='r', fill=False))
    
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.set_aspect('auto') 
    
    ax = annotate_axes(ax, '', '', fontsize=args.fontsize / args.fontfactor, label=label, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def setup_fullpage_figure(n_page, n_total_pages, n_lines, cmin, cmax, cmap, args, in_re_units=True):
    '''
    Initialises a full page figure for plotting all the emission line maps of individual galaxies along with the stacked maps
    Returns figure and axes handle
    '''
    fig, axes = plt.subplots(args.max_gal_per_page, n_lines, figsize=(1.5 * n_lines, 1.5 * args.max_gal_per_page))#, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.08, wspace=0., hspace=0.)

    # ------------adding colorbars at top of the page---------------
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    cmap = plt.get_cmap(cmap)

    sm = mpl_cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, location='top', shrink=0.95, pad=0.01, aspect=60)
    cbar.set_label('Surface brightness [ergs/s/cm^2/kpc^2]', fontsize=args.fontsize, labelpad=5)    
    cbar.ax.tick_params(labelsize=args.fontsize / args.fontfactor)

    # ------------adding axis labels and page numbers for the whole page---------------
    label = r'Offset (R$_e$)' if in_re_units else r'Offset (arcseconds)'
    fig.supxlabel(label, fontsize=args.fontsize)
    fig.supylabel(label, fontsize=args.fontsize)

    page_text = f'Page {n_page} of {n_total_pages}'
    fig.text(0.95, 0.01, page_text, ha='right', fontsize=args.fontsize / args.fontfactor, color='k')

    return fig, axes

# ----------declaring mass and SFR bins-------------------
delta_log_mass, delta_log_sfr = 1, 0.5
log_mass_bins = np.arange(7.5, 11.5 + delta_log_mass/2, delta_log_mass)
log_sfr_bins = np.arange(-0.5, 2.5 + delta_log_sfr/2, delta_log_sfr)
    
# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    # args.line_list = ['OIII', 'OII', 'NeIII-3867', 'Hb', 'OIII-4363', 'Ha', 'SII']
    args.line_list = ['OIII', 'OII', 'Hb', 'Ha', 'SII']
    n_lines = len(args.line_list) + 1 # one additional column for the direct image
    args.ids_in_thisbin = args.id
    args.fontfactor = 1.5

    # -----------define colorbar properties-----------
    cmin, cmax, cmap = -19.5, -16.5, 'cividis'
    
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
        output_dir = args.output_dir / args.field / 'stacking'
        if args.adaptive_bins: output_dir = Path(str(output_dir).replace('stacking', 'stacking_adaptive'))
        output_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir / 'plots'
        fig_dir.mkdir(parents=True, exist_ok=True)
        fits_dir = output_dir / 'maps'
        fits_dir.mkdir(parents=True, exist_ok=True)

        # ---------read the photometric catalog file--------------------
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        df['field'] = args.field

        # ---------crossmatch with cosmos-web to get stellar mass and SFR--------------------
        df = get_passage_masses_from_cosmos(df, args, id_col='id')

        # ---------merge with effective radius catalog--------------------
        df_re = Table.read(args.output_dir / f'catalogs/{args.field}_re_list.fits').to_pandas()
        if 'redshift' in df_re: df_re.drop(columns=['redshift'], axis=1, inplace=True)
        df_re['id'] = df_re['id'].astype(int)
        df_re = df_re[df_re['re_kpc'] > 0]
        df = pd.merge(df, df_re, on='id', how='inner')

        # ----------binning the dataframe in an adaptive way----------------
        if args.adaptive_bins:
            final_bins = get_adaptive_bins(df, (log_mass_bins[0],log_mass_bins[-1]), (log_sfr_bins[0],log_sfr_bins[-1]), max_n=args.max_gal_per_bin)

            df['adaptive_bin_id'] = -1
            df['mass_interval'] = None
            df['sfr_interval'] = None
            
            for i, b in enumerate(final_bins):
                mask = (df['log_mass'] >= b['m_min']) & (df['log_mass'] < b['m_max']) & (df['log_sfr'] >= b['s_min']) & (df['log_sfr'] < b['s_max'])
                df.loc[mask, 'adaptive_bin_id'] = i
                
                df.loc[mask, 'mass_interval'] = pd.Interval(left=b['m_min'], right=b['m_max'], closed='left')
                df.loc[mask, 'sfr_interval'] = pd.Interval(left=b['s_min'], right=b['s_max'], closed='left')
                df['bin_intervals'] = list(zip(df['mass_interval'], df['sfr_interval']))

            df = df[df['adaptive_bin_id'] != -1].copy()
            bin_list = pd.unique(df['bin_intervals'])
        
        # ----------binning the dataframe uniformly by mass and SFR bins-------------------
        else:
            df['mass_interval'] = pd.cut(df['log_mass'], bins=log_mass_bins)
            df['sfr_interval'] = pd.cut(df['log_sfr'], bins=log_sfr_bins)
            df = df.dropna(subset=['mass_interval', 'sfr_interval'])
            df['bin_intervals'] = list(zip(df['mass_interval'], df['sfr_interval']))

            all_mass_intervals = df['mass_interval'].cat.categories
            all_sfr_intervals = df['sfr_interval'].cat.categories
            bin_list = list(itertools.product(all_mass_intervals, all_sfr_intervals))
            if args.debug_bin: bin_list = [item for item in bin_list if (item[0].left == 8.5) & (item[0].right == 9.5) & (item[1].left == 1.0) & (item[1].right == 1.5)] # to choose the mass=8.5-9.5, sfr=1-1.5 bin for debugging purposes

        # ------------looping over each bin-----------------------
        list(bin_list).sort(key=lambda x: (x[0].left, x[1].left))
        nbin_good = 0

        for index2, this_mass_sfr_bin in enumerate(bin_list):
            if args.debug_bin and nbin_good > 0: break
            start_time3 = datetime.now()
            bin_text = f'logmassbin_{this_mass_sfr_bin[0].left}-{this_mass_sfr_bin[0].right}_logsfrbin_{this_mass_sfr_bin[1].left}-{this_mass_sfr_bin[1].right}'
            print(f'\tStarting bin ({index2 + 1}/{len(bin_list)}) {bin_text}..', end=' ')

            output_filename = fits_dir / f'stacked_maps_{bin_text}.fits'
            if output_filename.exists() and not args.clobber:
                print(f'Reading existing fits file {output_filename}. Re-run with --clobber to rewrite.')
                # -------reading previously saved stacked fits file------------
                line_dict = read_stacked_maps(output_filename, args)
                stacked_maps, stacked_maps_err, nobj_arr, constituent_ids_array = [], [], [], []
                for index4, this_line in enumerate(args.line_list):
                    if this_line.upper() in line_dict:
                        this_map = line_dict[this_line.upper()]
                        nobj = line_dict[f'{this_line.upper()}_nobj']
                        ids = line_dict[f'{this_line.upper()}_ids']
                        stacked_map, stacked_map_err = unp.nominal_values(this_map.data), unp.std_devs(this_map.data)
                    else:
                        stacked_map, stacked_map_err, nobj, ids = np.nan, np.nan, np.nan, 0
                    
                    stacked_maps.append(stacked_map)
                    stacked_maps_err.append(stacked_map_err)
                    nobj_arr.append(nobj)
                    constituent_ids_array.append(ids)                 
                
                nbin_good += 1
                nobj_good = 0
            else:
                # --------determine which objects fall in this bin----------
                mask = df['bin_intervals'] == this_mass_sfr_bin
                ids_in_thisbin = df.loc[mask, 'id'].tolist()
                if not args.do_all_obj:
                    ids_in_thisbin = np.array(ids_in_thisbin)[np.isin(ids_in_thisbin, args.ids_in_thisbin)] # subset of IDs of "good-looking" galaxies, just for testing

                nobj_good = 0
                ngal_this_bin = len(ids_in_thisbin)
                nrows_total = ngal_this_bin + 1 # one extra row for the stacked plots
                
                if ngal_this_bin > 0:
                    print(f'which has {ngal_this_bin} objects.')
                    line_maps_array = [[] for _ in range(len(args.line_list))]
                    line_maps_err_array = [[] for _ in range(len(args.line_list))] # each array must be len(args.line_list) long and each element will become a stack of 2D images
                    constituent_ids_array = [[] for _ in range(len(args.line_list))] # each array must be len(args.line_list) long and each element will become a stack of 2D images

                    # --------setup PDF for mammoth figures for all emission line map plots----------
                    if not args.debug_align:
                        fullplot_filename_orig = fig_dir / f'binmembers/binmembers_orig_maps_{bin_text}.pdf'
                        fullplot_filename_flux = fig_dir / f'binmembers/binmembers_flux_maps_{bin_text}.pdf'
                        fullplot_filename_err = Path(str(fullplot_filename_flux).replace('flux', 'err'))
                        
                        pdf_orig = PdfPages(fullplot_filename_orig)
                        pdf_flux = PdfPages(fullplot_filename_flux)
                        pdf_err = PdfPages(fullplot_filename_err)
                        n_total_pages = int(np.ceil(nrows_total / args.max_gal_per_page))

                    # ------------looping over chunks of args.max_gal_per_page-----------------------
                    for chunk_id in range(0, nrows_total, args.max_gal_per_page):

                        # --------setup figure for this chunk of the bin----------
                        chunk = ids_in_thisbin[chunk_id : chunk_id + args.max_gal_per_page]
                        if not args.debug_align:
                            n_page = chunk_id // args.max_gal_per_page + 1
                            n_useful_rows_in_page = len(chunk) + (1 if n_page == n_total_pages else 0) # +1 for stack on last page

                            fig_orig, axes_orig = setup_fullpage_figure(n_page, n_total_pages, n_lines,  cmin, cmax, cmap, args, in_re_units=False)
                            fig_flux, axes_flux = setup_fullpage_figure(n_page, n_total_pages, n_lines,  cmin, cmax, cmap, args)
                            fig_err, axes_err = setup_fullpage_figure(n_page, n_total_pages, n_lines,  cmin, cmax, cmap, args)
                        
                        # ----------------looping over the objects in this chunk-------------
                        for index3, this_galaxy in enumerate(chunk):
                            args.id = chunk[index3]
                            print(f'\t\tCommencing ({chunk_id + index3 + 1}/{ngal_this_bin}) ID {args.id}..')
                            obj = df.loc[df['id']==args.id].iloc[0]

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
                            args.distance = cosmo.luminosity_distance(args.z)
                            args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
                            imsize_arcsec = full_hdu['DSCI'].data.shape[0] * args.pix_size_arcsec
                            args.pix_size_re = (2 * args.re_limit) / args.npix_side
                            
                            args.EB_V = 0. # until it gets over-written, if both H alpha and H beta lines are present
                            args.semi_major = obj['a_image']
                            args.semi_minor = obj['b_image']
                            args.pa = obj['theta_image']
                            args.re_kpc, args.re_arcsec = obj['re_kpc'], obj['re_arcsec']
                            
                            # --------determining true center of object rom direct image---------------------
                            args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W', silent=not args.debug_align)

                            # ---------------segmentation map---------------
                            segmentation_map = full_hdu['SEG'].data
                            segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
                            segmentation_map = ndimage.shift(segmentation_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
                            args.segmentation_map = segmentation_map
                            
                            rotated_segmentation_map = rotate_line_map(segmentation_map, args)
                            deprojected_segmentation_map = deproject_line_map(rotated_segmentation_map, args)
                            rescaled_segmentation_map = rescale_line_map(deprojected_segmentation_map, args)

                            # ---------------direct image---------------
                            filter = 'F150W'                    
                            direct_image, exptime = get_direct_image(full_hdu, filter, args) # this is already offset corrected and trimmed
                            rotated_direct_image = rotate_line_map(direct_image, args)
                            deprojected_direct_image = deproject_line_map(rotated_direct_image, args)
                            rescaled_direct_image = rescale_line_map(deprojected_direct_image, args)

                            # ---------plotting direct image and the rescaled direct image------------------------------
                            direct_image_cmap, direct_image_cmin, direct_image_cmax = 'Greys_r', None, None
                            if not args.debug_align:
                                axes_orig[index3, 0] = plot_2D_map(direct_image, axes_orig[index3, 0], f'{filter} OG: {args.id}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False, segmentation_map=segmentation_map, in_re_units=False, seg_col='w')
                                axes_flux[index3, 0] = plot_2D_map(rescaled_direct_image, axes_flux[index3, 0], f'{filter}_rs: {args.id}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False)
                                axes_err[index3, 0] = plot_2D_map(rescaled_direct_image, axes_err[index3, 0], f'{filter}: {args.id}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False)

                            # ----------plotting the direct image: for debugging--------------
                            if args.debug_align:
                                re_pix = args.re_arcsec / args.pix_size_arcsec
                                fig_debug, axes_debug_2d = plt.subplots(2, len(args.line_list) + 1, figsize=(13, 6))
                                fig_debug.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)

                                axes = axes_debug_2d[0]
                                axes[0] = myimshow(direct_image, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=direct_image_cmap, col='w')
                                axes[1].set_visible(False) # no smoothed direct image
                                axes[2] = myimshow(rotated_direct_image, axes[2], contour=rotated_segmentation_map != args.id, re_pix=re_pix, label='Rotated', cmap=direct_image_cmap, col='w')
                                axes[3] = myimshow(deprojected_direct_image, axes[3], contour=deprojected_segmentation_map != args.id, re_pix=re_pix, label='Deprojected', cmap=direct_image_cmap, col='w')
                                axes[4] = myimshow(rescaled_direct_image, axes[4], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=direct_image_cmap, col='w')

                            # ---------looping over all emission lines------------------------------
                            get_rescaled_recentering = True
                            nlines_good = 0
                            for index4, this_line in enumerate(args.line_list):
                                if this_line not in args.available_lines:
                                    print(f'\t\t\tLine {this_line} not available for object {args.id}. So skipping this line..')
                                    if not args.debug_align:
                                        axes_flux[index3, index4 + 1].remove()
                                        axes_err[index3, index4 + 1].remove()
                                    continue
                                try:
                                    # -----------extracting the emission line map----------------
                                    line_map, line_map_err = get_emission_line_map(this_line, full_hdu, args, dered=True, silent=True)

                                    # ------smoothing the map before rotating--------
                                    smoothing_kernel = Box2DKernel(args.kernel_size, mode=args.kernel_mode)
                                    smoothed_line_map = convolve(line_map, smoothing_kernel)
                                    smoothed_line_map_err = convolve(line_map_err, smoothing_kernel)

                                    # --------rotating the line map---------------------
                                    #rotated_line_map = rotate_line_map(line_map, args)
                                    #rotated_line_map_err = rotate_line_map(line_map_err, args)
                                    rotated_line_map = rotate_line_map(smoothed_line_map, args)
                                    rotated_line_map_err = rotate_line_map(smoothed_line_map_err, args)

                                    # --------deprojecting the line map---------------------
                                    deprojected_line_map = deproject_line_map(rotated_line_map, args)
                                    deprojected_line_map_err = deproject_line_map(rotated_line_map_err, args)

                                    # --------scaling the line map by effective radius---------------------
                                    rescaled_line_map = rescale_line_map(deprojected_line_map, args)
                                    rescaled_line_map_err = rescale_line_map(deprojected_line_map_err, args)

                                    # --------computing new recentering offsets, after rescaling---------------------
                                    if get_rescaled_recentering:
                                        args.ndelta_xpix_rescaled, args.ndelta_ypix_rescaled = get_center_offsets(rescaled_line_map, args, silent=not args.debug_align)
                                        recentered_segmentation_map = ndimage.shift(rescaled_segmentation_map, [args.ndelta_xpix_rescaled, args.ndelta_ypix_rescaled], order=0, cval=np.nan)
                                        get_rescaled_recentering = False

                                        if args.debug_align:
                                            recentered_direct_image = ndimage.shift(rescaled_direct_image, [args.ndelta_xpix_rescaled, args.ndelta_ypix_rescaled], order=0, cval=np.nan)
                                            axes[5] = myimshow(recentered_direct_image, axes[5], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered', cmap=direct_image_cmap, col='w')
                                    
                                    # --------recentering the line map---------------------
                                    recentered_line_map = ndimage.shift(rescaled_line_map, [args.ndelta_xpix_rescaled, args.ndelta_ypix_rescaled], order=0, cval=np.nan)
                                    recentered_line_map_err = ndimage.shift(rescaled_line_map_err, [args.ndelta_xpix_rescaled, args.ndelta_ypix_rescaled], order=0, cval=np.nan)

                                    # -------plotting the intermediate steps: for debugging--------------
                                    if args.debug_align:
                                        fig_debug.text(0.05, 0.98, f'{args.field}: ID {args.id}: {this_line} map: Alignment diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')
                                        print(f'Deb345: shapes: {this_line} map = {np.shape(line_map)}, rotated = {np.shape(rotated_line_map)}, deprojected = {np.shape(deprojected_line_map)}, recaled = {np.shape(rescaled_line_map)}, redshift={args.z:.1f}, log_mass={obj["log_mass"]:.1f}, log_sfr={obj["log_sfr"]:.1f}, a={args.semi_major:.1f}, b={args.semi_minor:.1f}, a/b={args.semi_major/args.semi_minor:.1f}, pa={args.pa:.1f}, re={re_pix:.1f} pixels') ##
                                        
                                        axes = axes_debug_2d[1]
                                        axes[0] = myimshow(line_map, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=cmap, col='k')
                                        axes[1] = myimshow(smoothed_line_map, axes[1], contour=segmentation_map != args.id, re_pix=re_pix, label='Smoothed', cmap=cmap, col='k')
                                        axes[2] = myimshow(rotated_line_map, axes[2], contour=rotated_segmentation_map != args.id, re_pix=re_pix, label='Rotated', cmap=cmap, col='k')
                                        axes[3] = myimshow(deprojected_line_map, axes[3], contour=deprojected_segmentation_map != args.id, re_pix=re_pix, label='Deprojected', cmap=cmap, col='k')
                                        axes[4] = myimshow(rescaled_line_map, axes[4], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=cmap, col='k')
                                        axes[5] = myimshow(recentered_line_map, axes[5], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered', cmap=cmap, col='k')

                                        figname = fig_dir / f'debug_align_{args.id}.png'
                                        fig_debug.savefig(figname, transparent=args.fortalk)
                                        print(f'\nSaved figure as {figname}')

                                        plt.show(block=False)
                                        sys.exit(f'Exiting here because of --debug_align mode; if you want to run the full code as usual then remove the --debug_align option and re-run')

                                    # ---------------appending the line map-------------------------
                                    line_maps_array[index4].append(recentered_line_map)
                                    line_maps_err_array[index4].append(recentered_line_map_err)
                                    constituent_ids_array[index4].append(f'{args.field}-{args.id}')
                                    nlines_good += 1

                                    # -------------plotting this line map of this galaxy-------------------
                                    axes_orig[index3, index4 + 1] = plot_2D_map(line_map, axes_orig[index3, index4 + 1], f'{this_line}: {args.id} OG', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True, segmentation_map=segmentation_map, in_re_units=False)
                                    axes_flux[index3, index4 + 1] = plot_2D_map(recentered_line_map, axes_flux[index3, index4 + 1], f'{this_line}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True)
                                    axes_err[index3, index4 + 1] = plot_2D_map(recentered_line_map_err, axes_err[index3, index4 + 1], f'{this_line}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True)
                                
                                except Exception as e:
                                    print(f'Skipping {this_line} for obj {args.id} because it failed due to: {e}')
                                    continue
                                
                            print(f'\t\t\tFound {nlines_good} lines for {args.id}')
                            if nlines_good > 0: nobj_good += 1
                    
                        # -------------removing unused rows-------------------------
                        for r in range(n_useful_rows_in_page, args.max_gal_per_page):
                            for c in range(n_lines):
                                axes_orig[r, c].set_visible(False)
                                axes_flux[r, c].set_visible(False)
                                axes_err[r, c].set_visible(False)

                        # -------------finalising saving each page of the mammoth PDFs (unless last page)-------------------------
                        if n_page < n_total_pages:
                            pdf_orig.savefig(fig_orig)
                            pdf_flux.savefig(fig_flux)
                            pdf_err.savefig(fig_err)
                   
                    # -----stacking all line maps for all objects in this bin-----------
                    stacked_maps, stacked_maps_err, nobj_arr = [], [], []
                    for index4, this_line in enumerate(args.line_list):
                        stacked_map, stacked_map_err = weighted_stack_line_maps(line_maps_array[index4], line_maps_err_array[index4])
                        stacked_maps.append(stacked_map)
                        stacked_maps_err.append(stacked_map_err)
                        nobj_arr.append(len(constituent_ids_array[index4]))

                        # --------displaying stacked maps at the bottom of the mammoth figure---------
                        curr_row = (nrows_total % args.max_gal_per_page) - 1
                        if np.ndim(stacked_map) == 2:
                            axes_orig[curr_row, index4 + 1] = plot_2D_map(stacked_map, axes_orig[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0)
                            axes_flux[curr_row, index4 + 1] = plot_2D_map(stacked_map, axes_flux[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0)
                            axes_err[curr_row, index4 + 1] = plot_2D_map(stacked_map_err, axes_err[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0)

                    # -------writing out stacked line maps as fits files--------------
                    if args.do_all_obj: write_stacked_maps(stacked_maps, stacked_maps_err, constituent_ids_array, output_filename, args)
                    
                    # ------------finalising the mammoth PDFs------------
                    pdf_orig.savefig(fig_orig)
                    pdf_flux.savefig(fig_flux)
                    pdf_err.savefig(fig_err)

                    pdf_orig.close()
                    pdf_flux.close()
                    pdf_err.close()

                    if not args.debug_align: plt.close('all')                   
                    print(f'Saved {fullplot_filename_flux}')

                else:
                    print(f'which has no object. Skipping this bin.')
                    continue

            print(f'\nCompleted bin mass={this_mass_sfr_bin[0]}, sfr={this_mass_sfr_bin[1]} ({nobj_good} / {len(nobj_arr)} objects) in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')
            if nobj_good > 1: nbin_good += 1

        print(f'Completed field {field} ({nbin_good} / {len(bin_list)} bins) in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
