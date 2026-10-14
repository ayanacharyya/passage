'''
    Filename: stack_emission_maps.py
    Notes: Stacks 2D emission line maps (from grizli) for objects in a given field/s, in bins of mass and/or SFR
    Author : Ayan
    Created: 15-01-26
    Example: run stack_emission_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28 --do_all_obj
             run stack_emission_maps.py --field Par28 --debug_align --id 2822,2698,2583,2171,672
             run stack_emission_maps.py --field Par28 --do_all_obj --clobber --debug_bin --skip_deproject --skip_re_scaling
             run stack_emission_maps.py --field Par28 --do_all_obj --re_limit 2
             run stack_emission_maps.py --system ssd --field Par28 --do_all_obj --adaptive_bins --max_gal_per_bin 20
             run stack_emission_maps.py --system ssd --field Par28 --do_all_obj --adaptive_bins --bin_by_distance
             run stack_emission_maps.py --system ssd --field Par28 --do_all_obj --adaptive_bins --bin_by_distance_mass
             run stack_emission_maps.py --system ssd --do_all_fields --do_all_obj --adaptive_bins --bin_by_distance
             run stack_emission_maps.py --system ssd --field Par28 --id 1982 --clobber --debug_bin --debug_align --adaptive_bins --bin_by_distance_mass
             run stack_emission_maps.py --system ssd --do_all_fields --do_all_obj --adaptive_bins --bin_by_distance_mass
             run stack_emission_maps.py --system ssd --do_all_fields --do_all_obj --cut_z_flag 4 --adaptive_bins --bin_by_sfh_mass --skip_deproject
'''

from header import *
from util import *
setup_plot_style()
from make_diagnostic_maps import trim_image, get_dereddened_flux, myimshow, get_offsets_from_center, get_cutout
from make_sfms_bins import get_binned_df, required_lines, sfms

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
        image = trim_image(image, args, skip_re_trim=True)
        exptime = 1
    except:
        try:
            hdu = full_hdu['DSCI', f'{filter}-CLEAR']
            image = hdu.data
            image = trim_image(image, args, skip_re_trim=True)
            exptime = full_hdu[0].header[f'T_{filter.upper()}']
        except:
            try:
                hdu = full_hdu['DSCI', f'{filter.upper()}-{filter.upper()}-CLEAR']
                image = hdu.data
                image = trim_image(image, args, skip_re_trim=True)
                exptime = full_hdu[0].header[f'T_{filter.upper()}']
            except:
                try:
                    image = full_hdu['DSCI', 'CLEAR'].data
                    image = trim_image(image, args, skip_re_trim=True)
                    filter = 'CLEAR'
                    exptime = 1
                except:
                    full_field_filename = args.input_dir / f'{args.field}' / 'Products' / f'{args.field}_{filter.lower()}-clear_drz_sci.fits'
                
                    exptime = fits.open(full_field_filename)[0].header['EXPTIME']
                    print(f'{filter.upper()} not found in full_hdu extension. Therefore trying to get cutout from full field image {full_field_filename}')
                    pos = SkyCoord(full_hdu[0].header['RA'], full_hdu[0].header['DEC'], unit = 'deg')
                    size = 2 * args.arcsec_limit * u.arcsec
                    try: target_header = full_hdu['DSCI', 'F140W'].header
                    except: target_header = full_hdu['DSCI', 'CLEAR'].header
                    
                    temp1, temp2 = args.only_seg, args.vorbin
                    args.only_seg, args.vorbin = False, False
                    image = get_cutout(full_field_filename, pos, size, target_header, args, plot_test_axes=None, skip_re_trim=True)
                    args.only_seg, args.vorbin = temp1, temp2

    if image is not None:
        image = ndimage.shift(image, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
    
    return image, exptime

# --------------------------------------------------------------------------------------------------------------------
def correct_ha_F18(N2_plus_Ha_map, args):
    '''
    Extract the H-alpha flux from the N2+Halpha compund, following F18 N2/Ha calibration
    '''
    def F18_func(x, redshift, log_mass): # Eq 3 of Faisst+2018
        psi = x + 0.138 - 0.042 * (1 + redshift) **2 
        solve = 3.696 * psi + 3.236 * psi**(-1) + 0.729 * psi ** (-2) + 14.928 + 0.156 * (1 + redshift)**2 - log_mass
        return solve

    log_N2Ha = brentq(F18_func, -3, 0, args=(args.z, args.log_mass))
    N2Ha = 10 ** log_N2Ha
    Ha_map = N2_plus_Ha_map / (1 + 1.333 * N2Ha) # 1.333 factor is because F18 consider BOTH NII components, assuming NII 6548 = N II 6584 / 3

    if hasattr(N2_plus_Ha_map, "__len__"): # if it is an array
        Ha_map[Ha_map==0] = np.nan

    return Ha_map

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
        if line == 'OIII': # special treatment for OIII 5007 line
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            if not silent: print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'Ha': # special treatment for Ha line
            line_map = correct_ha_F18(line_map, args)
            line_map_err = correct_ha_F18(line_map_err, args)
            if not silent: print(f'Correcting Ha for NII component, using Faisst+2018 calibration based on galaxy mass and redshift')
            
    line_map = line_map * factor
    line_map_err = line_map_err * factor

    # ----------getting a smaller cutout around the object center-----------
    line_map = trim_image(line_map, args, skip_re_trim=True)
    line_map_err = trim_image(line_map_err, args, skip_re_trim=True)

    # ----------re-center line map-----------
    line_map = ndimage.shift(line_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
    line_map_err = ndimage.shift(line_map_err, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)

    # -------masking outside segmentation map (using NaNs for scipy compatibility)---------
    mask = (args.segmentation_map != args.id)
    line_map[mask] = np.nan
    line_map_err[mask] = np.nan

    # -----------getting the dereddened flux value-----------------
    if dered and args.EB_V != 0:
        line_map_quant = get_dereddened_flux(unp.uarray(line_map, line_map_err), line_wave, args.EB_V)
        line_map = unp.nominal_values(line_map_quant)
        line_map_err = unp.std_devs(line_map_quant)

    # ---------converting from flux to surface brightness units------------
    pixscale_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value # kpc
    line_map /= (pixscale_kpc ** 2) # in ergs/s/cm^2/kpc^2
    line_map_err /= (pixscale_kpc ** 2)

    # ---------getting integrated flux-----------------
    line_int = ufloat(np.nansum(line_map), np.sqrt(np.nansum(line_map_err ** 2)))

    return line_map, line_map_err, line_int

# --------------------------------------------------------------------------------------------------------------------
def transform_map_and_err(line_map, line_map_err, args):
    '''
    Transforms flux map using order=3 cubic spline and error map via variance (order=1).
    Applies dimensional scaling by Re_kpc^2 when stacking in dimensionless Re-space.
    '''
    # Dimensionless Re-space scaling: scale flux by Re_kpc^2
    if not args.skip_re_scaling:
        line_map = line_map * (args.re_kpc ** 2)
        line_map_err = line_map_err * (args.re_kpc ** 2)

    transformed_map = transform_map_single_pass(line_map, args)
    
    # Resample in variance space (err^2) to preserve statistical independence
    variance_map = line_map_err ** 2
    transformed_var = transform_map_single_pass(variance_map, args)
    transformed_err = np.sqrt(np.maximum(transformed_var, 0.0))

    return transformed_map, transformed_err

# --------------------------------------------------------------------------------------------------------------------
def transform_map_single_pass(image, args, is_seg=False):
    '''
    Combines rotation, axis deprojection, Re-rescaling, and shifts into a single affine transform.
    Executes in a single pass to eliminate multi-step low-pass blurring.
    '''
    out_shape = (args.npix_side, args.npix_side)
    
    # Input image was already centered by (ndelta_xpix, ndelta_ypix) in get_direct_image / get_emission_line_map
    shift_in = (0, 0)
    shift_out = (getattr(args, 'ndelta_xpix_direct_rescaled', 0), getattr(args, 'ndelta_ypix_direct_rescaled', 0))
    stretch_factor = 1.0 if args.skip_deproject else args.semi_major / args.semi_minor # do not do the minor/major axis ratio correction if deproject is turned off
    
    if args.skip_re_scaling:
        target_kpc_px = args.npix_side / (2 * args.kpc_limit)
        current_kpc_px = 1.0 / args.pix_size_kpc
        zoom_factor = target_kpc_px / current_kpc_px
    else:
        target_re_px = args.npix_side / (2 * args.re_limit)
        current_re_px = args.re_arcsec / args.pix_size_arcsec
        zoom_factor = target_re_px / current_re_px

    in_center = np.array(image.shape) / 2.0 - 0.5
    out_center = np.array(out_shape) / 2.0 - 0.5

    alpha = np.radians(90) + args.pa # obj['theta_image'] and therefore args.pa is already in radians
    # also in sky convention if PA is negative then the major is rotated CLOCKWISE from the horizontal line by that amount, therefore
    # to fix that we need to rotate COUNTER-CLOCKWISE by that amount, which means we use the same negative value of PA as alpha
    # because for positive alpha the code below rotates clockwise and for negative alpha it rotates anticlockwise
    # essentially, if we were to align the galaxy back to major axis as horizontal, we would leave alpha=PA because the difference 
    # in conventions do the job of reversing the sign for us atuomatically; however, we want to align the galaxy such that major axis is vertical
    # we therefore rotate by an amount which is 'how much is the galaxy;s major axis off from vertical in the anticlockwise direction (because that is positive PA)?'
    # which is equal to 90 + PA
    cos_a, sin_a = np.cos(alpha), np.sin(alpha)
    R = np.array([[cos_a, sin_a], [-sin_a, cos_a]])

    sy = zoom_factor * stretch_factor
    sx = zoom_factor
    S = np.array([[1.0 / sy, 0.0], [0.0, 1.0 / sx]])

    M = R @ S
    offset = (in_center - np.array(shift_in)) - M @ (out_center + np.array(shift_out))

    order = 0 if is_seg else 1
    transformed = ndimage.affine_transform(image, M, offset=offset, output_shape=out_shape, order=order, cval=np.nan)
    return transformed

# --------------------------------------------------------------------------------------------------------------------
def get_center_offsets_gemini(dir_img, args, silent=False):
    '''
    Accepts a 2D map, then computes the offset from the original center of image to the brightest pixel.
    Returns two integers (offset in x and y axes).
    '''
    smoothing_kernel = Box2DKernel(5, mode=args.kernel_mode)
    dir_img_smoothed = convolve(dir_img, smoothing_kernel)
    smooth_shape = dir_img_smoothed.shape
    ncells = 5
    
    cen_y, cen_x = smooth_shape[0] // 2, smooth_shape[1] // 2

    y_min = max(cen_y - ncells, 0)
    y_max = min(cen_y + ncells, smooth_shape[0])
    x_min = max(cen_x - ncells, 0)
    x_max = min(cen_x + ncells, smooth_shape[1])
    
    dir_img_smoothed_subarea = dir_img_smoothed[y_min:y_max, x_min:x_max]

    # Handle cases where subarea is all NaNs gracefully
    if np.all(np.isnan(dir_img_smoothed_subarea)):
        if not silent: print(f'For {args.field}:{args.id}: Direct image subarea is all NaNs, returning 0 offset')
        return 0, 0

    max_idx = np.nanargmax(dir_img_smoothed_subarea)
    rel_y, rel_x = np.unravel_index(max_idx, dir_img_smoothed_subarea.shape)

    brightest_y = y_min + rel_y
    brightest_x = x_min + rel_x

    ndelta_xpix = cen_x - brightest_x
    ndelta_ypix = cen_y - brightest_y

    if not silent: print(f'For {args.field}:{args.id}: Determined x and y offsets from direct image = {ndelta_xpix}, {ndelta_ypix}')

    return ndelta_xpix, ndelta_ypix

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
    brightest_coords = np.where(dir_img_smoothed == np.nanmax(dir_img_smoothed_subarea))
    brightest_x, brightest_y = brightest_coords[0][0], brightest_coords[1][0]
    cen_x, cen_y = int(np.shape(dir_img)[0] / 2), int(np.shape(dir_img)[1] / 2)
    ndelta_xpix = cen_x - brightest_x
    ndelta_ypix = cen_y - brightest_y
    if not silent: print(f'For {args.field}:{args.id}: Determined x and y offsets from direct image = {ndelta_xpix}, {ndelta_ypix}')

    return ndelta_xpix, ndelta_ypix

# --------------------------------------------------------------------------------------------------------------------
def weighted_stack_line_maps_old(line_maps_array, line_maps_err_array, additional_weights=None):
    '''
    Accepts N emission line maps and error maps and computes a weighted stack.
    Normalizes each galaxy's pixel weights by its spatial mean so that pixel-level S/N downweighting is preserved 
    while giving each galaxy equal total weight across the stack.
    '''
    line_maps_array = np.array(line_maps_array)
    line_maps_err_array = np.array(line_maps_err_array)
    ngal = line_maps_array.shape[0]

    # Pixel inverse variance
    raw_weights = 1.0 / (line_maps_err_array ** 2)

    # Normalize per galaxy: mean weight across valid pixels = 1
    normalized_weights = np.zeros_like(raw_weights)
    for i in range(ngal):
        valid_mask = ~np.isnan(raw_weights[i]) & ~np.isnan(line_maps_array[i])
        if np.any(valid_mask):
            mean_w = np.nanmean(raw_weights[i][valid_mask])
            if mean_w > 0:
                normalized_weights[i] = raw_weights[i] / mean_w

    if additional_weights is not None:
        add_w = np.array(additional_weights)[:, np.newaxis, np.newaxis]
        normalized_weights *= add_w

    weights_sum = np.nansum(normalized_weights, axis=0)
    stacked_line_map = np.nansum(line_maps_array * normalized_weights, axis=0) / weights_sum

    # Correct error propagation for normalized weighted mean: sigma_stack = sqrt(sum(w_i^2 * sigma_i^2)) / sum(w_i)
    variance_term = np.nansum((normalized_weights ** 2) * (line_maps_err_array ** 2), axis=0)
    stacked_line_map_err = np.sqrt(variance_term) / weights_sum

    return stacked_line_map, stacked_line_map_err

# --------------------------------------------------------------------------------------------------------------------
def weighted_stack_line_maps(line_maps_array, line_maps_err_array, master_weights):
    '''
    Accepts N emission line maps and error maps for a single line,
    along with a pre-computed N x Y x X master weight array shared across all lines.
    '''
    line_maps_array = np.array(line_maps_array)
    line_maps_err_array = np.array(line_maps_err_array)
    
    # Mask invalid pixels across line flux, error, or weights
    valid_mask = np.isfinite(line_maps_array) & np.isfinite(line_maps_err_array) & np.isfinite(master_weights)
    
    w = np.where(valid_mask, master_weights, 0.0)
    weights_sum = np.nansum(w, axis=0)
    
    # Avoid division by zero for unmapped pixels
    weights_sum_safe = np.where(weights_sum == 0, np.nan, weights_sum)
    
    # Stacked flux map
    stacked_line_map = np.nansum(line_maps_array * w, axis=0) / weights_sum_safe

    # Correct error propagation using shared weights: sigma_stack = sqrt(sum(w_i^2 * sigma_{i, line}^2)) / sum(w_i)
    variance_term = np.nansum((w ** 2) * (line_maps_err_array ** 2), axis=0)
    stacked_line_map_err = np.sqrt(variance_term) / weights_sum_safe

    return stacked_line_map, stacked_line_map_err

# --------------------------------------------------------------------------------------------------------------------
def quickplot_raw_line_map(objid, field='Par028', line='OIII', arcsec_limit=1, no_plot=False, cmin=None, cmax=None, takelog=False, col='w'):
    filename = args.input_dir / field / 'Products'/ 'maps' / f'{field}_{objid:05d}.maps.fits'
    full_hdu = fits.open(filename)
    pix_size_arcsec = full_hdu[5].header['PIXASEC']
    available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
    segmentation_map = full_hdu['SEG'].data

    if line in available_lines:
        line_hdu = full_hdu['LINE', line]
        line_map = line_hdu.data * 1e-17
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
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['CONTENT'] = 'Stacked Galaxy Emission Line Maps'
    primary_hdu.header['NPIX'] = (args.npix_side, 'Pixels per side')
    pix_per_re = args.npix_side / (2 * args.re_limit)
    primary_hdu.header['RE_PX'] = (pix_per_re, f'Scale: 1 Re = {pix_per_re} pixels')
    
    hdu_list = [primary_hdu]

    for index, this_line in enumerate(args.line_list):
        if type(stacked_maps[index]) != np.ndarray or type(stacked_maps_err[index]) != np.ndarray:
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
def plot_2D_map(image, ax, label, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True, segmentation_map=None, in_re_units=True, in_kpc_units=False, seg_col='k'):
    if takelog: image = np.log10(image.data) if hasattr(image, 'data') else np.log10(image)

    if in_kpc_units:
        offset = args.pix_size_kpc / 2
        args.extent = (-args.kpc_limit - offset, args.kpc_limit - offset, -args.kpc_limit - offset, args.kpc_limit - offset)
    elif in_re_units:
        offset = args.pix_size_re / 2
        args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)
    else:
        offset = args.pix_size_arcsec / 2
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
def setup_fullpage_figure(n_page, n_total_pages, n_lines, cmin, cmax, cmap, args, in_re_units=True, in_kpc_units=False):
    fig, axes = plt.subplots(args.max_gal_per_page, n_lines, figsize=(1.5 * n_lines, 1.5 * args.max_gal_per_page))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.995, bottom=0.08, wspace=0., hspace=0.)

    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    cmap = plt.get_cmap(cmap)

    sm = mpl_cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, location='top', shrink=0.95, pad=0.01, aspect=60)
    cbar.set_label('Surface brightness [ergs/s/cm^2/kpc^2]', fontsize=args.fontsize, labelpad=5)    
    cbar.ax.tick_params(labelsize=args.fontsize / args.fontfactor)

    units = 'kpc' if in_kpc_units else r'R$_e$' if in_re_units else r'arcseconds'
    label = f'Offset ({units})'
    fig.supxlabel(label, fontsize=args.fontsize)
    fig.supylabel(label, fontsize=args.fontsize)

    page_text = f'Page {n_page} of {n_total_pages}'
    fig.text(0.95, 0.01, page_text, ha='right', fontsize=args.fontsize / args.fontfactor, color='k')

    return fig, axes

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.ids_in_thisbin = args.id
    args.fontfactor = 1.5

    # -----------define colorbar properties-----------
    filter_for_direct_image = 'F150W'
    direct_image_cmin, direct_image_cmax, direct_image_cmap = -3, -1, 'Greys_r'
    cmin, cmax, cmap = -3., 0., 'cividis'
    
    # ------------reading and binning dataframe-------------
    df, bin_list, args = get_binned_df(args, required_lines=required_lines, sfms=sfms)
    n_lines = len(args.line_list) + 1

    if args.id is not None: df = df[df['id'].isin(args.id)]
    
    nbin_good = 0
    nobj_total_binned = 0

    # ------------looping over each bin-----------------------
    for index2, this_mass_sfr_bin in enumerate(bin_list):
        start_time3 = datetime.now()
        if args.bin_by_distance:
            bin_text = f'delta_sfms_bin_{this_mass_sfr_bin.left}-{this_mass_sfr_bin.right}'
        elif args.bin_by_distance_mass:
            this_delta_sfms_bin = this_mass_sfr_bin[0]
            this_mass_bin = this_mass_sfr_bin[1]
            bin_text = f'delta_sfms_bin_{this_delta_sfms_bin.left}-{this_delta_sfms_bin.right}_mass_bin_{this_mass_bin.left}-{this_mass_bin.right}'
        elif args.bin_by_sfh_mass:
            this_delta_sfms_bin = this_mass_sfr_bin[0]
            this_mass_bin = this_mass_sfr_bin[1]
            bin_text = f'tform_ratio_bin_{this_delta_sfms_bin.left}-{this_delta_sfms_bin.right}_mass_bin_{this_mass_bin.left}-{this_mass_bin.right}'
        else:
            bin_text = f'logmassbin_{this_mass_sfr_bin[0].left}-{this_mass_sfr_bin[0].right}_logsfrbin_{this_mass_sfr_bin[1].left}-{this_mass_sfr_bin[1].right}'
        print(f'\tStarting bin ({index2 + 1}/{len(bin_list)}) {bin_text}..', end=' ')

        output_filename = args.fits_dir / f'stacked_maps{args.deproject_text}{args.rescale_text}_{bin_text}.fits'
        if output_filename.exists() and not args.clobber:
            # -------reading previously saved stacked fits file------------
            print(f'Reading existing fits file {output_filename}. Re-run with --clobber to rewrite.')
            line_dict = read_stacked_maps(output_filename, args)
            stacked_maps, stacked_maps_err, constituent_ids_array = [], [], []
            for index4, this_line in enumerate(args.line_list):
                if this_line.upper() in line_dict:
                    this_map = line_dict[this_line.upper()]
                    ids = line_dict[f'{this_line.upper()}_ids']
                    stacked_map, stacked_map_err = unp.nominal_values(this_map.data), unp.std_devs(this_map.data)
                else:
                    stacked_map, stacked_map_err, ids = np.nan, np.nan, [0]
                
                stacked_maps.append(stacked_map)
                stacked_maps_err.append(stacked_map_err)
                constituent_ids_array.append(ids)                 
            
            nbin_good += 1
            nobj_good = 0
            ngal_this_bin = len(ids)
            nobj_no_scale_line = 0
        else:
            # --------determine which objects fall in this bin----------
            if args.bin_by_distance_mass or args.bin_by_sfh_mass:
                mask = (df['bin_intervals'] == this_delta_sfms_bin) & (df['mass_intervals'] == this_mass_bin)
            else:
                mask = df['bin_intervals'] == this_mass_sfr_bin
            
            ids_in_thisbin = df.loc[mask, 'id'].tolist()
            if not args.do_all_obj:
                ids_in_thisbin = np.array(ids_in_thisbin)[np.isin(ids_in_thisbin, args.ids_in_thisbin)]

            nobj_good = 0
            ngal_this_bin = len(ids_in_thisbin)
            nrows_total = ngal_this_bin + 1
            nobj_no_scale_line = 0
            
            if ngal_this_bin > 0:
                print(f'which has {ngal_this_bin} objects.')
                line_maps_array = [[] for _ in range(len(args.line_list))]
                line_maps_err_array = [[] for _ in range(len(args.line_list))]
                constituent_ids_array = [[] for _ in range(len(args.line_list))]
                master_weight_dict = {}

                # --------setup PDF for mammoth figures for all emission line map plots----------
                if not args.debug_align:
                    fullplot_filename_orig = args.fig_dir / f'binmembers/binmembers_orig_maps{args.deproject_text}{args.rescale_text}_{bin_text}.pdf'
                    fullplot_filename_flux = args.fig_dir / f'binmembers/binmembers_flux_maps{args.deproject_text}{args.rescale_text}_{bin_text}.pdf'
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
                        n_useful_rows_in_page = len(chunk) + (1 if n_page == n_total_pages else 0)

                        fig_orig, axes_orig = setup_fullpage_figure(n_page, n_total_pages, n_lines, -20, -17, cmap, args, in_re_units=False)
                        fig_flux, axes_flux = setup_fullpage_figure(n_page, n_total_pages, n_lines, cmin, cmax, cmap, args, in_kpc_units=args.skip_re_scaling)
                        fig_err, axes_err = setup_fullpage_figure(n_page, n_total_pages, n_lines, cmin, cmax, cmap, args, in_kpc_units=args.skip_re_scaling)
                    
                    # ----------------looping over the objects in this chunk-------------
                    for index3, this_galaxy in enumerate(chunk):
                        obj = df.loc[df['id'] == this_galaxy].iloc[0]
                        args.field = obj['field']
                        args.id = obj['id']
                        print(f'\n\t\tCommencing ({chunk_id + index3 + 1}/{ngal_this_bin}) ID {args.field}:{this_galaxy}..')

                        # ------determining directories and filenames---------
                        product_dir = args.input_dir / args.field / 'Products'
                        full_fits_file = product_dir / 'full_psf_matched' / f'{args.field}_{args.id:05d}.full.fits'
                        maps_fits_file = product_dir / 'maps_psf_matched' / f'{args.field}_{args.id:05d}.maps.fits'
                        od_filename = product_dir / 'spec1D' / f'{args.field}_{args.id:05d}.1D.fits'

                        if os.path.exists(maps_fits_file):
                            full_filename = maps_fits_file
                        elif os.path.exists(full_fits_file):
                            full_filename = full_fits_file
                        else:
                            print(f'Could not find {full_fits_file} or {maps_fits_file}, so skipping it.')
                            continue
                        if not os.path.exists(od_filename):
                            od_filename = Path(str(od_filename).replace('.1D.', '.spec1D.'))

                        # ------------read in maps files--------------------------------
                        full_hdu = fits.open(full_filename)

                        # ----------determining object parameters------------
                        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
                        args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines])
                        args.z = full_hdu[0].header['REDSHIFT']
                        args.distance = cosmo.luminosity_distance(args.z)
                        args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
                        args.pix_size_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value
                        args.pix_size_re = (2 * args.re_limit) / args.npix_side
                        
                        args.EB_V = 0.
                        args.semi_major = obj['a_image']
                        args.semi_minor = obj['b_image']
                        args.pa = obj['theta_image']
                        args.re_kpc, args.re_arcsec = obj['re_kpc'], obj['re_arcsec']
                        args.log_mass = obj['log_mass']
                        
                        # --------determining true center of object rom direct image---------------------
                        args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter=filter_for_direct_image, silent=not args.debug_align)

                        # ---------------segmentation map---------------
                        segmentation_map = full_hdu['SEG'].data
                        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
                        segmentation_map = ndimage.shift(segmentation_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
                        args.segmentation_map = segmentation_map
                        
                        rescaled_segmentation_map = transform_map_single_pass(segmentation_map, args, is_seg=True)

                        # -----------extracting integrated OIII line flux, for scaling----------------
                        if args.scaling_line in args.available_lines:
                            _, _, line_int = get_emission_line_map(args.scaling_line, full_hdu, args, dered=False, silent=True)
                            integrated_scaling_flux = obj[f'flux_{args.scaling_line}']
                            scaling_flux_err = obj.get(f'flux_err_{args.scaling_line}', 0.0)
                            if scaling_flux_err is None or np.isnan(scaling_flux_err):
                                scaling_flux_err = 0.0

                            if integrated_scaling_flux > 0:
                                print(f'Scaling all line fluxes of {args.id} by {args.scaling_line} flux = {integrated_scaling_flux}')
                            else:
                                print(f'Since integrated flux for {args.scaling_line} < 0 (= {integrated_scaling_flux}), for {args.id}, not putting this in stack')
                                nobj_no_scale_line += 1
                                continue
                        else:
                            print(f'Since {args.scaling_line} not available for {args.id}, not putting this in stack')
                            nobj_no_scale_line += 1
                            continue

                        # ---------------direct image---------------
                        direct_image, exptime = get_direct_image(full_hdu, filter_for_direct_image, args)
                        rescaled_direct_image = transform_map_single_pass(direct_image, args)

                        args.ndelta_xpix_direct_rescaled, args.ndelta_ypix_direct_rescaled = get_center_offsets(rescaled_direct_image, args, silent=not args.debug_align)
                        recentered_direct_image = ndimage.shift(rescaled_direct_image, [args.ndelta_xpix_direct_rescaled, args.ndelta_ypix_direct_rescaled], order=0, cval=np.nan)
                        recentered_segmentation_map = ndimage.shift(rescaled_segmentation_map, [args.ndelta_xpix_direct_rescaled, args.ndelta_ypix_direct_rescaled], order=0, cval=np.nan)

                        # ---------plotting direct image and the rescaled direct image------------------------------
                        if not args.debug_align:
                            axes_orig[index3, 0] = plot_2D_map(direct_image, axes_orig[index3, 0], f'{filter_for_direct_image} OG: {args.id}', args, cmap=direct_image_cmap, takelog=True, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False, segmentation_map=segmentation_map, in_re_units=False, seg_col='w')
                            axes_flux[index3, 0] = plot_2D_map(rescaled_direct_image, axes_flux[index3, 0], f'{filter_for_direct_image}_rs: {args.id}', args, cmap=direct_image_cmap, takelog=True, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False, in_kpc_units=args.skip_re_scaling)
                            axes_err[index3, 0] = plot_2D_map(rescaled_direct_image, axes_err[index3, 0], f'{filter_for_direct_image}: {args.id}', args, cmap=direct_image_cmap, takelog=True, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=False, in_kpc_units=args.skip_re_scaling)

                        # ----------plotting the direct image: for debugging--------------
                        if args.debug_align:
                            re_pix = args.re_arcsec / args.pix_size_arcsec
                            fig_debug, axes_debug_2d = plt.subplots(1 + len(args.line_list), 5, figsize=(10, 7))
                            fig_debug.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)
                            fig_debug.text(0.05, 0.98, f'{args.field}: ID {args.id}: Alignment diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

                            axes = axes_debug_2d[0]
                            axes[0] = myimshow(direct_image, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=direct_image_cmap, col='w')
                            axes[1] = myimshow(rescaled_direct_image, axes[1], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=direct_image_cmap, col='w')
                            axes[2] = myimshow(recentered_direct_image, axes[2], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered', cmap=direct_image_cmap, col='w')

                        # ---------looping over all emission lines------------------------------
                        nlines_good = 0
                        inv_vars = []
                        for index4, this_line in enumerate(args.line_list):
                            if this_line not in args.available_lines:
                                print(f'\t\t\tLine {this_line} not available for object {args.id}. So skipping this line..')
                                inv_vars.append(np.zeros((args.npix_side, args.npix_side)))
                                if not args.debug_align:
                                    axes_flux[index3, index4 + 1].remove()
                                    axes_err[index3, index4 + 1].remove()
                                continue
                            try:
                                # -----------extracting the emission line map----------------
                                line_map, line_map_err, _ = get_emission_line_map(this_line, full_hdu, args, dered=False, silent=True)

                                # ----------flux scaling with complete error propagation----------
                                line_map_err_scaled = np.sqrt((line_map_err / integrated_scaling_flux)**2 + (line_map * scaling_flux_err / (integrated_scaling_flux**2))**2)
                                line_map_scaled = line_map / integrated_scaling_flux

                                # -----------combined single-pass transform (Rotation, Deprojection, Re/kpc Rescaling, and Center Offsets)--------
                                rescaled_line_map, rescaled_line_map_err = transform_map_and_err(line_map_scaled, line_map_err_scaled, args)

                                # ------recentering emission line map---------
                                #args.ndelta_xpix_linemap_rescaled, args.ndelta_ypix_linemap_rescaled = args.ndelta_xpix_direct_rescaled, args.ndelta_xpix_direct_rescaled
                                if index4 == 0:
                                    args.ndelta_xpix_linemap_rescaled, args.ndelta_ypix_linemap_rescaled = get_center_offsets(rescaled_line_map, args, silent=not args.debug_align)
                                recentered_line_map = ndimage.shift(rescaled_line_map, [args.ndelta_xpix_linemap_rescaled, args.ndelta_ypix_linemap_rescaled], order=0, cval=np.nan)
                                recentered_line_map_err = ndimage.shift(rescaled_line_map_err, [args.ndelta_xpix_linemap_rescaled, args.ndelta_ypix_linemap_rescaled], order=0, cval=np.nan)
                                
                                # -------plotting the intermediate steps: for debugging--------------
                                if args.debug_align:
                                    print(f'Deb345: shapes: {this_line} map = {np.shape(line_map)}, recaled = {np.shape(recentered_line_map)}, redshift={args.z:.1f}, log_mass={obj["log_mass"]:.1f}, log_sfr={obj["log_sfr"]:.1f}, a={args.semi_major:.1f}, b={args.semi_minor:.1f}, a/b={args.semi_major/args.semi_minor:.1f}, pa={args.pa:.1f}, re={re_pix:.1f} pixels') ##
                                    
                                    axes = axes_debug_2d[index4 + 1]
                                    axes[0] = myimshow(line_map, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label=f'{this_line}: Original', cmap=cmap, col='w' if args.fortalk else 'k', cmin=integrated_scaling_flux * 1e-3, cmax=integrated_scaling_flux * 0.1)
                                    axes[1] = myimshow(rescaled_line_map, axes[1], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=cmap, col='w', cmin=0, cmax=0.5)
                                    axes[2] = myimshow(recentered_line_map, axes[2], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered', cmap=cmap, col='w', cmin=0, cmax=0.5)
                                    axes[3] = myimshow(rescaled_line_map_err, axes[3], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled u', cmap=cmap, col='w', cmin=0, cmax=0.05)
                                    axes[4] = myimshow(recentered_line_map_err, axes[4], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered u', cmap=cmap, col='w', cmin=0, cmax=0.05)

                                # ---------------appending the line map-------------------------
                                line_maps_array[index4].append(recentered_line_map)
                                line_maps_err_array[index4].append(recentered_line_map_err)
                                constituent_ids_array[index4].append(f'{args.field}-{args.id}')

                                # -------computing the inv variance map------------
                                valid = ~np.isnan(recentered_line_map_err) & (recentered_line_map_err > 0)
                                inv_var = np.zeros_like(recentered_line_map_err)
                                inv_var[valid] = 1.0 / (recentered_line_map_err[valid] ** 2)
                                inv_vars.append(inv_var)

                                nlines_good += 1

                                # -------------plotting this line map of this galaxy-------------------
                                if not args.debug_align:
                                    axes_orig[index3, index4 + 1] = plot_2D_map(line_map, axes_orig[index3, index4 + 1], f'{this_line}: {args.id} OG', args, cmap=cmap, takelog=True, vmin=-20, vmax=-17, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True, segmentation_map=segmentation_map, in_re_units=False)
                                    axes_flux[index3, index4 + 1] = plot_2D_map(recentered_line_map, axes_flux[index3, index4 + 1], f'{this_line}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True, in_kpc_units=args.skip_re_scaling)
                                    axes_err[index3, index4 + 1] = plot_2D_map(recentered_line_map_err, axes_err[index3, index4 + 1], f'{this_line}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index3 < n_useful_rows_in_page - 1, hide_yaxis=True, in_kpc_units=args.skip_re_scaling)
                            
                            except Exception as e:
                                print(f'Skipping {this_line} for obj {args.id} because it failed due to: {e}')
                                continue
                            
                        print(f'\t\t\tFound {nlines_good} lines for {args.id}')
                        if nlines_good > 0: nobj_good += 1
                        
                        if args.debug_align:
                            figname = args.fig_dir / f'debug_align_{args.id}{args.deproject_text}{args.rescale_text}.png'
                            fig_debug.savefig(figname, transparent=args.fortalk)
                            print(f'\nSaved figure as {figname}')

                            plt.show(block=False)
                            sys.exit(f'Exiting here after {timedelta(seconds=(datetime.now() - start_time).seconds)} because of --debug_align mode; if you want to run the full code as usual then remove the --debug_align option and re-run')

                        # ----------compute master invariance map for this galaxy based on ALL line var maps--------                            
                        raw_w = np.nanmean(np.array(inv_vars), axis=0) # Mean inverse variance across all lines at each pixel
                        #raw_w = inv_vars[0] # OR, just consider the OIII inverse variance map
                        master_weight = np.zeros_like(raw_w)                      
                        valid_px = ~np.isnan(raw_w) & (raw_w > 0)
                        
                        # Mean-normalize spatial weights to 1 for democratic galaxy weighting
                        if np.any(valid_px):
                            mean_w = np.nanmean(raw_w[valid_px])
                            if mean_w > 0:
                                master_weight = raw_w / mean_w
                        master_weight_dict.update({f'{args.field}-{args.id}': master_weight}) # creating a galaxy-master weightmap pair

                        # -------cleaning------
                        del full_hdu
                        gc.collect()
                    
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
                stacked_maps, stacked_maps_err = [], []
                for index4, this_line in enumerate(args.line_list):
                    master_weight_arr =[master_weight_dict[this_id] for this_id in constituent_ids_array[index4]]
                    stacked_map, stacked_map_err = weighted_stack_line_maps(line_maps_array[index4], line_maps_err_array[index4], master_weight_arr)
                    stacked_maps.append(stacked_map)
                    stacked_maps_err.append(stacked_map_err)

                    # --------displaying stacked maps at the bottom of the mammoth figure---------
                    if not args.debug_align:
                        curr_row = (nrows_total % args.max_gal_per_page) - 1
                        if np.ndim(stacked_map) == 2 and np.ndim(stacked_map_err) == 2:
                            axes_orig[curr_row, index4 + 1] = plot_2D_map(stacked_map, axes_orig[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0, in_kpc_units=args.skip_re_scaling)
                            axes_flux[curr_row, index4 + 1] = plot_2D_map(stacked_map, axes_flux[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0, in_kpc_units=args.skip_re_scaling)
                            axes_err[curr_row, index4 + 1] = plot_2D_map(stacked_map_err, axes_err[curr_row, index4 + 1], f'{this_line}: Stacked', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index4 > 0, in_kpc_units=args.skip_re_scaling)
                        
                        # --------removing the first subplot (corresponding to direct image) at the bottom of the mammoth figure---------
                        axes_orig[curr_row, 0].set_visible(False)
                        axes_flux[curr_row, 0].set_visible(False)
                        axes_err[curr_row, 0].set_visible(False)

                # -------writing out stacked line maps as fits files--------------
                if args.do_all_obj: write_stacked_maps(stacked_maps, stacked_maps_err, constituent_ids_array, output_filename, args)
                
                # ------------finalising the mammoth PDFs------------
                if not args.debug_align:
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

        nobj_total_binned += nobj_good
        print(f'\nCompleted bin {bin_text} ({nobj_good} / {ngal_this_bin} objects, {nobj_no_scale_line} skipped due to lack of {args.scaling_line}) in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')
        if nobj_good > 1: nbin_good += 1

    print(f'\nBinned total ({nobj_total_binned} / {len(df)}) objects, into {len(bin_list)} bins.')
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')