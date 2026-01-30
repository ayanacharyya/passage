'''
    Filename: stack_mock_maps.py
    Notes: Stacks a given number of 2D mock galaxy maps; this is for testing and debugging if the stacking algorithm is working fine
    Author : Ayan
    Created: 29-01-26
    Example: run stack_mock_maps.py --field Par999 --do_all_obj --clobber --debug_bin --skip_deproject --skip_re_scaling
             run stack_mock_maps.py --field Par999 --do_all_obj --kpc_limit 5 --clobber_mock
             run stack_mock_maps.py --field Par999 --re_limit 2 --do_all_obj 
             run stack_mock_maps.py --field Par999 --re_limit 2 --id 99908 --debug_align 
'''

from header import *
from util import *
from make_diagnostic_maps import trim_image, myimshow, get_offsets_from_center
from stack_emission_maps import rotate_line_map, deproject_line_map, get_center_offsets, rescale_line_map, weighted_stack_line_maps, plot_2D_map, get_direct_image, write_stacked_maps, get_emission_line_map, setup_fullpage_figure

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_niriss_psf(filter_name, supersampling_factor=1, fov_pixels=41):
    '''
    Computes NIRISS PSF in a given filter
    Returns 2D array
    '''
    if fov_pixels % 2 == 0: fov_pixels += 1
    niriss = webbpsf.NIRISS()
    niriss.filter = filter_name    
    psf = niriss.calc_psf(fov_pixels=fov_pixels, oversample=supersampling_factor)
    psf_array  = psf[0].data
     
    return psf_array

# --------------------------------------------------------------------------------------------------------------------
def generate_mock_galaxy_fits(params, output_filename, line_name='X', filter_name='F150W', box_size=50, psf_kernel=None):
    '''
    Generates a mock 2D galaxy image based on input quantities like semi-major (a_image) and semi-minor (b_image) axis, PA (theta_image), sersic index (sersic) 
    and offset from center (offset_xpix, offet_ypix), SNR (snr) level
    Creates a (noiseless) direct image, a segmentation map, a (noisy) flux map, a flux uncertanity map, and stores all of them in a multi-extension hdulist
    Saves the hdulist as output_filename
    '''
    if psf_kernel is None:
        psf_kernel = get_niriss_psf(filter_name=params['filter'], fov_pixels=box_size)

    # --------computing sersic properties-----------
    ellip = 1 - (params['b_image'] / params['a_image']) # Calculate eccentricity and orientation; Sersic2D uses 'ellip' = 1 - (b/a) and theta in radians
    theta_rad = np.radians(params['theta_image'])
    
    re_arcsec = params['re_kpc'] * cosmo.arcsec_per_kpc_proper(params['redshift']).value
    re_pix =re_arcsec / params['pix_size_arcsec']

    model = Sersic2D(amplitude=1, r_eff=re_pix, n=params['sersic'],
                   x_0=box_size//2 + params['offset_xpix'], 
                   y_0=box_size//2 + params['offset_ypix'],
                   ellip=ellip, theta=theta_rad) # Generate the Noiseless Sersic Profile
    
    # --------making the direct image map-----------
    y, x = np.mgrid[0: box_size, 0: box_size]
    raw_img = model(x, y)
    img_convolved = convolve(raw_img, psf_kernel, boundary='fill', fill_value=0)
    img_noiseless = (img_convolved / np.sum(img_convolved)) * params['tot_bright'] # Scale to total brightness

    # --------making the segmentation map-----------
    seg_map = (img_noiseless > (0.01 * np.max(img_noiseless))).astype(int) * params['id']

    # --------computing background noise-----------
    peak_signal = np.max(img_noiseless)
    sigma_bg = peak_signal / params['snr']
    bg_noise = np.random.normal(0, sigma_bg, size=img_noiseless.shape)
    
    # --------computing Poisson noise-----------
    gain = params['gain']
    #gain = 1e5
    img_electrons = img_noiseless * gain
    poisson_noise_electrons = np.random.poisson(np.maximum(img_electrons, 0)) - img_electrons # Use max(0) to ensure no negative lambdas for Poisson
    poisson_noise_flux = poisson_noise_electrons / gain

    # --------making noisy flux map-----------
    img_noisy = img_noiseless + bg_noise + poisson_noise_flux
    
    # --------making the flux err map-----------
    var_bg = np.full(np.shape(img_noiseless), sigma_bg ** 2)
    variance_map = var_bg + (np.maximum(img_noiseless, 0) / gain)
    weight_map = 1.0 / (variance_map + 1e-20)

    # --------setting up the huglist-----------
    hdul = fits.HDUList([fits.PrimaryHDU()])
    
    hdr = hdul[0].header
    hdr['OBJ_ID'] = params['id']
    hdr['REDSHIFT'] = params['redshift']
    hdr['RE_KPC'] = params['re_kpc']
    hdr['SERSIC'] = params['sersic']
    hdr['SNR_IN'] = params['snr']
    hdr['HASLINES'] = line_name

    hdr2 = fits.Header()
    hdr2['PIXASEC'] = obj['pix_size_arcsec']
    hdr2['EXTVER'] = filter_name
    hdul.append(fits.ImageHDU(img_noiseless, name='DSCI', header=hdr2))

    hdul.append(fits.ImageHDU(seg_map, name='SEG'))

    hdr3 = fits.Header()
    hdr3['EXTVER'] = line_name
    hdr3['RESTWAVE'] = 9999 # dummy value in Angstrom
    hdul.append(fits.ImageHDU(img_noisy, name='LINE', header=hdr3))
    hdul.append(fits.ImageHDU(weight_map, name='LINEWHT', header=hdr3))

    # --------optionally, saving the fits file-----------
    hdul.writeto(output_filename, overwrite=True)
    print(f'Mock galaxy {params["field"]}-{params["id"]} saved to {output_filename}')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.5

    deproject_text = '_nodeproject' if args.skip_deproject else ''
    rescale_text = '_norescale' if args.skip_re_scaling else ''

    # -----------define colorbar properties-----------
    cmin, cmax, cmap = -18.5, -15.5, 'cividis'
    stacked_cmap = 'viridis'
    direct_image_cmap, direct_image_cmin, direct_image_cmax = 'Greys_r', 0, 1

    # -----------define global properties-----------
    line_name = 'XX'
    filter_name = 'F150W'
    pix_size_arcsec = 0.04
    mock_box_size = 51
    kernel_loaded = False
    args.arcsec_limit = mock_box_size * pix_size_arcsec / 2.
    args.line_list = [line_name]

    # -----------define object properties-----------
    data = {
        'id': [99901, 99902, 99903, 99904, 99905, 99906, 99907, 99908, 99909, 99910], # IDs
        'redshift': [2.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], # Redshift
        'pix_size_arcsec': [pix_size_arcsec] * 10, # fixed pixel scale
        'a_image':        [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], # Major axis
        'b_image':        [5.0, 5.0, 2.5, 1.5, 5.0, 5.0, 5.0, 5.0, 5.0, 1.5], # Minor axis
        'theta_image':    [0,   0,   0,   40,   0,   0,   0,   0,   0,   40], # Position Angle
        're_kpc':         [3.0, 3.0, 3.0, 3.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0], # Effective radii in kpc
        'offset_xpix':    [0,   0,   0,   0,   0,   -1,   0,  0,   0,    0 ], # Centering tests
        'offset_ypix':    [0,   0,   0,   0,   0,   2,   0,  0,   0,    0 ], # Centering tests
        'tot_bright':     [100, 100, 100, 100, 100,  100, 800, 100, 100, 100], # Brightness weighting
        'snr':            [10,  10,  10,  20,  10,   10,  10,  2,  10,  2], # Noise weighting
        'sersic':         [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0], # Concentration
        'gain':           [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0], # Gain parameter for making Poisson noise more/less agressive
        'field':          [args.field] * 10, # fixed field name
        'notes':          ['fiducial', 'high-z', 'ellipse', 'ell. + incl.', 'small', 'off-center', 'bright', 'high bkg noise', 'bulgy', 'ell. + incl. + noisy'] # to keep track of the purpose of each model
        }
    df_object = pd.DataFrame(data)
    if args.id is not None: df_object = df_object[df_object['id'].isin(args.id)].reset_index(drop=True)
    nrows_total = len(df_object) + 1
    
    # ----------define different iterations------------------
    skip_ids = [[0], [99902], [99903], [99904], [99905], [99906], [99905, 99907], [99908], [99909], [99910], [99908, 99904], [99910, 99909], [99910, 99903], [99901, 99902, 99905, 99906, 99907, 99908, 99909]] # each id here is skipped in each iteration; in the first iteration, no id is skipped
    n_iterations = len(skip_ids)

    # -------------for determining plot limits--------------------
    if args.skip_re_scaling: args.pix_size = (2 * args.kpc_limit) / args.npix_side # kpc
    else: args.pix_size = (2 * args.re_limit) / args.npix_side # Re

    offset = args.pix_size / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
    if args.skip_re_scaling: args.extent = (-args.kpc_limit - offset, args.kpc_limit - offset, -args.kpc_limit - offset, args.kpc_limit - offset)
    else: args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)

    # ------determining field-specific paths, etc-----------
    product_dir = args.input_dir / args.field / 'Products' / 'maps'
    product_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.output_dir / args.field / 'stacking_mock'
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / f'plots{deproject_text}{rescale_text}'
    fig_dir.mkdir(parents=True, exist_ok=True)
    Path(fig_dir / 'binmembers').mkdir(parents=True, exist_ok=True)
    fits_dir = output_dir / f'maps{deproject_text}{rescale_text}'
    fits_dir.mkdir(parents=True, exist_ok=True)
    output_filename = fits_dir / f'mock_stacked_maps{deproject_text}{rescale_text}.fits'

    # -------setting up massive PDF figures-----------
    if args.debug_align:
        fig_debug, axes_debug_2d = plt.subplots(2, 6, figsize=(13, 6))
        fig_debug.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)
    else:
        fullplot_filename_orig = fig_dir / f'binmembers/mock_binmembers_orig_maps{deproject_text}{rescale_text}.pdf'
        fullplot_filename_flux = fig_dir / f'binmembers/mock_binmembers_flux_maps{deproject_text}{rescale_text}.pdf'
        fullplot_filename_err = Path(str(fullplot_filename_flux).replace('flux', 'err'))
        
        pdf_orig = PdfPages(fullplot_filename_orig)
        pdf_flux = PdfPages(fullplot_filename_flux)
        pdf_err = PdfPages(fullplot_filename_err)

        args.max_gal_per_page = nrows_total
        fig_orig, axes_orig = setup_fullpage_figure(1, 1, n_iterations + 1,  cmin, cmax, cmap, args, in_kpc_units=args.skip_re_scaling)
        fig_flux, axes_flux = setup_fullpage_figure(1, 1, n_iterations + 1,  cmin, cmax, cmap, args, in_kpc_units=args.skip_re_scaling)
        fig_err, axes_err = setup_fullpage_figure(1, 1, n_iterations + 1,  cmin, cmax, cmap, args, in_kpc_units=args.skip_re_scaling)

    # --------looping over varios test iterations-------------
    for index, objects_to_skip in enumerate(skip_ids):
        print(f'Commencing iteration ({index + 1}/{len(skip_ids)}) to skip ID {objects_to_skip}..')

        line_maps_array = [[] for _ in range(len(args.line_list))] # each array must be len(args.line_list) long and each element will become a stack of 2D images
        line_maps_err_array = [[] for _ in range(len(args.line_list))] # each array must be len(args.line_list) long and each element will become a stack of 2D images
        constituent_ids_array = [[] for _ in range(len(args.line_list))]
        total_brightness_array = [[] for _ in range(len(args.line_list))]

        # ------------looping over objects for this iteration---------------
        for index2, obj in df_object.iterrows():
            args.id = obj['id']
            print(f'\tCommencing ({index2 + 1}/{len(df_object)}) IDs: {args.id}..')

            # ------skipping object based on the conditions of this iteration----------
            if args.do_all_obj and args.id in objects_to_skip:
                print(f'\t\tSkipping object {args.id}..')
                if not args.debug_align:
                    axes_orig[index2, index + 1].set_visible(False)
                    axes_flux[index2, index + 1].set_visible(False)
                    axes_err[index2, index + 1].set_visible(False)
                continue
                
            # ----------generating the mock object------------
            full_filename = product_dir / f'{obj["field"]}_{obj["id"]:05d}_maps.fits'
            if not os.path.exists(full_filename) or (args.clobber_mock and index == 0): # only clobber if this is the first iteration
                if not kernel_loaded:
                    psf_kernel = get_niriss_psf(filter_name=filter_name, fov_pixels=mock_box_size)
                    kernel_loaded = True
                generate_mock_galaxy_fits(obj, full_filename, line_name=line_name, filter_name=filter_name, box_size=mock_box_size, psf_kernel=psf_kernel)
            else:
                print(f'\t\tReading from existing {full_filename}')

            # ------------read in maps files--------------------------------
            full_hdu = fits.open(full_filename)
            
            # ----------determining object parameters------------
            args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
            args.z = full_hdu[0].header['REDSHIFT']
            args.distance = cosmo.luminosity_distance(args.z)
            args.pix_size_arcsec = full_hdu['DSCI', filter_name].header['PIXASEC']
            args.pix_size_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value
            args.pix_size_re = (2 * args.re_limit) / args.npix_side
            
            args.EB_V = 0.
            args.semi_major = obj['a_image']
            args.semi_minor = obj['b_image']
            args.pa = obj['theta_image']
            args.re_kpc = obj['re_kpc']
            args.re_arcsec = obj['re_kpc'] * cosmo.arcsec_per_kpc_proper(args.z).value
            
            # --------determining true center of object rom direct image---------------------
            args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter=filter_name, silent=not args.debug_align)

            # ---------------segmentation map---------------
            segmentation_map = full_hdu['SEG'].data
            segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
            segmentation_map = ndimage.shift(segmentation_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
            args.segmentation_map = segmentation_map
            
            rotated_segmentation_map = rotate_line_map(segmentation_map, args)
            deprojected_segmentation_map = rotated_segmentation_map if args.skip_deproject else deproject_line_map(rotated_segmentation_map, args)
            rescaled_segmentation_map = rescale_line_map(deprojected_segmentation_map, args)

            # ---------------direct image---------------
            direct_image, exptime = get_direct_image(full_hdu, filter_name, args) # this is already offset corrected and trimmed
            rotated_direct_image = rotate_line_map(direct_image, args)
            deprojected_direct_image = rotated_direct_image if args.skip_deproject else deproject_line_map(rotated_direct_image, args)
            rescaled_direct_image = rescale_line_map(deprojected_direct_image, args)
            total_brightness = np.nansum(direct_image)

            # ----------plotting the direct image: for debugging--------------
            if args.debug_align:
                re_pix = args.re_arcsec / args.pix_size_arcsec
                axes = axes_debug_2d[0]
                axes[0] = myimshow(direct_image, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=direct_image_cmap, col='w')
                axes[1].set_visible(False) # no smoothed direct image
                axes[2] = myimshow(rotated_direct_image, axes[2], contour=rotated_segmentation_map != args.id, re_pix=re_pix, label='Rotated', cmap=direct_image_cmap, col='w')
                axes[3] = myimshow(deprojected_direct_image, axes[3], contour=deprojected_segmentation_map != args.id, re_pix=re_pix, label='Deprojected', cmap=direct_image_cmap, col='w')
                axes[4] = myimshow(rescaled_direct_image, axes[4], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=direct_image_cmap, col='w')
            
            # ---------plotting direct image and the rescaled direct image------------------------------
            else:
                axes_orig[index2, 0] = plot_2D_map(direct_image, axes_orig[index2, 0], f'{obj["notes"]}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=False, segmentation_map=segmentation_map, in_re_units=False, seg_col='w')
                axes_flux[index2, 0] = plot_2D_map(rescaled_direct_image, axes_flux[index2, 0], f'{obj["notes"]}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=False, in_kpc_units=args.skip_re_scaling)
                axes_err[index2, 0] = plot_2D_map(rescaled_direct_image, axes_err[index2, 0], f'{obj["notes"]}', args, cmap=direct_image_cmap, takelog=False, vmin=direct_image_cmin, vmax=direct_image_cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=False, in_kpc_units=args.skip_re_scaling)

            # -----------extracting the emission line map----------------
            for index3, this_line in enumerate(args.line_list):
                line_map, line_map_err = get_emission_line_map(args.line_list[index3], full_hdu, args, dered=True, silent=True)
                total_brightness = np.nansum(line_map)

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
                deprojected_line_map = rotated_line_map if args.skip_deproject else deproject_line_map(rotated_line_map, args)
                deprojected_line_map_err = rotated_line_map_err if args.skip_deproject else deproject_line_map(rotated_line_map_err, args)

                # --------scaling the line map by effective radius---------------------
                rescaled_line_map = rescale_line_map(deprojected_line_map, args)
                rescaled_line_map_err = rescale_line_map(deprojected_line_map_err, args)

                # --------computing new recentering offsets, after rescaling---------------------
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
                    fig_debug.text(0.05, 0.98, f'{args.field}: ID {args.id}: {line_name} map: Alignment diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')
                    print(f'\t\tDeb325: shapes: {line_name} map = {np.shape(line_map)}, rotated = {np.shape(rotated_line_map)}, deprojected = {np.shape(deprojected_line_map)}, recaled = {np.shape(rescaled_line_map)}, redshift={args.z:.1f}, a={args.semi_major:.1f}, b={args.semi_minor:.1f}, a/b={args.semi_major/args.semi_minor:.1f}, pa={args.pa:.1f}, re={re_pix:.1f} pixels') ##
                    
                    axes = axes_debug_2d[1]
                    axes[0] = myimshow(line_map, axes[0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=cmap, col='k')
                    axes[1] = myimshow(smoothed_line_map, axes[1], contour=segmentation_map != args.id, re_pix=re_pix, label='Smoothed', cmap=cmap, col='k')
                    axes[2] = myimshow(rotated_line_map, axes[2], contour=rotated_segmentation_map != args.id, re_pix=re_pix, label='Rotated', cmap=cmap, col='k')
                    axes[3] = myimshow(deprojected_line_map, axes[3], contour=deprojected_segmentation_map != args.id, re_pix=re_pix, label='Deprojected', cmap=cmap, col='k')
                    axes[4] = myimshow(rescaled_line_map, axes[4], contour=rescaled_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Rescaled', cmap=cmap, col='k')
                    axes[5] = myimshow(recentered_line_map, axes[5], contour=recentered_segmentation_map != args.id, re_pix=args.npix_side / (2 * args.re_limit), label='Recentered', cmap=cmap, col='k')

                    figname = fig_dir / f'debug_align_{args.id}{deproject_text}{rescale_text}.png'
                    fig_debug.savefig(figname, transparent=args.fortalk)
                    print(f'\nSaved figure as {figname}')

                    plt.show(block=False)
                    sys.exit(f'\nExiting here because of --debug_align mode; if you want to run the full code as usual then remove the --debug_align option and re-run')

                # ---------------appending the line map-------------------------
                line_maps_array[index3].append(recentered_line_map)
                line_maps_err_array[index3].append(recentered_line_map_err)
                constituent_ids_array[index3].append(f'{args.field}-{args.id}')
                total_brightness_array[index3].append(total_brightness)
                
                # -------------plotting this line map of this galaxy-------------------
                if not args.debug_align:
                    axes_orig[index2, index + 1] = plot_2D_map(line_map, axes_orig[index2, index + 1], f'{line_name}: {args.id} OG', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=True, segmentation_map=segmentation_map, in_re_units=False)
                    axes_flux[index2, index + 1] = plot_2D_map(recentered_line_map, axes_flux[index2, index + 1], f'{line_name}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=True, in_kpc_units=args.skip_re_scaling)
                    axes_err[index2, index + 1] = plot_2D_map(recentered_line_map_err, axes_err[index2, index + 1], f'{line_name}: {args.id}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index2 < nrows_total - 1, hide_yaxis=True, in_kpc_units=args.skip_re_scaling)
        
        # -----stacking all line maps for all objects in this bin-----------
        stacked_maps, stacked_maps_err, nobj_arr = [], [], []
        for index3, this_line in enumerate(args.line_list):
            stacked_map, stacked_map_err = weighted_stack_line_maps(line_maps_array[index3], line_maps_err_array[index3], additional_weights= 1 / np.array(total_brightness_array[index3]) ** 2)
            stacked_maps.append(stacked_map)
            stacked_maps_err.append(stacked_map_err)
            nobj_arr.append(len(constituent_ids_array[index3]))

            # --------displaying stacked maps at the bottom of the mammoth figure---------
            if not args.debug_align:
                curr_row = - 1

                axes_orig[curr_row, 0].set_visible(False)
                axes_flux[curr_row, 0].set_visible(False)
                axes_err[curr_row, 0].set_visible(False)

                if np.ndim(stacked_map) == 2:
                    axes_orig[curr_row, index + 1] = plot_2D_map(stacked_map, axes_orig[curr_row, index + 1], f'{line_name}: Stacked', args, cmap=stacked_cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index > 0, in_kpc_units=args.skip_re_scaling)
                    axes_flux[curr_row, index + 1] = plot_2D_map(stacked_map, axes_flux[curr_row, index + 1], f'{line_name}: Stacked', args, cmap=stacked_cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index > 0, in_kpc_units=args.skip_re_scaling)
                    axes_err[curr_row, index + 1] = plot_2D_map(stacked_map_err, axes_err[curr_row, index + 1], f'{line_name}: Stacked', args, cmap=stacked_cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index > 0, in_kpc_units=args.skip_re_scaling)

        # -------writing out stacked line maps as fits files--------------
        if args.do_all_obj: write_stacked_maps(stacked_maps, stacked_maps_err, constituent_ids_array, output_filename, args)
    
    # -------------finalising saving each page of the mammoth PDFs (unless last page)-------------------------
    if not args.debug_align:
        pdf_orig.savefig(fig_orig)
        pdf_flux.savefig(fig_flux)
        pdf_err.savefig(fig_err)
        
        pdf_orig.close()
        pdf_flux.close()
        pdf_err.close()

        plt.close('all')                   
        print(f'\nSaved {fullplot_filename_flux}')

    print(f'\nCompleted in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
