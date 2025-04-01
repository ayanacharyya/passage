'''
    Filename: extract_spectrum_from_grism.py
    Notes: Extracts 1D spectra and 2D emission line maps from calibrated grism data, following Vulcani+15 method
    Author : Ayan
    Created: 13-03-25
    Example: run extract_spectrum_from_grism.py --field Par028 --drv 0.5 --id 1303,1934,2867,300 --arcsec_limit 3 --plot_circle_at_radius 1 --test_cutout --debug_zero_order
'''
from header import *
from util import *
from make_diagnostic_maps import annotate_PAs
from scipy.interpolate import UnivariateSpline

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_fits_data(filename, ext=0):
    '''
    Opens a fits file and returns the data in the given extension as well as the corresponding header as WCS
    '''
    hdu = fits.open(filename)
    data = hdu[ext].data
    header = hdu[ext].header
    wcs = pywcs.WCS(header)
    pixscale = header['PIXSCALE'] if 'PIXSCALE' in header else np.nan
    pixscale = u.pixel_scale(pixscale * u.arcsec/u.pixel)

    return data, pixscale, wcs

# --------------------------------------------------------------------------------------------------------------------
def make_ax_labels(ax, xlabel, ylabel, fontsize, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, p=None):
    '''
    Puts the x and y labels on a given axis handle for a given fontsize, modifies the tick label sizes to the same size
    Returns axis handle
    '''
    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.tick_params(axis='x', which='major', labelsize=fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)

    if not hide_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def cutout_grism_spectra(ra, dec, filter, orient, data, wcs, pixscale, ax, args, cmap='Greys', cont=None, show_log=False, mask=None):
    '''
    Cutout a rectangular region corresponding to the 2D spectra from the grism image and plot it on a given axis
    Returns axis handle and the 2D cutout spectra and corresponding 1D wavelength array
    '''
    # -----------------------------
    row = df_filter_prop[df_filter_prop['filter'] == filter]
    first_order_shift_pix = row['npix_offset_first'].values[0] # this should be a function of the filter
    if orient == 'r': first_order_shift_pix+= 6
    lambda_extent_pix = row['npix_extent'].values[0] # this should be a function of the filter wavelength extent
    spatial_cutout_extent_pix = int((args.arcsec_limit/2 * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit
    spatial_extract_extent_pix = int((args.extract_arcsec * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit
    spatial_extract_offset_pix = row[f'spatial_offset_{orient}'].values[0] # to correct for a weird offset in the grism images

    pos_sky = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    pos_pix = PixCoord.from_sky(pos_sky, wcs)

    if orient == 'r':
        spectra_start_pix_y = pos_pix.y - first_order_shift_pix
        spectra_end_pix_y = spectra_start_pix_y - lambda_extent_pix
        spectra_center_pix_x = pos_pix.x - spatial_extract_offset_pix
        spectra_center_pix_y = (spectra_start_pix_y + spectra_end_pix_y) / 2

        cutout = Cutout2D(data, [spectra_center_pix_x, spectra_center_pix_y], [lambda_extent_pix, spatial_cutout_extent_pix]) # cutout2D size = (ny, nx)

    elif orient == 'c':
        spectra_start_pix_x = pos_pix.x - first_order_shift_pix
        spectra_end_pix_x = spectra_start_pix_x - lambda_extent_pix
        spectra_center_pix_x = (spectra_start_pix_x + spectra_end_pix_x) / 2
        spectra_center_pix_y = pos_pix.y - spatial_extract_offset_pix

        cutout = Cutout2D(data, [spectra_center_pix_x, spectra_center_pix_y], [spatial_cutout_extent_pix, lambda_extent_pix]) # cutout2D size = (ny, nx)

    cutout_data = cutout.data
    if cont is not None:
        try: cutout_data = cutout_data - cont
        except: cutout_data = (cutout_data.T - cont).T
    if show_log: cutout_data = np.log10(cutout_data)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout_data), cutout_data))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)

    if orient == 'r': cutout_data = cutout_data.T
    #p = ax.imshow(cutout_data, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    p = ax.imshow(cutout_data, cmap=cmap, origin='lower', vmin=0 if not show_log else vmin, vmax=3 if args.test_cutout and not show_log else 0.07 if not show_log else vmax)

    ax.scatter(lambda_extent_pix / 2, spatial_cutout_extent_pix / 2, marker='x', s=10, c='r')
    ax.add_patch(plt.Rectangle((0, spatial_cutout_extent_pix/2 - spatial_extract_extent_pix/2), lambda_extent_pix, spatial_extract_extent_pix, color='r', fill=False, lw=0.5))
    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)#, bbox=dict(facecolor='k', edgecolor='black', alpha=0.5))

    ax = make_ax_labels(ax, '', '', args.fontsize, p=p, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

    if show_log: cutout_data = 10 ** cutout_data # to bring spec2d back to flux units
    cutout_data = cutout_data.T # this transpose is necessary to make dim_x = wavelength and dim_y = flux
    cutout_2d_spectra = cutout_data[:, int(spatial_cutout_extent_pix/2 - spatial_extract_extent_pix/2):int(spatial_cutout_extent_pix/2 + spatial_extract_extent_pix/2)]
    wave_arr = np.linspace(row['lambda_low'].values[0], row['lambda_high'].values[0], lambda_extent_pix)[:np.shape(cutout_2d_spectra)[0]] # in microns, observed frame

    if mask is not None:
        mask_wave_pix = [np.where(wave_arr >= item)[0][0] for item in mask]
        ax.add_patch(plt.Rectangle((mask_wave_pix[0], spatial_cutout_extent_pix / 2 - spatial_extract_extent_pix / 2), mask_wave_pix[1] - mask_wave_pix[0], spatial_extract_extent_pix, color='k', fill=False, lw=0.5))

    return ax, cutout_2d_spectra, wave_arr

# --------------------------------------------------------------------------------------------------------------------
def cutout_grism_image(ra, dec, filter, orient, data, wcs, pixscale, ax, args, cmap='Greys', hide_xaxis=False, hide_yaxis=False, hide_cbar=False, debug_zero_order=False):
    '''
    Cutout a square region around (ra, dec) corresponding to direct image of an object from the grism image and plot it on a given axis
    Returns axis handle
    '''
    row = df_filter_prop[df_filter_prop['filter'] == filter]
    zero_order_shift_pix = row['npix_offset_zero'].values[0] # this should be a function of the filter
    spatial_extract_offset_pix = row[f'spatial_offset_{orient}'].values[0] # to correct for a weird offset in the grism images

    # ------------determining ra and dec of zeroth order image------------------------
    if debug_zero_order: size = args.arcsec_limit * 10 * u.arcsec
    else: size = args.arcsec_limit * 1 * u.arcsec
    pos_sky = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    pos_pix = PixCoord.from_sky(pos_sky, wcs)

    if orient == 'r':
        zero_order_shift_pix_x = -spatial_extract_offset_pix
        zero_order_shift_pix_y = zero_order_shift_pix
    elif orient == 'c':
        zero_order_shift_pix_x = zero_order_shift_pix
        zero_order_shift_pix_y = -spatial_extract_offset_pix

    zero_order_shift_arcsec_x = (zero_order_shift_pix_x * u.pixel).to(u.arcsec, pixscale)
    zero_order_shift_arcsec_y = (zero_order_shift_pix_y * u.pixel).to(u.arcsec, pixscale)

    pos_pix_zero = pos_pix + PixCoord(zero_order_shift_pix_x, zero_order_shift_pix_y)
    pos_sky_zero = pywcs.utils.pixel_to_skycoord(pos_pix_zero.x, pos_pix_zero.y, wcs)
    if debug_zero_order: cutout = Cutout2D(data, pos_sky, size, wcs=wcs)
    else: cutout = Cutout2D(data, pos_sky_zero, size, wcs=wcs)
    cutout_data = np.log10(cutout.data)

    try:
        good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout_data), cutout_data))
        vmin = np.percentile(good_data, 1)
        vmax = np.percentile(good_data, 99)
    except IndexError:
        vmin, vmax = None, None

    p = ax.imshow(cutout_data, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax, extent=args.extent)
    ax.scatter(0, 0, marker='x', s=10, c='r')
    # ax.scatter((pos_sky_zero.data.lon.value - ra) * 3600, (pos_sky_zero.data.lat.value - dec) * 3600, marker='x', s=10, c='k')#, transform=ax.get_transform('fk5'))
    if debug_zero_order: ax.scatter(zero_order_shift_arcsec_x.value, zero_order_shift_arcsec_y.value, marker='x', s=20, c='pink')
    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)

    if 'pa_arr' in args: ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize)
    if args.plot_circle_at_arcsec is not None: ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5))#, transform=ax.get_transform('fk5')))

    ax.coords['ra'].set_ticks(number=3)
    ax.coords['dec'].set_ticks(number=3)

    ax.coords['ra'].set_ticklabel(visible=False)
    ax.coords['dec'].set_ticklabel(visible=False)

    ax.coords.grid(color='k', alpha=0.5, linestyle='solid')

    # -------just for testing----------
    if False: #debug_zero_order:
        # -----------------------------
        row = df_filter_prop[df_filter_prop['filter'] == filter]
        first_order_shift_pix = row['npix_offset_first'].values[0] # this should be a function of the filter
        lambda_extent_pix = row['npix_extent'].values[0] # this should be a function of the filter wavelength extent
        spatial_cutout_extent_pix = int((args.arcsec_limit * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit

        pos_pix = PixCoord.from_sky(pos_sky, wcs)

        print(f'\nDeb159: {filter}-{orient}: object pos (sky)=({pos_sky.data.lon.value:.4f}, {pos_sky.data.lat.value:.4f}), \
                \nobject pos (pix)=({pos_pix.x}, {pos_pix.y}), \
                \noffset to 1st order (pix)={first_order_shift_pix}')

        if orient == 'r':
            spectra_start_pix_y = pos_pix.y - first_order_shift_pix
            spectra_end_pix_y = spectra_start_pix_y - lambda_extent_pix
            spectra_center_pix_x = pos_pix.x - spatial_extract_offset_pix
            spectra_center_pix_y = (spectra_start_pix_y + spectra_end_pix_y) / 2

            spectra_start_arcsec_y = ((spectra_start_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)
            spectra_end_arcsec_y = ((spectra_end_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)
            spectra_center_arcsec_x = ((spectra_center_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)

            ax.plot([spectra_center_arcsec_x.value, spectra_center_arcsec_x.value], [spectra_start_arcsec_y.value, spectra_end_arcsec_y.value], lw=1, c='k')

        elif orient == 'c':
            spectra_start_pix_x = pos_pix.x - first_order_shift_pix
            spectra_end_pix_x = spectra_start_pix_x - lambda_extent_pix
            spectra_center_pix_x = (spectra_start_pix_x + spectra_end_pix_x) / 2
            spectra_center_pix_y = pos_pix.y - spatial_extract_offset_pix

            spectra_start_arcsec_x = ((spectra_start_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)
            spectra_end_arcsec_x = ((spectra_end_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)
            spectra_center_arcsec_y = ((spectra_center_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)

            ax.plot([spectra_start_arcsec_x.value, spectra_end_arcsec_x.value], [spectra_center_arcsec_y.value, spectra_center_arcsec_y.value], lw=1, c='k')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def cutout_direct_image(ra, dec, data, wcs, ax, args, cmap='Greys', hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Cutout a square region around (ra, dec) from the direct/clear image and plot it on a given axis
    Returns axis handle
    '''
    size = args.arcsec_limit * u.arcsec
    pos_sky = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    cutout = Cutout2D(data, pos_sky, size, wcs=wcs).data
    cutout = np.log10(cutout)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout), cutout))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)

    p = ax.imshow(cutout, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    ax.scatter(0, 0, marker='x', s=10, c='r')

    if 'pa_arr' in args: ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize)
    if args.plot_circle_at_arcsec is not None: ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5))

    ax.coords['ra'].set_major_formatter('d.ddddd')
    ax.coords['dec'].set_major_formatter('d.ddddd')
    ax.coords['ra'].set_axislabel('RA', fontsize=args.fontsize)
    ax.coords['dec'].set_axislabel('Dec', fontsize=args.fontsize)

    ax.coords['ra'].set_ticks(number=3)
    ax.coords['dec'].set_ticks(number=3)

    ax.coords['ra'].set_ticklabel(size=args.fontsize)
    ax.coords['dec'].set_ticklabel(size=args.fontsize, rotation=90)

    ax.coords.grid(color='k', alpha=0.5, linestyle='solid')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def fit_cont(wave, spec, order=3):
    '''
    Fits the given continuum using spline
    Returns fitted continuum
    '''
    #contfit = np.ma.polyfit(wave, spec, order)
    #cont = np.poly1d(contfit)(wave)

    contfit = UnivariateSpline(wave, spec, k=order)
    cont = contfit(wave)

    return cont

# --------------------------------------------------------------------------------------------------------------------
def plot_cont(wave, spec2d, ax, masks, fit_order=3, orient='c', hide_xaxis=True, hide_yaxis=True, hide_cbar=True):
    '''
    Plots and fits the continuum corresponding to a grism image, on a given axis
    Returns axis handle and fit parameters
    '''
    spec1d = np.mean(spec2d, axis=1)
    ax.plot(wave, spec1d, lw=1, c='salmon')

    for mask in masks: # masking out intended line/s
        spec1d = np.ma.masked_where((wave >= mask[0]) & (wave <= mask[1]), spec1d)
        ax.fill_betweenx([ax.get_ylim()[0], ax.get_ylim()[1]], mask[0], mask[1], color='grey', alpha=0.3, lw=0)

    cont = fit_cont(wave, spec1d, order=fit_order)
    ax.plot(wave, cont, lw=1, c='cornflowerblue')
    ax.plot(wave, spec1d.data - cont, lw=1, c='darkseagreen')

    ax.axhline(0, c='k', lw=0.5, ls='dashed')
    ax.set_ylim(-0.01, 1 if args.test_cutout else 0.07)
    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)
    ax = make_ax_labels(ax, 'Observed wavelength (microns)', 'Flux', args.fontsize,  hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax, cont


# --------------------------------------------------------------------------------------------------------------------
def plot_emission_map_rc(wave, spec2d, ax, mask, orient='c', cmap='viridis', hide_xaxis=True, hide_yaxis=True, hide_cbar=True, show_log=True):
    '''
    Cutout a region of dimensions wavelength mask x extract_arcsec from the continuum subtracted 2D grism spectra and plot it on a given axis
    Returns axis handle and 2D emission map
    '''
    good_wave_pix = (wave >= mask[0]) & (wave <= mask[1])
    line_map = spec2d[good_wave_pix, :]
    line_map = line_map.T
    x_shape, y_shape = np.shape(line_map)
    line_map = rebin(line_map, (min(x_shape, y_shape), min(x_shape, y_shape)))

    if show_log: line_map = np.log10(line_map)
    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(line_map), line_map))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)
    p = ax.imshow(line_map, cmap=cmap, origin='lower', vmin=0, vmax=3 if args.test_cutout else 0.07)

    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)
    ax = make_ax_labels(ax, '', '', args.fontsize, p=p, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    if show_log: line_map = 10 ** line_map

    return ax, line_map

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_map_combined(line_map_c, line_map_r, ax, args, cmap='viridis', hide_xaxis=True, hide_yaxis=True, hide_cbar=True, show_log=True):
    '''
    Combine 2D line maps from 2 grism orients and and plot it on a given axis
    Returns axis handle and the combined 2D map
    '''

    line_map = np.mean([line_map_c, line_map_r], axis=0)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(line_map), line_map))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)
    p = ax.imshow(line_map, cmap=cmap, origin='lower',  vmin=0, vmax=3 if args.test_cutout else 0.07)

    ax.text(0.98, 0.98, 'comb', fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)
    ax = make_ax_labels(ax, '', '', args.fontsize, p=p, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax, line_map

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------setting up hardcoded values---------------------
    filter = 'f200w'
    filter_line_dict = {'f115w':[['OII', 1.20, 1.24]], 'f150w':[['OIII', 1.70, 1.74]], 'f200w':[['Ha', 2.08, 2.12]]}
    dirimg_cmap, grism_dir_cmap, grism_2d_cmap = 'cividis', 'winter', 'winter'
    contfit_order = 5

    # --------declaring filter-specific properties ---------
    # --from https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#gsc.tab=0-------
    df_filter_prop = pd.DataFrame({'filter': ['f115w', 'f150w', 'f200w'],
                                   'lambda_cen': [1.148, 1.501, 1.989],
                                   'lambda_low': [1.013, 1.330, 1.751],
                                   'lambda_high': [1.283, 1.671, 2.226],
                                   'npix_extent': [60, 80, 120],
                                   'npix_offset_first': [-6, 55, 147],
                                   'npix_offset_zero': [214, 220, 206],
                                   'spatial_offset_c': [2, 2, 1],
                                   'spatial_offset_r': [2, 2, 2]})

    # --------setting up global values---------------------
    args.id_arr = np.atleast_1d(args.id)
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
    product_dir = args.input_dir / f'{args.field}' / 'Products'
    catalog_file = args.output_dir / 'catalogs' / 'Par028_v0.5_venn_OIII,Ha,OII,Hb,SII,SNR>2.0_df.txt'

    args.lines = [item[0] for item in filter_line_dict[filter]]
    args.line_masks = [[item[1], item[2]] for item in filter_line_dict[filter]]

    if catalog_file.exists(): df_cat = pd.read_csv(catalog_file)
    else: df_cat = None

    # --------------read in the full direct image file-----------------------
    direct_image_file = product_dir / f'{args.field}_{filter}-clear_drz_sci.fits'
    direct_data, direct_pixscale, direct_wcs = get_fits_data(direct_image_file)

    # --------------read in the full grism files-----------------------
    grism_c_file = product_dir / f'{args.field}_{filter}-gr150c_drz_sci.fits'
    grism_r_file = product_dir / f'{args.field}_{filter}-gr150r_drz_sci.fits'
    grism_c_data, grism_c_pixscale, grism_c_wcs = get_fits_data(grism_c_file)
    grism_r_data, grism_r_pixscale, grism_r_wcs = get_fits_data(grism_r_file)

    # --------looping over objects---------------------
    for index, args.id in enumerate(args.id_arr):
        print(f'Doing object {args.id} which is {index + 1} out of {len(args.id_arr)} objects..')

        # ---------setting up filenames-----------------------------
        plots_dir = args.output_dir / f'{args.field}' / f'non_grizli_extractions'
        plots_dir.mkdir(parents=True, exist_ok=True)

        # --------setting up full figure----------------------------------
        nrow, ncol = 4, 9
        fig = plt.figure(figsize=(12, 6), layout='constrained') # layout = 'tight' or 'constrained'

        ax_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=2, rowspan=2, projection=direct_wcs)
        ax_grism_dir_c = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 2), colspan=1, rowspan=1, projection=grism_c_wcs)
        ax_grism_dir_r = plt.subplot2grid(shape=(nrow, ncol), loc=(1, 2), colspan=1, rowspan=1, projection=grism_r_wcs)
        ax_grism_2d_c = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 3), colspan=3, rowspan=1)
        ax_grism_2d_r = plt.subplot2grid(shape=(nrow, ncol), loc=(1, 3), colspan=3, rowspan=1)

        ax_contfit_c = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 6), colspan=3, rowspan=1)
        ax_contfit_r = plt.subplot2grid(shape=(nrow, ncol), loc=(1, 6), colspan=3, rowspan=1)

        ax_grism_sub_c = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 3), colspan=3, rowspan=1)
        ax_grism_sub_r = plt.subplot2grid(shape=(nrow, ncol), loc=(3, 3), colspan=3, rowspan=1)
        ax_em_c = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 6), colspan=1, rowspan=1)#, projection=grism_c_wcs)
        ax_em_r = plt.subplot2grid(shape=(nrow, ncol), loc=(3, 6), colspan=1, rowspan=1)#, projection=grism_r_wcs)
        ax_em = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 7), colspan=2, rowspan=2)

        # -------get object properites from photcat---------------------
        if df_cat is not None:
            ra = df_cat[df_cat['objid'] == args.id]['ra'].values[0]
            dec = df_cat[df_cat['objid'] == args.id]['dec'].values[0]
            redshift = df_cat[df_cat['objid'] == args.id]['redshift'].values[0]
        else:
            sys.exit(f'Photcat file does not exist in {phot_catalog_file}. So where do I get the RA/DEC from?')

        if args.test_cutout: ra, dec, redshift = 150.0893979, 2.4278231, 0.001 # coords of a bright point source, for testing purposes
        # --------cutout and plot the direct image---------------
        ax_dirimg = cutout_direct_image(ra, dec, direct_data, direct_wcs, ax_dirimg, args, cmap=dirimg_cmap, hide_cbar=True)

        # ----------cutout and plot the grism images of r & c orients---------------
        if args.debug_zero_order:
            fig2 = plt.figure(figsize=(12, 6))
            ax2 = fig2.add_subplot(111, projection=grism_c_wcs)
            ax2 = cutout_grism_image(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax2, args, cmap=grism_dir_cmap, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, debug_zero_order=args.debug_zero_order)
        ax_grism_dir_c = cutout_grism_image(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_dir_c, args, cmap=grism_dir_cmap, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)
        ax_grism_dir_r = cutout_grism_image(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_dir_r, args, cmap=grism_dir_cmap, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

        # ----------cutout and plot the grism 2D spectra of r & c orients---------------
        ax_grism_2d_c, spec2d_grism_c, wave_grism_c = cutout_grism_spectra(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_2d_c, args, cmap=grism_2d_cmap, show_log=False)
        ax_grism_2d_r, spec2d_grism_r, wave_grism_r = cutout_grism_spectra(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_2d_r, args, cmap=grism_2d_cmap, show_log=False)

        # --------model the 2D continuum and plot---------------
        ax_contfit_c, cont_grism_c = plot_cont(wave_grism_c, spec2d_grism_c, ax_contfit_c, args.line_masks, fit_order=contfit_order, orient='c', hide_xaxis=True, hide_yaxis=False, hide_cbar=True)
        ax_contfit_r, cont_grism_r = plot_cont(wave_grism_r, spec2d_grism_r, ax_contfit_r, args.line_masks, fit_order=contfit_order, orient='r', hide_xaxis=False, hide_yaxis=False, hide_cbar=True)

        # --------plot the continuum subtracted 2D spectra of r & c orients---------------
        ax_grism_sub_c, spec2d_contsub_grism_c, _ = cutout_grism_spectra(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_sub_c, args, cmap=grism_2d_cmap, show_log=False, cont=cont_grism_c, mask=args.line_masks[0])
        ax_grism_sub_r, spec2d_contsub_grism_r, _ = cutout_grism_spectra(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_sub_r, args, cmap=grism_2d_cmap, show_log=False, cont=cont_grism_r, mask=args.line_masks[0])

        # --------plot the emission maps of r & c orients---------------
        ax_em_c, line_map_c = plot_emission_map_rc(wave_grism_c, spec2d_contsub_grism_c, ax_em_c, args.line_masks[0], cmap=grism_dir_cmap, orient='c', hide_xaxis=True, hide_yaxis=True, hide_cbar=True, show_log=False)
        ax_em_r, line_map_r = plot_emission_map_rc(wave_grism_r, spec2d_contsub_grism_r, ax_em_r, args.line_masks[0], cmap=grism_dir_cmap, orient='r', hide_xaxis=True, hide_yaxis=True, hide_cbar=True, show_log=False)

        # --------rotate and plot the combined the emission map---------------
        ax_em, line_map = plot_emission_map_combined(line_map_c, line_map_r, ax_em, args, cmap=grism_dir_cmap, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

        # ------------saving the full figure--------------------------
        figname = plots_dir / f'{args.id:05d}_extraction.png'
        fig.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure as {figname}')
        plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
