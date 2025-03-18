'''
    Filename: extract_spectrum_from_grism.py
    Notes: Extracts 1D spectra and 2D emission line maps from calibrated grism data, following Vulcani+15 method
    Author : Ayan
    Created: 13-03-25
    Example: run extract_spectrum_from_grism.py --field Par028 --drv 0.5 --id 1303 --arcsec_limit 3 --plot_circle_at_radius 1
'''
from header import *
from util import *
from make_diagnostic_maps import annotate_PAs

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
def cutout_grism_spectra(ra, dec, filter, orient, data, wcs, pixscale, ax, args, cmap='Greys', cont=None):
    '''
    Cutout a rectangular region corresponding to the 2D spectra from the grism image and plot it on a given axis
    Returns axis handle and the 2D cutout spectra and corresponding 1D wavelength array
    '''
    # --------declaring filter-specific properties ---------
    # --from https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#gsc.tab=0-------
    df_filter_prop = pd.DataFrame({'filter':['f115w', 'f150w', 'f200w'],
                                   'lambda_cen':[1.148, 1.501, 1.989],
                                   'lambda_low':[1.013, 1.330, 1.751],
                                   'lambda_high':[1.283, 1.671, 2.226],
                                   'npix_extent':[60, 80, 120],
                                   'npix_offset_first':[-6, 55, 147],
                                   'npix_offset_zero': [-210, -220, -230]})

    # -----------------------------
    row = df_filter_prop[df_filter_prop['filter'] == filter]
    first_order_shift_pix = row['npix_offset_first'].values[0] # this should be a function of the filter
    lambda_extent_pix = row['npix_extent'].values[0] # this should be a function of the filter wavelength extent
    spatial_cutout_extent_pix = int((args.arcsec_limit/2 * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit
    spatial_extract_extent_pix = int((args.plot_circle_at_arcsec * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit

    pos_sky = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    pos_pix = PixCoord.from_sky(pos_sky, wcs)

    if orient == 'r':
        spectra_start_pix_y = pos_pix.y - first_order_shift_pix
        spectra_end_pix_y = spectra_start_pix_y - lambda_extent_pix
        spectra_center_pix_x = pos_pix.x
        spectra_center_pix_y = (spectra_start_pix_y + spectra_end_pix_y) / 2

        cutout = Cutout2D(data, [spectra_center_pix_x, spectra_center_pix_y], [lambda_extent_pix, spatial_cutout_extent_pix]) # cutout2D size = (ny, nx)

    elif orient == 'c':
        spectra_start_pix_x = pos_pix.x - first_order_shift_pix
        spectra_end_pix_x = spectra_start_pix_x - lambda_extent_pix
        spectra_center_pix_x = (spectra_start_pix_x + spectra_end_pix_x) / 2
        spectra_center_pix_y = pos_pix.y

        cutout = Cutout2D(data, [spectra_center_pix_x, spectra_center_pix_y], [spatial_cutout_extent_pix, lambda_extent_pix]) # cutout2D size = (ny, nx)

    cutout_data = cutout.data
    if cont is not None:
        try: cutout_data = cutout_data - cont
        except: cutout_data = (cutout_data.T - cont).T
    #cutout_data = np.log10(cutout_data)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout_data), cutout_data))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)

    if orient == 'r': cutout_data = cutout_data.T
    ax.imshow(cutout_data, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)

    ax.scatter(lambda_extent_pix / 2, spatial_cutout_extent_pix / 2, marker='x', s=10, c='r')
    ax.add_patch(plt.Rectangle((0, spatial_cutout_extent_pix/2 - spatial_extract_extent_pix/2), lambda_extent_pix, spatial_extract_extent_pix, color='r', fill=False, lw=0.5))
    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)#, bbox=dict(facecolor='k', edgecolor='black', alpha=0.5))

    ax = make_ax_labels(ax, '', '', args.fontsize,  hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

    #cutout_data = 10 ** cutout_data # to bring spec2d back to flux units
    cutout_data = cutout_data.T # this transpose is necessary to make dim_x = wavelength and dim_y = flux
    cutout_2d_spectra = cutout_data[0:lambda_extent_pix, int(spatial_cutout_extent_pix/2 - spatial_extract_extent_pix/2):int(spatial_cutout_extent_pix/2 + spatial_extract_extent_pix/2)]
    wave_arr = np.linspace(row['lambda_low'].values[0], row['lambda_high'].values[0], lambda_extent_pix)

    return ax, cutout_2d_spectra, wave_arr

# --------------------------------------------------------------------------------------------------------------------
def cutout_grism_image(ra, dec, filter, orient, data, wcs, pixscale, ax, args, cmap='Greys', hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Cutout a square region around (ra, dec) corresponding to direct image of an object from the grism image and plot it on a given axis
    Returns axis handle
    '''
    # --------declaring filter-specific properties ---------
    # --from https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr150-grisms#gsc.tab=0-------
    df_filter_prop = pd.DataFrame({'filter': ['f115w', 'f150w', 'f200w'],
                                   'lambda_cen': [1.148, 1.501, 1.989],
                                   'lambda_low': [1.013, 1.330, 1.751],
                                   'lambda_high': [1.283, 1.671, 2.226],
                                   'npix_extent': [60, 80, 120],
                                   'npix_offset_first': [-6, 55, 147],
                                   'npix_offset_zero': [214, 220, 230]})
    row = df_filter_prop[df_filter_prop['filter'] == filter]
    zero_order_shift_pix = row['npix_offset_zero'].values[0] # this should be a function of the filter

    # ------------determining ra and dec of zeroth order image------------------------
    size = args.arcsec_limit * 1 * u.arcsec
    pos_sky = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    pos_pix = PixCoord.from_sky(pos_sky, wcs)

    if orient == 'r':
        zero_order_shift_arcsec_x = (0 * u.arcsec)
        zero_order_shift_arcsec_y = (zero_order_shift_pix * u.pixel).to(u.arcsec, pixscale)
        zero_order_shift_pix_x = 0
        zero_order_shift_pix_y = zero_order_shift_pix
    elif orient == 'c':
        zero_order_shift_arcsec_x = (zero_order_shift_pix * u.pixel).to(u.arcsec, pixscale)
        zero_order_shift_arcsec_y = (0 * u.arcsec)
        zero_order_shift_pix_x = zero_order_shift_pix
        zero_order_shift_pix_y = 0

    pos_pix_zero = pos_pix + PixCoord(zero_order_shift_pix_x, zero_order_shift_pix_y)
    pos_sky_zero = pywcs.utils.pixel_to_skycoord(pos_pix_zero.x, pos_pix_zero.y, wcs)
    cutout = Cutout2D(data, pos_sky_zero, size, wcs=wcs)
    cutout_data = np.log10(cutout.data)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout_data), cutout_data))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)

    p = ax.imshow(cutout_data, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    ax.scatter(0, 0, marker='x', s=10, c='r')
    #ax.scatter((pos_sky_zero.data.lon.value - ra) * 3600, (pos_sky_zero.data.lat.value - dec) * 3600, marker='x', s=10, c='k')#, transform=ax.get_transform('fk5'))
    #ax.scatter(zero_order_shift_arcsec_x.value, zero_order_shift_arcsec_y.value, marker='x', s=20, c='pink')
    ax.text(0.98, 0.98, orient, fontsize=args.fontsize, c='r', ha='right', va='top', transform=ax.transAxes)#, bbox=dict(facecolor='k', edgecolor='black', alpha=0.5))

    if 'pa_arr' in args: ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize)
    if args.plot_circle_at_arcsec is not None: ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5))#, transform=ax.get_transform('fk5')))

    ax.coords['ra'].set_ticks(number=3)
    ax.coords['dec'].set_ticks(number=3)

    ax.coords['ra'].set_ticklabel(visible=False)
    ax.coords['dec'].set_ticklabel(visible=False)

    ax.coords.grid(color='k', alpha=0.5, linestyle='solid')

    # # -------just for testing----------
    # # -----------------------------
    # row = df_filter_prop[df_filter_prop['filter'] == filter]
    # first_order_shift_pix = row['npix_offset'].values[0] # this should be a function of the filter
    # lambda_extent_pix = row['npix_extent'].values[0] # this should be a function of the filter wavelength extent
    # spatial_cutout_extent_pix = int((args.arcsec_limit * u.arcsec).to(u.pixel, pixscale).value) # this should be a function of arcsec.limit
    #
    # pos_pix = PixCoord.from_sky(pos_sky, wcs)
    #
    # print(f'\nDeb159: {filter}-{orient}: object pos (sky)=({pos.data.lon.value:.4f}, {pos.data.lat.value:.4f}), \
    #         \nobject pos (pix)=({pos_pix.x}, {pos_pix.y}), \
    #         \noffset to 1st order (pix)={first_order_shift_pix}')
    #
    # if orient == 'r':
    #     spectra_start_pix_y = pos_pix.y - first_order_shift_pix
    #     spectra_end_pix_y = spectra_start_pix_y - lambda_extent_pix
    #     spectra_center_pix_x = pos_pix.x
    #     spectra_center_pix_y = (spectra_start_pix_y + spectra_end_pix_y) / 2
    #
    #     spectra_start_arcsec_y = ((spectra_start_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)
    #     spectra_end_arcsec_y = ((spectra_end_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)
    #     spectra_center_arcsec_x = ((spectra_center_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)
    #
    #     ax.plot([spectra_center_arcsec_x.value, spectra_center_arcsec_x.value], [spectra_start_arcsec_y.value, spectra_end_arcsec_y.value], lw=1, c='k')
    #
    # elif orient == 'c':
    #     spectra_start_pix_x = pos_pix.x - first_order_shift_pix
    #     spectra_end_pix_x = spectra_start_pix_x - lambda_extent_pix
    #     spectra_center_pix_x = (spectra_start_pix_x + spectra_end_pix_x) / 2
    #     spectra_center_pix_y = pos_pix.y
    #
    #     spectra_start_arcsec_x = ((spectra_start_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)
    #     spectra_end_arcsec_x = ((spectra_end_pix_x - pos_pix.x) * u.pixel).to(u.arcsec, pixscale)
    #     spectra_center_arcsec_y = ((spectra_center_pix_y - pos_pix.y) * u.pixel).to(u.arcsec, pixscale)
    #
    #     ax.plot([spectra_start_arcsec_x.value, spectra_end_arcsec_x.value], [spectra_center_arcsec_y.value, spectra_center_arcsec_y.value], lw=1, c='k')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def cutout_direct_image(ra, dec, data, wcs, ax, args, cmap='Greys', hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Cutout a square region around (ra, dec) from the direct/clear image and plot it on a given axis
    Returns axis handle
    '''
    size = args.arcsec_limit * u.arcsec
    pos = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    cutout = Cutout2D(data, pos, size, wcs=wcs).data
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
    Returns fit parameters
    '''
    contfit = np.polyfit(wave, spec, order)

    return contfit

# --------------------------------------------------------------------------------------------------------------------
def plot_cont(wave, spec2d, ax, redshift=0, fit_order=3, hide_xaxis=True, hide_yaxis=True, hide_cbar=True):
    '''
    Plots and fits the continuum corresponding to a grism image, on a given axis
    Returns axis handle and fit parameters
    '''
    spec1d = np.mean(spec2d, axis=1)
    wave = 1e4 * wave / (1 + redshift) # now in Angstrom
    ax.plot(wave, np.log10(spec1d), lw=1, c='salmon')

    contfit = fit_cont(wave, spec1d, order=fit_order)
    cont = np.poly1d(contfit)(wave)
    ax.plot(wave, np.log10(cont), lw=1, c='cornflowerblue')

    ax = make_ax_labels(ax, 'Restframe wavelength (A)', 'log flux', args.fontsize,  hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax, cont

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------setting up global values---------------------
    args.id_arr = np.atleast_1d(args.id)
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
    product_dir = args.input_dir / args.drv / f'{args.field}' / 'Products'
    catalog_file = args.output_dir / 'catalogs' / 'Par028_v0.5_venn_OIII,Ha,OII,Hb,SII,SNR>2.0,mass_df.txt'

    filter = 'f115w'
    line = 'Ha'
    dirimg_cmap, grism_dir_cmap, grism_2d_cmap = 'cividis', 'winter', 'winter'

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
        extract_dir = args.output_dir / f'{args.field}' / f'{args.id:05d}'
        extract_dir.mkdir(parents=True, exist_ok=True)

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

        #ra, dec, redshift = 150.0893979, 2.4278231, 0.001 # coords of a bright point source, for testing purposes
        # --------cutout and plot the direct image---------------
        ax_dirimg = cutout_direct_image(ra, dec, direct_data, direct_wcs, ax_dirimg, args, cmap=dirimg_cmap, hide_cbar=True)

        # ----------cutout and plot the grism images of r & c orients---------------
        #fig2 = plt.figure(figsize=(12, 6))
        #ax2 = fig2.add_subplot(111, projection=grism_c_wcs)
        #ax2 = cutout_grism_image(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax2, args, cmap=grism_dir_cmap, hide_xaxis=False, hide_yaxis=False, hide_cbar=False)
        ax_grism_dir_c = cutout_grism_image(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_dir_c, args, cmap=grism_dir_cmap, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)
        ax_grism_dir_r = cutout_grism_image(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_dir_r, args, cmap=grism_dir_cmap, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

        # ----------cutout and plot the grism 2D spectra of r & c orients---------------
        ax_grism_2d_c, spec2d_grism_c, wave_grism_c = cutout_grism_spectra(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_2d_c, args, cmap=grism_2d_cmap)
        ax_grism_2d_r, spec2d_grism_r, wave_grism_r = cutout_grism_spectra(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_2d_r, args, cmap=grism_2d_cmap)

        # --------model the 2D continuum and plot---------------
        ax_contfit_c, cont_grism_c = plot_cont(wave_grism_c, spec2d_grism_c, ax_contfit_c, redshift=redshift, fit_order=3, hide_xaxis=True, hide_yaxis=False, hide_cbar=True)
        ax_contfit_r, cont_grism_r = plot_cont(wave_grism_r, spec2d_grism_r, ax_contfit_r, redshift=redshift, fit_order=3, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)

        # --------plot the continuum subtracted 2D spectra of r & c orients---------------
        ax_grism_sub_c, _, _ = cutout_grism_spectra(ra, dec, filter, 'c', grism_c_data, grism_c_wcs, grism_c_pixscale, ax_grism_sub_c, args, cmap=grism_2d_cmap, cont=cont_grism_c)
        ax_grism_sub_r, _, _ = cutout_grism_spectra(ra, dec, filter, 'r', grism_r_data, grism_r_wcs, grism_r_pixscale, ax_grism_sub_r, args, cmap=grism_2d_cmap, cont=cont_grism_r)

        # --------plot the emission maps of r & c orients---------------

        # --------rotate and plot the combined the emission map---------------

        # ------------saving the full figure--------------------------
        figname = extract_dir / f'{args.id:05d}_extraction.png'
        fig.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure as {figname}')
        plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
