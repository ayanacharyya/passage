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
    wcs = pywcs.WCS(hdu[ext].header)

    return data, wcs

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
def cutout2d_and_plot(ra, dec, direct_data, direct_wcs, ax, args, cmap='Greys', hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Cutout a rectangular region around (x,y) from a bigger image (direct_data) and plot it on a given axis
    Returns axis handle
    '''
    size = args.arcsec_limit * u.arcsec  # 3 arcsec so that it leads to shape 125 x 125 pixels
    pos = SkyCoord(*(np.array([ra, dec]) * u.deg), frame='fk5')
    cutout = Cutout2D(direct_data, pos, size, wcs=direct_wcs).data
    cutout = np.log10(cutout)

    good_data = np.ma.compressed(np.ma.masked_where(~np.isfinite(cutout), cutout))
    vmin = np.percentile(good_data, 1)
    vmax = np.percentile(good_data, 99)

    p = ax.imshow(cutout, cmap=cmap, origin='lower', extent=args.extent, alpha=1, vmin=vmin, vmax=vmax)
    ax.scatter(0, 0, marker='x', s=10, c='grey')

    if 'pa_arr' in args: ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize)
    if args.plot_circle_at_arcsec is not None: ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5))

    ax = make_ax_labels(ax, 'RA (")', 'Dec (")', args.fontsize, p=p, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------setting up global values---------------------
    args.id_arr = np.atleast_1d(args.id)
    args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
    product_dir = args.input_dir / args.drv / f'{args.field}' / 'Products'
    phot_catalog_file = product_dir / f'{args.field}_photcat.fits'

    filter = 'f200w'
    line = 'Ha'
    dirimg_cmap = 'cividis'

    if phot_catalog_file.exists(): df_photcat = Table.read(phot_catalog_file).to_pandas()
    else: df_photcat = None

    # --------------read in the full direct image file-----------------------
    direct_image_file = product_dir / f'{args.field}_{filter}-clear_drz_sci.fits'
    direct_data, direct_wcs = get_fits_data(direct_image_file)

    # --------------read in the full grism files-----------------------
    grism_c_file = product_dir / f'{args.field}_{filter}-gr150c_drz_sci.fits'
    grism_r_file = product_dir / f'{args.field}_{filter}-gr150r_drz_sci.fits'
    grism_c_data, grism_c_wcs = get_fits_data(grism_c_file)
    grism_r_data, grism_r_wcs = get_fits_data(grism_r_file)

    # --------looping over objects---------------------
    for index, args.id in enumerate(args.id_arr):
        print(f'Doing object {args.id} which is {index + 1} out of {len(args.id_arr)} objects..')

        # ---------setting up filenames-----------------------------
        extract_dir = args.output_dir / f'{args.field}' / f'{args.id:05d}'
        extract_dir.mkdir(parents=True, exist_ok=True)

        # --------setting up full figure----------------------------------
        nrow, ncol = 4, 8
        fig = plt.figure(figsize=(12, 6), layout='constrained') # layout = 'tight' or 'constrained'

        ax_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=2, rowspan=2)
        ax_grism_raw_c = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 2), colspan=3, rowspan=1)
        ax_grism_raw_r = plt.subplot2grid(shape=(nrow, ncol), loc=(1, 2), colspan=3, rowspan=1)
        ax_contfit = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 5), colspan=3, rowspan=2)

        ax_grism_sub_c = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 2), colspan=3, rowspan=1)
        ax_grism_sub_r = plt.subplot2grid(shape=(nrow, ncol), loc=(3, 2), colspan=3, rowspan=1)
        ax_em_c = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 5), colspan=1, rowspan=1)
        ax_em_r = plt.subplot2grid(shape=(nrow, ncol), loc=(3, 5), colspan=1, rowspan=1)
        ax_em = plt.subplot2grid(shape=(nrow, ncol), loc=(2, 6), colspan=2, rowspan=2)

        # -------get object properites from photcat---------------------
        if df_photcat is not None:
            ra = df_photcat[df_photcat['id'] == args.id]['ra'].values[0]
            dec = df_photcat[df_photcat['id'] == args.id]['dec'].values[0]
        else:
            sys.exit(f'Photcat file does not exist in {phot_catalog_file}. So where do I get the RA/DEC from?')

        # --------cutout and plot the direct image---------------
        ax_dirimg = cutout2d_and_plot(ra, dec, direct_data, direct_wcs, ax_dirimg, args, cmap=dirimg_cmap)

        # --------cutout and plot the grism images of r & c orients---------------

        # --------model the 2D continuum and plot---------------

        # --------plot the continuum subtracted 2D spectra of r & c orients---------------

        # --------plot the emission maps of r & c orients---------------

        # --------rotate and plot the combined the emission map---------------

        # ------------saving the full figure--------------------------
        figname = extract_dir / f'{args.id:05d}_extraction.png'
        fig.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure as {figname}')
        plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
