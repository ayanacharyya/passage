'''
    Filename: make_rgb_image_cy5.py
    Notes: Makes RGB images using COSMOS-Web filters, for JWST Cycle 5 proposal targets
    Author : Ayan
    Created: 15-10-25
    Example: run make_rgb_image_cy5.py --id 119887
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def annotate_axes(ax, xlabel, ylabel, args, label='', clabel='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''
    ax.text(0.05, 0.9, label, c='k', fontsize=args.fontsize/args.fontfactor, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9), transform=ax.transAxes)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(xlabel, fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(ylabel, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    if not hide_cbar and p is not None:
        cax = inset_axes(ax, width="5%", height="100%", loc='right', bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.set_label(clabel, fontsize=args.fontsize)

        cbar.locator = ticker.MaxNLocator(integer=cticks_integer, nbins=4)#, prune='both')
        cbar.update_ticks()
        if hide_cbar_ticks:
            cbar.ax.set_yticklabels([])
        else:
            cbar.ax.tick_params(labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_cutout(filename, pos, size, ext=0):
    '''
    Return a cutout from a given filename of a fits image, around a given position within a given extent,
    Returns the 2D cutout as a 2D array
    '''
    data = fits.open(filename)

    image = data[ext].data
    source_header = data[ext].header
    wcs = pywcs.WCS(source_header)

    cutout = Cutout2D(image, pos, size, wcs=wcs)
    cutout_data = cutout.data

    return cutout_data

# --------------------------------------------------------------------------------------------------------------------
def get_direct_image(filename, ra, dec, cutout_width_arcsec, ext=0):
    '''
    Loads the direct image for a given filter for a given object
    Returns the image
    '''
    #exptime = fits.open(filename)[0].header['EXPTIME']
    pos = SkyCoord(ra, dec, unit = 'deg')
    size = cutout_width_arcsec * u.arcsec
    image = get_cutout(filename, pos, size, ext=ext)
 
    return image

# --------------------------------------------------------------------------------------------------------------------
def plot_rgb_image(obj, filters, ax, cutout_width_arcsec, args, colors=['r', 'g', 'b'], hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, annotate_fov_arcsec=None, hide_filter_names=False, zcol='redshift'):
    '''
    Plots the direct image as an RGB image, combining the given list of filters, for a given object, in a given axis
    Returns the axis handle
    '''
    # -------get image for each filter---------
    image_arr = []
    for index, filter in enumerate(filters):
        print(f'Reading in direct image for filter {filter}..')
        filename = args.bg_image_dir / f'mosaic_nircam_{filter.lower()}_COSMOS-Web_60mas_{obj["tile"]}_v1.0_sci.fits.gz'
        image = get_direct_image(filename, obj['ra'], obj['dec'], cutout_width_arcsec)
        image_arr.append(image)
        if not hide_filter_names: ax.text(0.05, 0.95 - index * 0.1, f'{filter}', c=colors[index], fontsize=args.fontsize / args.fontfactor, ha='left', va='top', transform=ax.transAxes, zorder=10)

    # -------create RGB image---------
    pctl, maximum = 99.9, 0.
    for img in image_arr:
        val = np.percentile(img, pctl)
        if val > maximum: maximum = val

    rgb_image = make_rgb(image_arr[0], image_arr[1], image_arr[2], interval=ManualInterval(vmin=0, vmax=maximum), stretch=SqrtStretch())

    # -------plot RGB image---------
    offset = 0 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
    extent = (-cutout_width_arcsec/2 - offset, cutout_width_arcsec/2 - offset, -cutout_width_arcsec/2 - offset, cutout_width_arcsec/2 - offset)
    p = ax.imshow(rgb_image, origin='lower', extent=extent, alpha=1)
    ax.set_aspect('auto') 

    #ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.text(0.05, 0.95, f'COSMOS2025 ID {obj[idcol]}', c='w', fontsize=args.fontsize / args.fontfactor, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.95, 0.05, f'z={obj[zcol]:.2f}', c='w', fontsize=args.fontsize / args.fontfactor, ha='right', va='bottom', transform=ax.transAxes)

    # ----------annotate axis---------------
    if annotate_fov_arcsec is not None:
        offset_x, offset_y, angle = rot_off_dict[obj[idcol]]
        limit = annotate_fov_arcsec / 2
        ax.add_patch(plt.Rectangle((offset_x - limit, offset_y - limit), 2 * limit, 2 * limit, lw=2, color='orangered', fill=False, angle=angle))

    ax.set_xlim(-cutout_width_arcsec/2, cutout_width_arcsec/2)  # arcsec
    ax.set_ylim(-cutout_width_arcsec/2, cutout_width_arcsec/2)  # arcsec
    ax.set_aspect('equal')

    if not skip_annotate: ax = annotate_axes(ax, 'arcsec', 'arcsec', args, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return ax, rgb_image

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ----------values to change----------
    cutout_width_arcsec = 4
    #rgb_filters = ['F444W', 'F277W', 'F150W', 'F115W']
    rgb_filters = ['F444W', 'F277W', 'F150W']
    #rgb_filters = ['F277W', 'F150W', 'F115W']
    #rgb_filters = ['F444W', 'F277W', 'F115W']
    #rgb_filters = ['F444W', 'F150W', 'F115W']
    args.fontfactor = 1
    rot_off_dict = {432033: [0, 0, 0], 549577: [0, 0, 0], 270854: [0, 0, 0], 119887: [0, 0, 0]} # [offset in ra (arcsec), offset in dec (arcsec), PA (deg)]

    # ----------nomenclature----------
    idcol = 'ID_c2025' # column name of ID in the photometry catalog
    target_list = args.input_dir / 'COSMOS/zCOSMOS_2025.fits'
    fig_dir = HOME / 'Documents/writings/telescope_proposals/JWSTCy5/Figures' # this is where the ETC related outputs would be stored
    
    if args.bg_image_dir is None: args.bg_image_dir = args.input_dir / 'COSMOS/imaging'
    else: args.bg_image_dir = Path(args.bg_image_dir)
    
    # --------reading in files-------------------
    df_targets = Table.read(target_list).to_pandas()    
    df_targets['tile'] = df_targets['tile'].astype(str)
    
    # ------declaring the figure--------
    nrow, ncol = 2, 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(6, 6), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.98, wspace=0.01, hspace=0.01)

    # ------looping over targets---------
    for index, id in enumerate(args.id):
        print(f'Doing id {id}: {index+1} of {len(args.id)}..')
        ax = axes[int(index/ncol)][index % ncol]
        obj = df_targets[df_targets[idcol] == id].iloc[0]
        ax, rgb_image = plot_rgb_image(obj, rgb_filters, ax, cutout_width_arcsec, args, zcol='redshift_zcosmos', skip_annotate=False, hide_xaxis=True, hide_yaxis=True, annotate_fov_arcsec=3, hide_filter_names=True)

    # -------saving the figure---------
    figname = fig_dir / f'all_targets_rgb.png'
    fig.savefig(figname)
    print(f'Saved {figname}\n')

    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
