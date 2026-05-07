'''
    Filename: make_rgb_image.py
    Notes: Makes RGB images using multi-filter imaging, for a given filter
    Author : Ayan
    Created: 07-05-26
    Example: run make_rgb_image.py --system ayan_ssd --field Par034
             run make_rgb_image.py --system ayan_ssd --field Par016,Par019,Par021,Par034,Par040,Par042,Par044
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ---------------looping over fields--------------------
    for index2, args.field in enumerate(args.field_arr):
        print(f'\nDoing field ({index2 + 1}/{len(args.field_arr)}) {args.field}')
        
        # --------get filenames-------------------
        image_dir = args.input_dir / args.field
        if not os.path.exists(image_dir):
            print(f'\t{image_dir} does not exist, so skipping this field..')
            continue
        
        filters = ['F814W',  # Red
                'F625W',  # Green
                'F475W',   # Blue
                ]
        files = [glob.glob(str(image_dir) + f'/*{args.field}*{filter}*_drc_sci.fits')[0] for filter in filters]
        if len(files) < 3:
            print(f'\tOnly {len(files)} (< 3) have been found for this field, so skipping it..')
            continue

        # --------loading redest image as reference-----------------
        with fits.open(files[0]) as hdu_r:
            image_r = hdu_r[0].data
            header_r = hdu_r[0].header
            photflam_r = header_r['PHOTFLAM']
            wcs_r = pywcs.WCS(header_r)
            target_shape = image_r.shape

        # ----------setting up figure-------------
        fig = plt.figure(figsize=(8, 7), layout='constrained')
        ax = fig.add_subplot(111, projection=wcs_r)
        colors = ['r', 'g', 'b']

        # -------get image for each filter---------
        images = {}
        for index, filename in enumerate(files):
            filter = filters[index]
            
            print(f'\tReading in direct image {filename}..')
            with fits.open(filename) as hdul:
                #data, _ = reproject_interp(hdul[0], wcs_r, shape_out=target_shape)
                data = hdul[0].data.copy()
                header = hdul[0].header

                # If header['UNIT'] is 'COUNTS', divide by EXPTIME, If it's already 'ELECTRONS/S', skip this
                if header.get('BUNIT', '').upper() == 'COUNTS':
                    exptime = header.get('EXPTIME', 1.0)
                    data /= exptime

                photflam = header.get('PHOTFLAM', 1.0)
                data = data * photflam / photflam_r
            
            bkg = np.nanmedian(data)
            images[filter] = data - bkg

            ax.text(0.05, 0.95 - index * 0.05, f'{filter}', c=colors[index], fontsize=args.fontsize, ha='left', va='top', transform=ax.transAxes, zorder=10)

        # -------create RGB image---------
        rgb_image = make_lupton_rgb(images[filters[0]], images[filters[1]], images[filters[2]], Q=5, stretch=0.1)
        
        # -------plot RGB image---------
        p = ax.imshow(rgb_image, origin='lower')
        ax.set_aspect('auto')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # -------saving the figure---------
        figname = image_dir / f'rgb.png'
        save_fig(fig, image_dir, f'{args.field}_rgb_{",".join(filters)}.png', args, dpi=1000)
    
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
