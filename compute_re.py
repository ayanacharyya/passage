'''
    Filename: compute_re.py
    Notes: Computes effective radius for all objects in a given PASSAGE field/s, and stores the result as a FITS table
    Author : Ayan
    Created: 16-01-26
    Example: run compute_re.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run compute_re.py --field Par28
'''

from header import *
from util import *
from make_diagnostic_maps import get_re, get_offsets_from_center

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_niriss_psf(pix_size_arcsec, filter):
    '''
    Computes NIRISS PSF for a given filter, at a given pixel scale
    Returns stacked 2D line map and 2D uncertainty map
    '''
    niriss = webbpsf.NIRISS()
    niriss.filter = filter
    niriss.pixelscale = pix_size_arcsec
    psf = niriss.calc_psf(fov_arcsec=0.5, oversample=1)
    print('\n\n')
    
    return psf

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    filter_for_re = 'F150W'

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
        print(f'Starting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / 'catalogs'
        outfilename = output_dir / f'{args.field}_re_list.fits'

        if os.path.exists(outfilename) and not args.clobber:
            print(f'Result file for {args.field} already exists as {outfilename}, so skipping this field.')
            continue

        df_re = pd.DataFrame(columns=['id', 'redshift', 're_kpc', 're_arcsec'])

        # ---------read the photometric catalog file--------------------
        df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
        id_arr = df['id'].values
                
        # ------------looping over the objects-----------------------
        get_psf = True
        for index2, this_id in enumerate(id_arr):
            args.id = this_id
            print(f'\tCommencing ({index2 + 1}/{len(id_arr)}) ID {args.id}..')

            # ------determining directories and filenames---------
            full_fits_file = product_dir / 'full' / f'{args.field}_{args.id:05d}.full.fits'
            maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'

            if os.path.exists(maps_fits_file): # if the fits files are in maps/
                full_filename = maps_fits_file
            elif os.path.exists(full_fits_file): # if the fits files are in full/
                full_filename = full_fits_file
            else:
                print(f'Could not find {full_fits_file} or {maps_fits_file}, so skipping it.')
                continue

            # ------------read in maps files--------------------------------
            full_hdu = fits.open(full_filename)

            # ----------determining object parameters------------
            args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
            args.z = full_hdu[0].header['REDSHIFT']

            # --------determining true center of object---------------------
            args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W')

            # --------determining effective radius of object---------------------
            if get_psf:
                psf = get_niriss_psf(args.pix_size_arcsec, filter_for_re) # computing the PSF once (assuming the pixel scale is same for the full sample)
                get_psf = False
            
            try:
                re_kpc, re_arcsec = get_re(full_hdu, args, filter=filter_for_re, psf=psf)
            except:
                print(f'Could not compute R_e for object {args.id}. Skipping this object.')
                continue
        
            # -------appending to dataframe--------------
            df_re.loc[len(df_re)] = [args.id, args.z, re_kpc, re_arcsec]
        
        # ----------write out the dataframe--------------
        Table.from_pandas(df_re).write(outfilename, format='fits', overwrite=True)

        print(f'Written {outfilename}, completing field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
