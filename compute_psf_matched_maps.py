'''
    Filename: compute_psf_matched_maps.py
    Notes: Computes effective radius for all objects in a given PASSAGE field/s, and stores the result as a FITS table
    Author : Ayan
    Created: 16-01-26
    Example: run compute_psf_matched_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run compute_psf_matched_maps.py --field Par28 --id 2822 --clobber
             run compute_psf_matched_maps.py --field Par28 --do_all_obj --clobber
             run compute_psf_matched_maps.py --system ssd --do_all_fields --do_all_obj --clobber --write_file
'''

from header import *
from util import *
from make_diagnostic_maps import get_re_from_extension, get_offsets_from_center
from make_sfms_bins import read_passage_sed_catalog

webbpsf.setup_logging(level='WARNING') # to suppress chatty-ness of webbpsf
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_niriss_psf(filter_name, fov_pixels, supersampling_factor=1):
    '''
    Computes NIRISS PSF in a given filter
    Returns 2D array
    '''
    if fov_pixels % 2 == 0: fov_pixels += 1
    niriss = webbpsf.NIRISS()

    if type(filter_name) == str:
        print(f'\n\t\tCreating PSF at wavelength {filter_name} filter with fov_pixels={fov_pixels}..')
        niriss.filter = filter_name    
        psf = niriss.calc_psf(fov_pixels=fov_pixels, oversample=supersampling_factor)
    elif type(filter_name) == np.float64 or type(filter_name) == float: # here 'filter' is actually the observed wavelength in microns, not the filter
        print(f'\n\t\tCreating PSF at wavelength {filter_name} microns with fov_pixels={fov_pixels}..')
        psf = niriss.calc_psf(monochromatic=filter_name * 1e-6, fov_pixels=fov_pixels, oversample=supersampling_factor)
    else:
        sys.exit(f'Unrecignised data type for filter_name={filter_name}; it can only be str or float')

    psf_array  = psf[0].data
    
    return psf_array

# --------------------------------------------------------------------------------------------------------------------
def match_to_psf(image, obs_wave, psf_cache, target_psf, args, supersampling_factor=1):
    '''
    Function to match a given 2D image from a given source PSF to taget_PSF of F200W
    Returns convolved image fo same shape
    '''
    fov_pixels = 2 * args.arcsec_limit / args.pix_size_arcsec

    if type(obs_wave) == np.float64 or type(obs_wave) == float:
        obs_wave = np.round(obs_wave, 3)
    
    if obs_wave in psf_cache:
        obs_psf = psf_cache[obs_wave]
    else:
        obs_psf = get_niriss_psf(obs_wave, fov_pixels, supersampling_factor=supersampling_factor)
        psf_cache[obs_wave] = obs_psf
    
    kernel = matching.create_matching_kernel(obs_psf, target_psf, window=window)
    psf_matched_image = convolve(image, kernel, boundary='extend')

    return psf_matched_image, psf_cache

# --------------------------------------------------------------------------------------------------------------------
def process_fits_extensions(full_hdu, target_psf, args):
    '''
    Processes each extension of a given fits file and PSF-matches it to a given target PSF
    Returns hdu list of new modified fits file, which can then be written to disk
    '''
    #match_list = ['DSCI', 'DWHT', 'LINE', 'CONTAM', 'CONTINUUM']
    match_list = ['LINE', 'CONTAM', 'CONTINUUM']
    new_hdul = fits.HDUList()

    psf_cache = {}

    for index, ext in enumerate(full_hdu):
        ext_name = ext.name.upper()
        if any(substring in ext_name for substring in match_list):
            print(f'\t\tPSF-matching extension: {ext.header["EXTVER"]} {ext_name} ({index + 1}/{len(full_hdu)})')
            data_to_match = ext.data.copy()
            
            # ----------re-center line map-----------
            data_to_match = ndimage.shift(data_to_match, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)

            # -----------PSF-matching-------------
            obs_line_wave = ext.header['WAVELEN'] / 1e4 # to get wavelength in microns
            psf_mached_image, psf_cache = match_to_psf(data_to_match, obs_line_wave, psf_cache, target_psf, args)
            
            new_ext = fits.ImageHDU(data=psf_mached_image, header=ext.header, name=ext_name)            
        else:
            print(f'\t\tCopying untouched: {ext_name} ({index + 1}/{len(full_hdu)})')
            if isinstance(ext, fits.PrimaryHDU):
                new_ext = fits.PrimaryHDU(data=ext.data, header=ext.header)
            elif isinstance(ext, fits.BinTableHDU):
                new_ext = fits.BinTableHDU(data=ext.data, header=ext.header, name=ext_name)
            else:
                new_ext = fits.ImageHDU(data=ext.data, header=ext.header, name=ext_name)

        new_hdul.append(new_ext)
    
    return new_hdul

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    filter_for_re = 'F150W'

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / 'Par[0-9][0-9][0-9]') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr

    # -------------setup PSF calculations------------------
    target_psf_wave = 1.984 # corresponding to F200W pivot wavelength, in microns
    window = matching.CosineBellWindow(alpha=0.35)
    get_target_psf = True

    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'Starting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'

        if not os.path.exists(product_dir):
            print(f'Could not find {product_dir}, so skipping this field {args.field}.')
            continue

        # ---------read the photometric catalog file--------------------
        if args.do_all_obj:
            df = GTable.read(product_dir / f'{args.field}_photcat.fits').to_pandas()
            id_arr = df['id'].values
        else:
            id_arr = args.id
                
        # ------------looping over the objects-----------------------
        for index2, this_id in enumerate(id_arr):
            args.id = this_id
            print(f'\t\tCommencing ({index2 + 1}/{len(id_arr)}) ID {args.id}..')

            # ------determining directories and filenames---------
            full_fits_file = product_dir / 'full' / f'{args.field}_{args.id:05d}.full.fits'
            maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'

            if os.path.exists(maps_fits_file): # if the fits files are in maps/
                full_filename = maps_fits_file
            elif os.path.exists(full_fits_file): # if the fits files are in full/
                full_filename = full_fits_file
            else:
                print(f'\t\tCould not find {full_fits_file} or {maps_fits_file}, so skipping it.')
                continue

            outfilename = Path(str(full_filename.parent).replace('maps', 'maps_psf_matched').replace('full', 'full_psf_matched')) / full_filename.name
            if os.path.exists(outfilename) and not args.clobber:
                print(f'\t\tResult file for {args.field}, {args.id} already exists as {outfilename}, so skipping this object.')
                continue

            # ------------read in maps files--------------------------------
            full_hdu = fits.open(full_filename)

            # ----------determining object parameters------------
            args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
            args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
            args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
            args.z = full_hdu[0].header['REDSHIFT']
            args.ra = full_hdu[0].header['RA']
            args.dec = full_hdu[0].header['DEC']
            
            # --------determining true center of object---------------------
            args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter=filter_for_re)

            # -----------determine PSFs for PSF matching-----------------
            if get_target_psf:
                fov_pixels = 2 * args.arcsec_limit / args.pix_size_arcsec
                target_psf = get_niriss_psf(target_psf_wave, fov_pixels=fov_pixels)
                get_target_psf = False

            # --------compute PSF-matched maps---------------------
            try:
                new_hdul = process_fits_extensions(full_hdu, target_psf, args)
                new_hdul.writeto(outfilename, overwrite=True)
                print(f'\n\t\tSuccessfully saved matched FITS to: {outfilename}')
            
            except Exception as e:
                print(f'\t\tCould not produce PSF-matched fits file for object {args.id} due to {e}. Skipping this object.')
                continue
            
        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')
    
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
