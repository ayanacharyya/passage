'''
    Filename: get_scaling_data.py
    Notes: Extracts scaling data for OIII emission line from FITS files
    Author : Ayan
    Created: 25-08-26
    Example: run get_scaling_data.py --system ssd --do_all_fields --do_all_obj --cut_z_flag 4
'''

from header import *
from util import *
setup_plot_style()
from make_diagnostic_maps import trim_image, get_offsets_from_center
from stack_emission_maps import get_emission_line_map
from make_sfms_bins import get_binned_df, required_lines, sfms

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    
    # ------------reading and binning dataframe-------------
    df, bin_list, args = get_binned_df(args, required_lines=required_lines, sfms=sfms)
    scaling_data = []
    scaling_line = 'OIII'

    # ------determining directories and filenames---------
    for index, obj in df.iterrows():
        args.field = obj['field']
        args.id = obj['id']
        print(f'\nCommencing ID {args.field}:{args.id} ({index+1}/{len(df)})..')
        product_dir = args.input_dir / args.field / 'Products'
        maps_fits_file = product_dir / 'maps_psf_matched' / f'{args.field}_{args.id:05d}.maps.fits'
        full_hdu = fits.open(maps_fits_file)

        args.z = full_hdu[0].header['REDSHIFT']
        args.distance = cosmo.luminosity_distance(args.z)
        args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
        args.pix_size_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value
        args.pix_size_re = (2 * args.re_limit) / args.npix_side
        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
        
        # --------determining true center of object rom direct image---------------------
        args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W', silent=not args.debug_align)

        # ---------------segmentation map---------------
        segmentation_map = full_hdu['SEG'].data
        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
        segmentation_map = ndimage.shift(segmentation_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
        args.segmentation_map = segmentation_map        

        # ----------getting scaling line flux----------
        if scaling_line in args.available_lines:
            _, _, line_int = get_emission_line_map(scaling_line, full_hdu, args, dered=False, silent=True)
        else:
            print(f'Since {scaling_line} not available for {args.id}, not putting this in stack')
            line_int = ufloat(-99, 0)
        scaling_data.append({'field': args.field, 'id': args.id, 'scaling_flux_sum': line_int.n, 'scaling_flux_speccat': obj[f'flux_{scaling_line}']})

    df_scaling = pd.DataFrame(scaling_data)
    args.stacking_dir = args.output_dir / f'stacking{args.required_lines_text}{args.cut_z_flag_text}'
    df_scaling.to_csv(args.stacking_dir / f'scaling_data.csv', index=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
