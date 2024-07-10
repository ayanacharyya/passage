'''
    Filename: make_flt_from_preprocessed2.py
    Notes: Makes grism flt files from pre-processed *rate.fits files, for a given field
           This script is heavily based on grizli-notebooks/JWST/grizli-niriss-2023.ipynb (NB2)
    Author : Ayan
    Created: 10-07-24
    Example: run make_flt_from_preprocessed2.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par008
             run make_flt_from_preprocessed2.py --Par008
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    args.raw_dir = args.input_dir / args.field / 'RAW'
    args.work_dir = args.input_dir / args.field / 'Prep'
    root = args.field
    os.chdir(args.work_dir)
    files = glob.glob(str(args.raw_dir / '*rate.fits'))
    files.sort()

    # ------making mosaics----------
    res = visit_processor.res_query_from_local(files=files)
    is_grism = np.array(['GR' in filt for filt in res['filter']])

    hdu = utils.make_maximal_wcs(files=files, pixel_scale=0.04, pad=4, get_hdu=True, verbose=False)
    ref_wcs = pywcs.WCS(hdu.header)
    _ = visit_processor.cutout_mosaic(args.field,
                                      res=res[~is_grism], # Pass the exposure information table for the direct images
                                      ir_wcs=ref_wcs,
                                      half_optical=False, # Otherwise will make JWST exposures at half pixel scale of ref_wcs
                                      kernel='square',  # Drizzle parameters
                                      pixfrac=0.8,
                                      clean_flt=False, # Otherwise removes "rate.fits" files from the working directory!
                                      s3output=None,
                                      make_exptime_map=False,
                                      weight_type='jwst',
                                      skip_existing=True,
                                      )

    # ----------making source catalog----------------
    _cat = prep.make_SEP_catalog('indef-01571-292.0-nis-f200w-clear', threshold=1.2) # this part still does not work

    # ----------making grism models------------------
    pad = 800
    un = utils.Unique(res['pupil']) # Process by blocking filter
    all_grism_files = glob.glob(str(args.input_dir / args.field / 'Extractions') + '*GrismFLT.fits')
    all_grism_files.sort()

    if len(all_grism_files) == 0:
        grp = {}
        for filt in un.values:
            grism_files = ['{dataset}_rate.fits'.format(**row) for row in res[is_grism & un[filt]]]
            reference_file = glob.glob(f'{root}-{filt.lower()}_dr*_sci.fits*')[0]
            grp[filt] = auto_script.grism_prep(field_root=root, PREP_PATH='../Prep', EXTRACT_PATH='../Extractions/', gris_ref_filters={'GR150R': ['F115W', 'F150W', 'F200W'], 'GR150C': ['F115W', 'F150W', 'F200W']},
                                               refine_niter=0, refine_mag_limits=[18, 22], prelim_mag_limit=24, init_coeffs=[1], pad=pad, force_ref=reference_file,
                                               files=grism_files, model_kwargs={'compute_size': False, 'size': 48}, subtract_median_filter=False, use_jwst_crds=False, cpu_count=1)
    else:
        os.chdir(os.path.join(args.input_dir, args.field, 'Extractions'))

        file_filters = {}
        for f in all_grism_files:
            with fits.open(f) as im:
                filt = im[1].header['DFILTER'].split('-')[0]
                if filt not in file_filters:
                    file_filters[filt] = []

                file_filters[filt].append(f)

        filts = un.values
        grism_files = []
        for filt in filts:
            if filt in file_filters:
                grism_files += file_filters[filt]

        catalog = glob.glob(f'{root}-*.cat.fits')[0]
        seg_file = glob.glob(f'{root}-*_seg.fits')[0]

        grp = {}
        for filt in filts:
            if filt not in file_filters: continue

            print(filt, len(file_filters[filt]), ' '.join(file_filters[filt]))
            print('\n-------------------------\n')

            grp[filt] = multifit.GroupFLT(grism_files=file_filters[filt], direct_files=[], ref_file=None,  # automatically compute based on the filter
                                          seg_file=seg_file, catalog=catalog, cpu_count=1, sci_extn=1, pad=pad)

    # ------------drizzle combined images--------------
    utils.set_warnings()
    for k in grp:
        grp[k].drizzle_grism_models(root=root, kernel='square', scale=0.06, pixfrac=0.75)


    # -----------------------------------------------------
    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
