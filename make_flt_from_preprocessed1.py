'''
    Filename: make_flt_from_preprocessed1.py
    Notes: Makes grism flt files from pre-processed *rate.fits files, for a given field
           This script is heavily based on grizli-notebooks/glass-niriss-wfss.ipynb (NB1)
    Author : Ayan
    Created: 10-07-24
    Example: run make_flt_from_preprocessed1.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par008
             run make_flt_from_preprocessed1.py --Par008
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

    # --------drizzle full mosaics-----
    if len(glob.glob('*failed')) > 0: os.remove('*failed')
    kwargs = get_yml_parameters()

    mosaic_args = kwargs['mosaic_args']
    mosaic_args['fill_mosaics'] = False

    # Set the mosaic pixel scale here
    mosaic_args['wcs_params']['pixel_scale'] = 0.04 # native 0.065 for NIRISS
    mosaic_args['half_optical_pixscale'] = True

    mosaic_args['ir_filters'] = ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']
    # mosaic_args['optical_filters'] = ['F115W','F150W','F200W'] # NIRCam

    mosaic_args['wcs_params']['pad_reference'] = 6  # small padding around edge, arcsec
    kwargs['mosaic_drizzle_args']['static'] = False

    auto_script.make_combined_mosaics(root, mosaic_args=mosaic_args, mosaic_drizzle_args=kwargs['mosaic_drizzle_args'])

    #!imsize -d -n 5 {root}-f*_drz_sci.fits # Show image dimensions
    '''
    # -----------making rgb mosaic (if 3 filters present)------------------
    slx, sly, rgb_filts, fig = auto_script.field_rgb(root=root, scl=8, HOME_PATH=None, xsize=6, force_rgb=['f200w', 'f150w', 'f115w'],
                                                     show_ir=True, rgb_scl=[1.65, 1.33, 1], gzext='*', output_format='png', suffix='.field', full_dimensions=4 # 4x4 in png relative to f444
                                                     )
    '''

    # -------catalog detection image--------------
    auto_script.make_filter_combinations(root, weight_fnu=True, min_count=1, filter_combinations={'ir': ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']})

    # -------source detection and aperture photometry--------------
    utils.set_warnings()
    phot = auto_script.multiband_catalog(field_root=root, detection_filter='ir', get_all_filters=True)
    '''
    # ------compare astrometry (cannot do without another source catalog)-----------------
    inp = utils.read_catalog('../abell2744_ip_2008_20220620_g3sw.cat', sextractor=True)

    idx, dr, dx, dy = inp.match_to_catalog_sky(phot, get_2d_offset=True)
    inp = inp[idx]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].hist(dr.value, bins=np.linspace(0, 1, 64))

    hasm = dr.value < 1

    axes[1].scatter(dx.value[hasm], dy.value[hasm], alpha=0.1)
    for ax in axes: ax.grid()
    '''
    # ----------make grism contamination model: only works for ONE filter, rather than combining different filers-----------------------------
    files = glob.glob('*GrismFLT.fits')

    if len(files) == 0:
        grism_prep_args = kwargs['grism_prep_args'] # Which filter to use as direct image?  Will try in order of the list until a match is found.
        grism_prep_args['refine_niter'] = 0 # For now, turn off refining contamination model with polynomial fits; to turn on set to 2
        grism_prep_args['init_coeffs'] = [1.0] # Flat-flambda spectra
        grism_prep_args['mask_mosaic_edges'] = False
        grism_prep_args['gris_ref_filters'] = {'GR150R': ['F115W', 'F150W', 'F200W'], 'GR150C': ['F115W', 'F150W', 'F200W']}
        grism_prep_args['force_ref'] = glob.glob(f'{root}-*_dr*_sci.fits*')[0]

        # Fairly bright for speedup, these can be adjusted based on how deep the spectra/visits are
        grism_prep_args['refine_mag_limits'] = [18, 22]
        grism_prep_args['prelim_mag_limit'] = 24

        grp = auto_script.grism_prep(field_root=root, pad=800, **grism_prep_args)

        grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), catalog='{0}-ir.cat.fits'.format(root), cpu_count=-1, sci_extn=1, pad=800)
    else:
        os.chdir(args.input_dir / args.field / 'Extractions')
        grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), catalog='{0}-ir.cat.fits'.format(root), cpu_count=-1, sci_extn=1, pad=800)

    # ---------------------------------------
    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
