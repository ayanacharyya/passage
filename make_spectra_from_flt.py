'''
    Filename: make_spectra_from_flt.py
    Notes: Makes 1D and 2D spectra from grism flt files, for a given object in a given field
           This script is heavily based on grizli-notebooks/JWST/grizli-niriss-2023.ipynb (NB2)
    Author : Ayan
    Created: 11-06-24
    Example: run make_spectra_from_flt.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_spectra_from_flt.py --line_list OII,Hb,OIII,Ha,Ha+NII,PaA,PaB --id 3
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    args.raw_dir = args.input_dir / args.field / 'RAW'
    args.work_dir = args.input_dir / args.field / 'Extractions'
    os.chdir(args.work_dir)
    root = args.field

    # ---------determine file and filter names--------------------
    files = glob.glob(str(args.raw_dir / '*rate.fits'))
    res = visit_processor.res_query_from_local(files=files)
    un = utils.Unique(res['pupil']) # Process by blocking filter

    all_grism_files = glob.glob(os.path.join(args.input_dir, args.field, 'Extractions', '*GrismFLT.fits'))
    all_grism_files.sort()

    # ----------------read in the grism flt files--------------------
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
        grp[filt] = multifit.GroupFLT(grism_files=file_filters[filt], direct_files=[], ref_file=None, seg_file=seg_file, catalog=catalog, cpu_count=1, sci_extn=1, pad=800)

   # -------initialise redshift fit defaults-------------------
    pline = {'kernel': 'square', 'pixfrac': 0.5, 'pixscale': 0.04, 'size': 8, 'wcs': None}
    if args.line_list == 'all': args.linelist = ['Lya', 'OII', 'Hb', 'OIII', 'Ha', 'NII','Ha+NII', 'SII', 'SIII', 'PaD','PaG','PaB','HeI-1083','PaA']

    fit_args = auto_script.generate_fit_params(field_root=root,
                                           zr=[0.02, 8], # full redshift range to fit
                                           dz=[0.004, 0.0004], # two-pass redshift grid steps in (1+z)
                                           sys_err=0.03,
                                           include_photometry=False,
                                           fit_trace_shift=False, # this can help with some trace misalignment
                                           dscale=0.01,
                                           fwhm=500. * u.km / u.second, # velocity width of emission line templates
                                           full_line_list=args.line_list, # Make line maps of these
                                           min_sens=1.e-6, min_mask=1.e-6,  # masking around sensitivity edges
                                           mask_resid=True,
                                           fcontam=0.2,
                                           # include a contribution of the contamination model in the uncertainties
                                           pline=pline,  # line map parameters
                                           )

    # ------------prep for beam files--------------------------------
    size = 48
    id_arr = grp[list(grp.keys())[0]].catalog['NUMBER'] if args.do_all_obj else args.id

    # ------------make beams files--------------------------------
    for index, this_id in enumerate(id_arr):
        start_time2 = datetime.now()
        print(f'\nMaking beam file for id {this_id} which is {index+1} out of {len(id_arr)}..')

        if os.path.exists(f'{root}_{this_id:05d}.beams.fits') and not args.clobber:
            print(f'Skipping id {this_id} due to existing file {root}_{this_id:05d}.beams.fits')
            continue

        beams = []
        mbf = {}
        beam_groups = {}

        for filt in grp:
            beams_i = grp[filt].get_beams(this_id, size=size, min_mask=fit_args['min_mask'], min_sens=fit_args['min_sens'], mask_resid=False)
            msg = f'{filt} {len(beams_i)}'
            beams += beams_i

        print(this_id, len(beams))
        if len(beams) == 0: continue
        mb = multifit.MultiBeam(beams, **fit_args)

        mb.write_master_fits()

        for b in mb.beams: print(b.grism.filter, b.grism.pupil, b.beam.conf.conf_file)

        # -----fitting redshift---------------
        output_subdirectory = args.output_dir / args.field / f'{this_id:05d}'
        output_subdirectory.mkdir(parents=True, exist_ok=True)

        _fit = fitting.run_all_parallel(this_id, zr=[args.zmin, args.zmax], verbose=~args.silent, get_output_data=True, group_name=str(output_subdirectory / args.field))
        print(f'Completed id {this_id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
