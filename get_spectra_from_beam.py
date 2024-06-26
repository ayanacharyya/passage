'''
    Filename: get_spectra_from_beam.py
    Notes: Makes 1D and 2D spectra from beam files, for a given object in a given field
           This script is heavily based on grizli-notebooks/NIRISS-Demo-grizli.ipynb and grizli-notebooks/JWST/glass-niriss-wfs.ipynb
    Author : Ayan
    Created: 11-06-24
    Last modified: 11-06-24
    Example: run get_spectra_from_beam.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run get_spectra_from_beam.py --id 3617 --line_list OII,Hb,OIII,Ha,Ha+NII,PaA,PaB
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    args.input_dir = args.input_dir / args.field / 'beams'
    args.output_dir = args.output_dir / args.field / f'{args.id:05d}'
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.field
    fit_args_file = args.output_dir / f'tfit_{args.id:05d}.npy'
    os.chdir(args.input_dir)

    # -----generating fit_params dict-----------
    pline = {'kernel': 'square', 'pixfrac': 0.5, 'pixscale': 0.04, 'size': 8, 'wcs': None}
    dummy_args = auto_script.generate_fit_params(pline=pline, field_root=root, min_sens=0.0, min_mask=0.0, include_photometry=args.include_photometry, save_file=str(fit_args_file), full_line_list=args.line_list)
    '''
    #### the following part is from Peter, and atm it is redundant #######
    # ----loading multibeam file---------------
    #beams = grp.get_beams(args.id, size=50, min_mask=0, min_sens=0, show_exception=True, beam_id='A')

    mb = MultiBeam(str(args.input_dir / f'{root}_{args.id:05d}.beams.fits'), group_name=root, fcontam=0.2, min_sens=0.0, min_mask=0.0)
    mb.fit_trace_shift()

    fig_1d = mb.oned_figure()
    fig_2d = mb.drizzle_grisms_and_PAs(size=32, scale=1.0, diff=False)

    mb.write_master_fits()

    # -----fitting redshift---------------
    mb_out, st, fit_table, template_fit, line_hdu = fitting.run_all_parallel(args.id, zr=[args.zmin, args.zmax], verbose=~args.silent, get_output_data=True, args_file=str(fit_args_file), group_name=str(args.output_dir / args.field))

    # -------------regenerate the emission line maps----------------
    line_hdu = mb.drizzle_fit_lines(template_fit, pline, force_line=utils.DEFAULT_LINE_LIST, save_fits=False, mask_lines=True, min_line_sn=1, mask_sn_limit=3, verbose=~args.silent, get_ir_psfs=False)
    ##############################################
    '''
    # -----fitting redshift---------------
    _fit = fitting.run_all_parallel(args.id, zr=[args.zmin, args.zmax], verbose=~args.silent, get_output_data=True, args_file=str(fit_args_file), group_name=str(args.output_dir / args.field))

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
