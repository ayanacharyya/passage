'''
    Filename: make_spectra_from_beams.py
    Notes: Extracts 1D and 2D spectra from existing beams files, for a given object/s in a given field
           This script is heavily based on grizli-notebooks/JWST/grizli-niriss-2023.ipynb (NB2)
    Author : Ayan
    Created: 12-07-24
    Example: run make_spectra_from_beams.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_spectra_from_beams.py --line_list OII,Hb,OIII,Ha,Ha+NII,PaA,PaB --field Par9 --id 3 --zmin 1.5 --zmax 5
'''

from header import *
from util import *
from grizli.utils import GTable

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    args.work_dir = args.input_dir / args.field / 'Extractions'
    os.chdir(args.work_dir)

    # ------------read in catalog file--------------------------------
    if args.do_all_obj:
        catalog_file = args.input_dir / args.field / 'Prep' / f'{args.field}-ir.cat.fits'
        catalog = GTable.read(catalog_file)
        id_arr = catalog['NUMBER']
    else:
        id_arr = args.id

    # ------------extract spectra--------------------------------
    for index, this_id in enumerate(id_arr):
        start_time2 = datetime.now()
        print(f'\nFitting spectra for id {this_id} which is {index+1} out of {len(id_arr)}..')

        # -----fitting redshift---------------
        pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
        output_subdirectory = args.output_dir / args.field / f'{this_id:05d}{pixscale_text}'
        output_subdirectory.mkdir(parents=True, exist_ok=True)

        _fit = fitting.run_all_parallel(this_id, zr=[args.zmin, args.zmax], pline={'pixscale': args.pixscale}, verbose=~args.silent, full_line_list=args.line_list, get_output_data=True, group_name=str(output_subdirectory / args.field))
        print(f'Completed id {this_id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
