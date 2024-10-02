'''
    Filename: run_passagepipe.py
    Notes: Runs Vihang's PASSAGEPipe from start to finish, for a given field
    Author : Ayan
    Created: 15-07-24
    Example: run run_passagepipe.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run run_passagepipe.py --Par 009 --zmin 1.5 --zmax 5 --start_id 70
             run run_passagepipe.py --field Par6 --zmin 0.5 --zmax 3.8 --filters F115W,F150W --start_step 5 --start_id 121
             run run_passagepipe.py --field Par61 --zmin 0.1 --zmax 8 --filters F115W,F150W
             run run_passagepipe.py --field Par27 --zmin 0.1 --zmax 8 --filters F115W,F150W,F200W
             run run_passagepipe.py --field Par40 --zmin 0.1 --zmax 8 --filters F115W,F150W,F200W
             run run_passagepipe.py --field Par61 --zmin 0.1 --zmax 8 --filters F115W,F150W --remake_figures
'''

from header import *
from util import *
start_time = datetime.now()

import passagepipe as pp
import passagepipe.pipeline as pipe
from passagepipe.post_process import makeSummaryFigures, createTarballs

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ---------determining config file parameters---------------------
    if args.mag_lim is None: args.mag_lim = 28
    filters_text = '\n      '.join([f'- {item}' for item in args.filters])
    force_ref_file = f'{args.field}-ir_drz_sci.fits'

    if args.start_step <=1:
        args.do_download = True
        args.download_from_mast = True
        args.redo_level1 = True
    if args.start_step <=2:
        args.do_prep = True
    if args.start_step <=3:
        args.do_image = True
    if args.start_step <=4:
        args.do_grism = True
    if args.start_step <=5:
        args.do_extract = True
    if args.start_step <=6:
        args.do_post = True

    # ----------replacing keywords in config file template to make the actual config file---------
    config_template = args.input_dir / 'config_files/passage_config_template.txt'
    config_filename = args.input_dir / f'config_files/passage_config_{args.field}.yml'

    replacements = {'FIELD': args.field, 'FILTERS': filters_text, 'DATA_DIR': args.input_dir, \
                    'DO_DOWNLOAD': args.do_download, 'DO_PREP': args.do_prep, 'DO_IMAGE': args.do_image, 'DO_GRISM': args.do_grism, 'DO_EXTRACT': args.do_extract, 'DO_POST': args.do_post, 'DO_UPLOAD': args.do_upload, \
                    'DOWNLOAD_FROM_MAST': args.download_from_mast, 'CLOBBER_DOWNLOAD': args.clobber_download, 'REDO_LEVEL1': args.redo_level1, \
                    'FORCE_REF': force_ref_file, 'MAGLIM': args.mag_lim, 'MAGMIN': args.magmin, 'MAGMAX': args.mag_lim, \
                    'CLOBBER_EXTRACT': args.clobber, 'ZMIN': args.zmin, 'ZMAX': args.zmax}  # keywords to be replaced in template config file

    with open(config_template) as infile, open(config_filename, 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(str(src), str(target))
            outfile.write(line)  # replacing and creating new jobscript file

    # ------------reading in config file-----------------------
    config = pp.PassageConfig()
    config.updateConfig(config_filename, initialize=True)

    try:
        objList = pp.utils.readPhotCatalogDefault(CONFIG=config)["NUMBER"]
        objList = objList[args.start_id : args.stop_id]
    except FileNotFoundError:
        objList = None

    # ------------running the pipeline-----------------------
    if args.dry_run:
        print('\nThis is just a dry-run, NOT calling the PASSAGEPipe. Omit the --dry_run option in order to actually run the pipeline. Print config if you want to check the loaded configurations.\n')
    elif args.remake_figures: # to only remake the figures and tarball them
        print('\n Only re-making the figures in Products/plots/ and re-making the tarball, assuming all other steps have already been done and all other files are already present.')
        product_subdirs = ['plots']
        makeSummaryFigures(CONFIG=config)
        createTarballs(CONFIG=config)
    else:
        pipe.run_pipeline(CONFIG=config, objList=objList)

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
