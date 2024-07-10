'''
    Filename: preprocess_raw_data1.py
    Notes: Pre process raw *rate.fits files, for a given field
           This script is heavily based on grizli-notebooks/glass-niriss-wfss.ipynb
    Author : Ayan
    Created: 10-07-24
    Example: run preprocess_raw_data1.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par008
             run preprocess_raw_data1.py --Par008
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------determining directories and global variables---------
    subdirectories_required = ['Prep', 'RAW', 'MAST', 'Products', 'Extractions']
    for this_dir in subdirectories_required:
        sub_dir = args.input_dir / args.field / this_dir
        sub_dir.mkdir(parents=True, exist_ok=True)
    args.download_dir = args.input_dir / args.field / 'MAST'
    args.raw_dir = args.input_dir / args.field / 'RAW'
    args.work_dir = args.input_dir / args.field / 'Prep'
    root = args.field

    # -----copy files over to working directory-----------
    os.chdir(args.download_dir)
    files = glob.glob('*rate.fits')
    for file in files: shutil.copyfile(os.path.join(args.download_dir, file), os.path.join(args.raw_dir, file)) # copy files from MAST/ to Prep/ otherwise MAST/ files will be overwritten
    os.chdir(args.work_dir)
    files = glob.glob(str(args.raw_dir / '*rate.fits'))
    files.sort()

    # ------changing some keywords and fixing bad pixels in files----------
    for file in files:
        _ = jwst_utils.set_jwst_to_hst_keywords(file, reset=True)

        im = fits.open(file)
        bad = im['SCI'].data > 1.e8
        print(f'Total bad pixels in {file} = {bad.sum()}')
        im['SCI'].data[bad] = 0
        im['DQ'].data[bad] = 2048

        im.writeto(file, overwrite=True)
        im.close()

    #print(utils.get_flt_info(files=files))

    # --------parse visits--------------------
    visits, all_groups, info = auto_script.parse_visits(field_root=root, files=files, combine_same_pa = False)

    # ----------------------------------------
    kwargs = get_yml_parameters()
    # Parameter lists
    visit_prep_args = kwargs['visit_prep_args']
    preprocess_args = kwargs['preprocess_args']

    # Maximum shift for "tweakshifts" relative alignment
    IS_PARALLEL = True

    tweak_max_dist = (5 if IS_PARALLEL else 1)
    if 'tweak_max_dist' not in visit_prep_args: visit_prep_args['tweak_max_dist'] = tweak_max_dist

    # Fit and subtract a SExtractor-like background to each visit
    visit_prep_args['imaging_bkg_params'] = {'bh': 256, 'bw': 256, 'fh': 3, 'fw': 3, 'pixel_scale': 0.1, 'get_median': False}

    # Alignment reference catalogs, searched in this order
    visit_prep_args['reference_catalogs'] = ['LS_DR9', 'PS1', 'GAIA', 'SDSS', 'WISE']

    # ----------Preview what will be processed------------------
    colors = [c['color'] for c in plt.rcParams['axes.prop_cycle']]

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    for i, v in enumerate(visits):
        sr = utils.SRegion(v['footprint'])

        ax.scatter(*sr.centroid[0], marker='.', c=colors[i])

        for patch in sr.patch(ec=colors[i], fc='None', alpha=0.5, label=v['product']):
            ax.add_patch(patch)

    ax.set_aspect(1. / np.cos(ax.get_ylim()[0] / 180 * np.pi))  # square with cos(dec)
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.legend()
    ax.grid()

    # ------pre-processing files----------
    '''
    inp = utils.read_catalog('../abell2744_ip_2008_20220620_g3sw.cat', sextractor=True)
    inp['ra'] = inp['ALPHAWIN_J2000']
    inp['dec'] = inp['DELTAWIN_J2000']
    radec_file = os.path.join(args.input_dir, 'abell2744_ip_2008_20220620_g3sw.radec')
    prep.table_to_radec(inp, radec_file)
    preprocess_args['parent_radec'] = radec_file # First visit will be aligned to this, subsequent visits that overlap will be aligned to each other.
    '''
    preprocess_args['master_radec'] = None # Query all-sky catalogs (LS DR9, PS1, GAIA, etc) for first visit. Subsequent visits that overlap will be aligned to each other.

    im = fits.open(files[0])
    os.environ['CRDS_CONTEXT'] = im[0].header['CRDS_CTX']

    auto_script.preprocess(field_root=root, HOME_PATH=args.input_dir, visit_prep_args=visit_prep_args, **preprocess_args)

    # -----------Update the visits file with the new exposure footprints-----
    visit_file = auto_script.find_visit_file(root=root)
    print('Update exposure footprints in {0}'.format(visit_file))
    _ = auto_script.get_visit_exposure_footprints(root=root, check_paths=['./'])

    #! cat *shifts.log # Tweakshifts
    #! cat *wcs.log # These derived shifts are non-zero because we shifted the reference catalog above

    # -----------visualise the pre-processed files---------------------------
    files = glob.glob('*sci.fits')
    for file in files:
        im = fits.open(file)
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        vm = np.nanpercentile(im[0].data, [5, 95])
        ax.imshow(im[0].data, vmin=-0.1 * vm[1], vmax=vm[1], cmap='gray_r')
        ax.set_aspect(1)
        ax.set_title(file)
        ax.axis('off')

    # -----------------------------------------------------------------------
    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
