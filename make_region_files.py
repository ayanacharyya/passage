'''
    Filename: make_region_files.py
    Notes: Generate region files of all PASSAGE field footprints
    Author : Ayan
    Created: 26-06-24
    Last modified: 26-06-24
    Example: run make_region_files.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_region_files.py --survey zcosmos
             run make_region_files.py --survey cosmos2020 --split_regions_by_fields
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def df_to_region(df, outfilename, shape='box', label='', label_ra=150.47, label_dec=2.55, color='red'):
    '''
    Converts the footprint list into a (box) region file and saves the file
    '''
    output_dir = Path(os.path.split(outfilename)[0])
    output_dir.mkdir(parents=True, exist_ok=True) # creating the directory in whihc to place the output, unless already existing

    with open(outfilename, 'w') as f:
        f.write('fk5\n'
                'global color=' + color + '\n'
                f'text {label_ra} {label_dec} {{{label}}}\n')

    df['shape'] = shape

    # -----to specify field id names as text--------
    if 'passage' in label.lower(): df['id'] = '# text={' + df['id'].str[3:].astype(int).astype(str) + '}'
    else: df['id'] = '# text={' + df['id'].astype(str) + '}'

    if shape == 'box':
        df['ra_width'] = df['ra_width'].astype(str) + '\"' # to specify arcseconds for the width
        df['dec_width'] = df['dec_width'].astype(str) + '\"' # to specify arcseconds for the width
        df = df[['shape', 'ra', 'dec', 'ra_width', 'dec_width', 'pa', 'id']]
    elif shape == 'circle':
        df['radius'] = df['radius'].astype(str) + '\"'  # to specify arcseconds for the width
        df = df[['shape', 'ra', 'dec', 'radius', 'id']]

    df.to_csv(outfilename, header=False, index=False, sep='\t', mode='a', doublequote=False, quotechar="'")
    print(f'Saved region file {outfilename}')

# -------------------------------------------------------------------------------------------------------
def get_passage_footprints(args=None, filename=None):
    '''
    Reads and returns the footprint list for PASSAGE fields
    '''
    if filename is None: filename = args.input_dir / 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'

    NIRISS_WFSS_fov = [133, 133] # FoV in arcseconds

    df = pd.read_csv(filename)
    df = df[['Par#', 'Target RA (J2000)', 'Target Dec (J2000)']].dropna().reset_index(drop=True)
    df = df.rename(columns={'Par#':'id', 'Target RA (J2000)':'ra', 'Target Dec (J2000)':'dec'})

    df['ra_width'] = NIRISS_WFSS_fov[0] # arcseconds
    df['dec_width'] = NIRISS_WFSS_fov[1]
    df['pa'] = 0 # dummy value; in degrees

    return df # df only contains id, ra, dec, ra_width, dec_width, pa (orienation) for each field

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    output_dir = args.input_dir / 'footprints/region_files'

    if args.survey == 'passage':
        df = get_passage_footprints(args=args) # reading the passage footprints
        df_to_region(df, output_dir / 'PASSAGE_fields.reg', shape='box', label='JWST-PASSAGE', label_ra=150.47, label_dec=2.55, color='red')  # making the region files

    elif args.survey == 'zcosmos':
        df = read_zCOSMOS_catalog(args=args)
        df['radius'] = 1  # arcseconds # choice of how big a circle to draw around each RA,DEC combo found in the zCOSMOS catalog
        df_to_region(df, output_dir / 'zcosmos_objects.reg', shape='circle', label='zCOSMOS', label_ra=150.47, label_dec=2.55, color='cyan') # making the region files

    elif args.survey == 'cosmos2020':
        df = read_COSMOS2020_catalog(args=args)
        df['radius'] = 1  # arcseconds # choice of how big a circle to draw around each RA,DEC combo found in the zCOSMOS catalog
        df_to_region(df, output_dir / 'cosmos2020_objects.reg', shape='circle', label='COSMOS2020', label_ra=150.47, label_dec=2.55, color='cyan') # making the region files

        # -----------to save multiple region files split by whatever is contained within the PASSAGE fields------
        if args.split_regions_by_fields:
            dummy_dataset = regions.make_example_dataset(data='simulated') # creating dummy data to generate dummy wcs, needed for sky_region.contains()
            dummy_wcs = dummy_dataset.wcs

            sky_regions = sky_regions = Regions.read(output_dir / 'PASSAGE_fields.reg', format='ds9') # if this file does not exist, first run this script with "--survey passage" option, before attemtping --split_regions_by_fields
            df = read_COSMOS2020_catalog(args=args)

            sky_regions = np.array(sky_regions)[passage_fields_in_cosmos] # this is to reduce run time, based on pre-calculated index array above; for a fresh run, comment out this line

            for index, sky_region in enumerate(sky_regions):
                if not type(sky_region) == regions.shapes.text.TextSkyRegion: # consider it if it is not just text
                    label_text = f'Par{sky_region.meta["text"]:03}'
                    print(f'Splitting by fields. Doing field {label_text}, which is {index + 1} out of {len(sky_regions)}..')  #

                    contained_ids = sky_region.contains(SkyCoord(df['ra'], df['dec'], unit='deg'), dummy_wcs)
                    n_sources = np.sum(contained_ids)

                    if n_sources > 0:
                        df_contained = df[contained_ids]
                        outfilename = args.input_dir / 'COSMOS' / f'{args.survey}_objects_in_{label_text}.fits'
                        Table.from_pandas(df_contained).write(outfilename, overwrite=True)
                        print(f'Saved contained dataframe as {outfilename}')

                        df_contained['radius'] = 1  # arcseconds # choice of how big a circle to draw around each RA,DEC combo found in the zCOSMOS catalog
                        df_to_region(df_contained, output_dir / f'{args.survey}_objects_in_{label_text}.reg', shape='circle', label='COSMOS2020', label_ra=150.47, label_dec=2.55, color='cyan')  # making the region files
                        print(f'Region {label_text} contains {n_sources} {args.survey} sources') #


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
