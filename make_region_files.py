'''
    Filename: make_region_files.py
    Notes: Generate region files of all PASSAGE field footprints
    Author : Ayan
    Created: 26-06-24
    Last modified: 26-06-24
    Example: run make_region_files.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_region_files.py --survey zcosmos
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
def get_passage_footprints(filename):
    '''
    Reads and returns the footprint list for PASSAGE fields
    '''
    NIRISS_WFSS_fov = [133, 133] # FoV in arcseconds

    df = pd.read_csv(filename)
    df = df[['Par#', 'Target RA (J2000)', 'Target Dec (J2000)']].dropna().reset_index(drop=True)
    df = df.rename(columns={'Par#':'id', 'Target RA (J2000)':'ra', 'Target Dec (J2000)':'dec'})

    df['ra_width'] = NIRISS_WFSS_fov[0] # arcseconds
    df['dec_width'] = NIRISS_WFSS_fov[1]
    df['pa'] = 0 # dummy value; in degrees

    return df # df only contains id, ra, dec, ra_width, dec_width, pa (orienation) for each field

# -------------------------------------------------------------------------------------------------------
def get_zCOSMOS_footprints(filename):
    '''
    Reads and returns the footprint list for all zCOSMOS galaxies (with spectra)
    '''

    data = fits.open(filename)
    table = Table(data[1].data)
    df = table.to_pandas()

    df = df.rename(columns={'OBJECT_ID':'id', 'RAJ2000':'ra', 'DEJ2000':'dec'})
    df = df[['ra', 'dec', 'id']]
    df['radius'] = 1 # arcseconds # choice of how big a circle to draw around each RA,DEC combo found in the zCOSMOS catalog

    return df # df only contains id, ra, dec, radius, for each object

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    if args.survey == 'passage':
        infilename = 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'
        outfilename = os.path.splitext(infilename.replace(' ', '-'))[0] + '.reg'

        df = get_passage_footprints(args.input_dir / infilename) # reading the passage footprints
        df_to_region(df, args.input_dir / 'footprints/region_files' / outfilename, shape='box', label='JWST-PASSAGE', label_ra=150.47, label_dec=2.55, color='red')  # making the region files

    elif args.survey == 'zcosmos':
        infilename = 'zCOSMOS_VIMOS_BRIGHT_DR3_CATALOGUE.fits'
        outfilename = os.path.splitext(infilename.replace(' ', '-'))[0] + '.reg'

        df = get_zCOSMOS_footprints(args.input_dir / 'zCOSMOS-DR3' / infilename) # reading the zCOSMOS footprints
        df_to_region(df, args.input_dir / 'footprints/region_files' / outfilename, shape='circle', label='zCOSMOS', label_ra=150.47, label_dec=2.55, color='cyan') # making the region files

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
