'''
    Filename: make_region_files.py
    Notes: Generate region files of all PASSAGE field footprints
    Author : Ayan
    Created: 26-06-24
    Last modified: 26-06-24
    Example: run make_region_files.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_region_files.py
'''

from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def df_to_region(df, outfilename, color='cyan'):
    '''
    Converts the footprint list into a region file and saves the file
    '''
    with open(outfilename, 'w') as f: f.write('fk5\nglobal color=' + color + '\n')

    df['shape'] = 'box'
    df['ra_width'] = df['ra_width'].astype(str) + '\"' # to specify arcseconds for the width
    df['dec_width'] = df['dec_width'].astype(str) + '\"' # to specify arcseconds for the width
    df['id'] = '# text={' + df['id'].str[3:].astype(int).astype(str) + '}' # to specify field id names as text
    df = df[['shape', 'ra', 'dec', 'ra_width', 'dec_width', 'pa', 'id']]

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


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    infilename = 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'
    outfilename = os.path.splitext(infilename.replace(' ', '-'))[0] + '.reg'

    df = get_passage_footprints(args.input_dir / infilename) # reading the passage footprints
    df_to_region(df, args.input_dir / 'footprints/region_files' / outfilename) # making the region files

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
