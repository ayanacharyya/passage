'''
    Filename: read_line_finding_catalog.py
    Notes: Reads in Kalina's (before reconciliation) line finding catalog file
    Author : Ayan
    Created: 02-05-25
    Example: run read_line_finding_catalog.py --field Par028
'''
from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_number(item):
    '''
    Extracts numbers (floats or ints) from strings
    Returns list of numbers
    '''
    if item == 0 or item == 1 or item == -1: return int(item)
    else: return float(item)

# --------------------------------------------------------------------------------------------------------------------
def read_line_finding_catalog(filename):
    '''
    Reads the catalog produced by Kalina's line finding code (for a given user, before reconciliation)
    Returns a datafrmae
    '''
    outfilename = filename.parent / Path(str(filename.stem) + '.csv')

    if not os.path.exists(outfilename) or args.clobber:
        print(f'Not found {outfilename}, so making a new one..')
        f = open(filename, 'r')
        lines = f.readlines()

        # ----------to get all column labels-----------
        column_labels, counter = [], 0
        for line in lines:
            if line[0] == '#':
                start_index = line.find(next(filter(str.isalpha, line)))
                end_index = line.find('\n') - 1
                label = line[start_index: end_index]
                column_labels.append(label)
                counter += 1
            else:
                break
        
        # ----------to insert a new column label at a specific location-----------
        insert_idx = np.where(np.array(column_labels) == 'fwhm_g141_error')[0][0]
        column_labels.insert(insert_idx, 'comp_fit_flag')

        # ----------making new dataframe-----------
        df = pd.DataFrame()
        for line in lines[counter:]: # looping over each line in file
            if line[0] == '#': continue
            this_row = []
            for item in line.split(): # looping over each item in each line
                if 'e+' in item or 'e-' in item: # exponential notations require special treatment...
                    if len(item.split('e')[1]) > 3: # ...if they have an additional digit at the end that was supposed to be a different column
                        this_row.append(get_number(item[:-1]))
                        this_row.append(get_number(item[-1]))
                    else:
                        this_row.append(get_number(item))
                elif 'False' in item: # items with False or True combined with digits in the same string require special handling to split them
                    items = item.split('False')
                    this_row.append(get_number(items[0]))
                    this_row.append(False)
                elif 'True' in item:
                    items = item.split('True')
                    this_row.append(get_number(items[0]))
                    this_row.append(True)
                else:
                    this_row.append(get_number(item))
            this_df = pd.DataFrame(dict(zip(column_labels, this_row)), index=[0])
            df = pd.concat([df, this_df])

        # ------writing out the dataframe---------
        df['objid'] = df['objid'].astype(int)
        df.to_csv(outfilename, index=None)
        print(f'Saved as new table in {outfilename}')
    else:
        print(f'Reading existing csv linelist from {outfilename}')
    
    df = pd.read_csv(outfilename)

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -------reading in catalog-------
    user = 'scarlata'
    catalog_filename = args.input_dir / 'data' / args.field / 'line_finding' / f'{args.field}lines_catalog_{user}.dat'
    df = read_line_finding_catalog(catalog_filename)

    # -----comparing with PJW linelist-------
    dfp = Table(fits.open(str(catalog_filename).replace('.dat', '.fits'))[1].data).to_pandas()
    dfp['objid'] = dfp['objid'].astype(int)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
