'''
    Filename: download_raw_data.py
    Notes: Download the raw *rate.fits files, for a given field, from MAST
           This script is based on Vihang Mehta's PASSAGE_Pipeline/passagepipe/utils.py
    Author : Ayan
    Created: 10-07-24
    Example: run download_raw_data.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par008
             run download_raw_data.py --field Par008
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def download_from_mast_old(args, prop_id='1571'):
    '''
    Fetch the individual visits associated to a field
    Query MAST and then trim the full MAST-queried table to a subset of specified observations, and download the trimmed set
    '''
    # ---------read google sheet-----------------
    sheetPath = '15BNOXSjA8xBRO6la8jf5OCID81QZR1Y7thfZesgOsV4'
    sheetName = 'Cy1+Executed'
    df = pd.read_csv(f'https://docs.google.com/spreadsheets/d/{sheetPath}/gviz/tq?tqx=out:csv&sheet={sheetName}')

    # -------determine obs ids from sheet--------------
    all_par_ids = np.where(~df['Par#'].isnull())[0]
    this_par_id = np.where(df['Par#']==args.field)[0][0]
    next_par_id = all_par_ids[np.where(all_par_ids == this_par_id)[0][0]+1]
    obs_ids = df.iloc[this_par_id : next_par_id]['Obs ID']
    print(f'Found these obs IDs for field {args.field}', obs_ids.values)

    # --------build MAST query---------------
    query.DEFAULT_QUERY['project'] = ['JWST']
    query.DEFAULT_QUERY['obs_collection'] = ['JWST']

    queryList = query.run_query(box=None, proposal_id=[prop_id], base_query=query.DEFAULT_QUERY)
    queryList['obs_id_num'] = pd.Series(queryList['obs_id']).map(lambda x: int(x.split('_')[1][1:]) if 'niriss' in x else int(x.split('_')[0][7:-3]))

    # --------trim MAST query---------------
    cond = np.in1d(queryList['obs_id_num'], obs_ids)
    queryList = queryList[cond]

    # --------modify MAST query, to download only the rate files---------------
    ratefile_list = []
    for index, url in enumerate(queryList['dataURL']):
        new_url = url.replace('_rateints', '_rate').replace('_cal', '_rate')
        ratefile_list.append(new_url)
    queryList['dataURL'] = ratefile_list

    # --------download from MAST---------------
    print(f"About to download the following {len(queryList['obs_id'])} rate files..", queryList['obs_id'])
    mastutils.download_from_mast(queryList, path=args.download_dir, force_rate=True, overwrite=args.clobber)

# --------------------------------------------------------------------------------------------------------------------
def download_from_mast(args, prop_id='1571'):
    '''
    Fetch the individual visits associated to a field
    Query MAST and then trim the full MAST-queried table to a subset of specified observations, and download the trimmed set
    '''

    # ---------read google sheet-----------------
    sheetPath = '15BNOXSjA8xBRO6la8jf5OCID81QZR1Y7thfZesgOsV4'
    sheetName = 'Cy1+Executed'
    df = pd.read_csv(f'https://docs.google.com/spreadsheets/d/{sheetPath}/gviz/tq?tqx=out:csv&sheet={sheetName}')

    # -------determine obs ids from sheet--------------
    all_par_ids = np.where(~df['Par#'].isnull())[0]
    this_par_id = np.where(df['Par#']==args.field)[0][0]
    next_par_id = all_par_ids[np.where(all_par_ids == this_par_id)[0][0]+1]
    obs_ids = df.iloc[this_par_id : next_par_id]['Obs ID']
    print(f'Found these obs IDs for field {args.field}', obs_ids.values)

    # --------build MAST query---------------
    query.DEFAULT_QUERY['project'] = ['JWST']
    query.DEFAULT_QUERY['obs_collection'] = ['JWST']

    queryList = query.run_query(box=None, proposal_id=[prop_id], base_query=query.DEFAULT_QUERY)
    subqueryList = Observations.get_product_list(queryList)

    cond = (subqueryList['calib_level'] == 1) & (subqueryList['productType'] == 'SCIENCE') & (subqueryList['productSubGroupDescription'] == 'UNCAL')
    uncalList = subqueryList[cond]
    uncalList['obs_id_num'] = pd.Series(uncalList['obs_id']).map(lambda x: int(x.split('_')[1][1:]) if 'niriss' in x else int(x.split('_')[0][7:-3]))

    # --------trim MAST query---------------
    cond = np.in1d(uncalList['obs_id_num'], obs_ids)
    uncalList = uncalList[cond]

    # --------modify MAST query, to download only the rate files---------------
    ratefile_list = []
    for index, url in enumerate(uncalList['dataURL']):
        new_url = url.replace('_rateints', '_rate').replace('_cal', '_rate')
        ratefile_list.append(new_url)
    uncalList['dataURL'] = ratefile_list

    # --------download from MAST---------------
    print(f"About to download the following {len(uncalList['obs_id'])} rate files..", uncalList['obs_id'])
    mastutils.download_from_mast(uncalList, path=args.download_dir, force_rate=True, overwrite=args.clobber)

# --------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()

    args.download_dir = args.input_dir / args.field / 'MAST'
    args.download_dir.mkdir(parents=True, exist_ok=True)
    download_from_mast(args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
