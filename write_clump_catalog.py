'''
    Filename: write_clump_catalog.py
    Notes: Writes out combined clump catalog from several individual catalogs (for Ramona)
    Author : RA, AA
    Created: 06-09-24
    Example: run write_clump_catalog.py
'''
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import os
from astropy.table import Table

import warnings
warnings.filterwarnings('ignore')

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def fill_na(df, columns, fill_with='N/A'):
    '''
    Function to fill given columns of a given dataframe with fill_na
    Returns dataframe
    '''
    if columns == 'all': columns = df.columns()
    for column in columns: df[column] = fill_with
    return df

# --------------------------------------------------------------------------------------------------------------------
def copy_columns(target_df, source_df, columns):
    '''
    Function to copy given columns of a source dataframe to target dataframe
    Returns target dataframe
    '''
    if columns == 'all': columns = source_df.columns()
    for column in columns: target_df[column] = source_df[column]
    return target_df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # -------------modify following if needed------------------------------------
    make_full = True # keeping this on will make the fully combined df
    do_survivability = False # keeping this on will add survivability columns
    do_strictly = False # keeping this on will check if lengths of two dataframes are same before merging them, and therefore ensure it is being done correctly
                        # but keeping it False will make sure the whole code runs, but the result may not be reliable

    input_path = '/Users/acharyya/Downloads/tables/'
    output_path = input_path + 'outputs/'
    os.makedirs(output_path, exist_ok=True)

    halos = ['002392', '002878', '004123', '005016', '005036', '008508']
    levels = ['1', '2', '3', '4']

    cols_to_grab_from_metfrac = ['envmets', 'metalratio']
    cols_to_grab_from_survivability = ['Tcool_edge', 'Tshear', 'survivability']
    cols_to_grab_from_tree = ['sanitycheckmasses', 'avgdensity', 'avgtemp', 'cid']

    full_outfilename = output_path + 'fullclumpstable_master2.h5'
    basic_outfilename = output_path + 'basicclumpstable_master2.h5'
    metfrac_outfilename = output_path + 'metfracclumpstable_master2.h5'
    survivability_outfilename = output_path + 'survivabilityclumpstable_master2.h5'

    # -------------should not need to modify below------------------------------------
    combined_all_df = pd.DataFrame()
    combined_basic_df = pd.DataFrame()
    combined_metfrac_df = pd.DataFrame()
    combined_survivability_df = pd.DataFrame()

    combined_basic_df.to_hdf(basic_outfilename, key='combined_basic_df', mode='w', format='table')
    combined_metfrac_df.to_hdf(metfrac_outfilename, key='combined_metfrac_df', mode='w', format='table')
    if do_survivability: combined_survivability_df.to_hdf(survivability_outfilename, key='combined_survivability_df', mode='w', format='table')
    if make_full: combined_all_df.to_hdf(full_outfilename, key='combined_all_df', mode='w', format='table')

    # ----------reading in the treefile----------------
    tree_infilename = input_path + 'treetable.csv'
    if os.path.exists(tree_infilename):
        print(f'Reading master tree file {tree_infilename}..')
        master_tree_df = pd.read_csv(tree_infilename)
        master_tree_df = master_tree_df.rename(columns={'fullhalo': 'halo', 'fulllevel': 'level'})

    # ----------looping over levels and halos----------------
    for index in range(len(halos) * len(levels)):
        thishalo = halos[int(index / len(levels))]
        thislevel = levels[index % len(levels)]
        print(f'\nStarting iteration {index+1} out of {len(halos) * len(levels)}: level {thislevel} of halo {thishalo}..')

        # -------------do halo stuff------------------------------------
        this_basic_infilename = input_path + f'shell{thislevel}/halo_{thishalo}_nref11c_nref9f_RD0027_RD0027_box1_clump_measurements.fits'
        if os.path.exists(this_basic_infilename):
            print(f'Reading halo file {this_basic_infilename}..')
            this_basic_df = Table.read(this_basic_infilename, format='fits').to_pandas()

            this_basic_df['halo'] = thishalo
            this_basic_df['level'] = thislevel

            combined_basic_df = pd.concat([combined_basic_df, this_basic_df])  # appending to master df
            this_basic_df.to_hdf(basic_outfilename, mode='a', append=True, key='combined_basic_df', format='table')  # also saving to file, just in case

            # -------------do metfrac stuff------------------------------------
            this_metfrac_infilename = input_path + f'metfracs/halo_{thishalo}_nref11c_nref9f_RD0027_level{thislevel}_clump_mets.fits'
            if os.path.exists(this_metfrac_infilename):
                print(f'Reading metfrac file {this_metfrac_infilename}..')
                this_metfrac_df = Table.read(this_metfrac_infilename, format='fits').to_pandas()
                this_metfrac_df['metalratio'] = this_metfrac_df['metallicity'] - this_metfrac_df['envmets']  # adding new column temporarily; subtracting to get 'ratio', because log space?

                if not do_strictly or len(this_basic_df) == len(this_metfrac_df):
                    this_basic_df = copy_columns(this_basic_df, this_metfrac_df, cols_to_grab_from_metfrac)
                else:
                    print(
                    f'Concatenating N/A instead of metfrac stuff because df lengths do not match: {len(this_basic_df)} and {len(this_metfrac_df)}, set do_strictly = False to bypass this check.')
                    this_basic_df = fill_na(this_basic_df, cols_to_grab_from_metfrac)
            else:
                print(f'Could not find {this_metfrac_infilename}, so concatenating N/A instead of metfrac stuff.')
                this_basic_df = fill_na(this_basic_df, cols_to_grab_from_metfrac)

            combined_metfrac_df = pd.concat([combined_metfrac_df, this_basic_df])  # appending to master df
            this_basic_df.to_hdf(metfrac_outfilename, key='combined_metfrac_df', mode='a', append=True, format='table')  # also saving to file, just in case

            # -------------do survivability stuff if needed------------------------------------
            if do_survivability:
                # combined_survivability_df = pd.concat([combined_metfrac_df, this_basic_df])
                this_survivability_infilename = input_path + f'survivabilityfiles/halo_{thishalo}_nref11c_nref9f_RD0027_RD0027_level{thislevel}_clump_survivability.fits'
                if os.path.exists(this_survivability_infilename):
                    print(f'Reading survivability file {this_survivability_infilename}..')
                    this_survivability_df = Table.read(this_survivability_infilename, format='fits').to_pandas()

                    if not do_strictly or len(this_basic_df) == len(this_survivability_df):
                        this_basic_df = copy_columns(this_basic_df, this_survivability_df, cols_to_grab_from_survivability)
                    else:
                        print(f'Concatenating N/A instead of survivability stuff because df lengths do not match: {len(this_basic_df)} and {len(this_survivability_df)}, set do_strictly = False to bypass this check.')
                        this_basic_df = fill_na(this_basic_df, cols_to_grab_from_survivability)
                else:
                    print(
                    f'Could not find {this_survivability_infilename}, so concatenating N/A instead of survivability stuff.')
                    this_basic_df = fill_na(this_basic_df, cols_to_grab_from_survivability)

                combined_survivability_df = pd.concat([combined_survivability_df, this_basic_df])  # appending to master df
                this_basic_df.to_hdf(survivability_outfilename, key='combined_survivability_df', mode='a', append=True, format='table')  # also saving to file, just

            # -------------do tree file stuff if needed------------------------------------
            if make_full:
                if os.path.exists(tree_infilename):
                    print(f'Taking subset of tree file {tree_infilename}..')
                    this_tree_df = master_tree_df[(master_tree_df['halo'] == thishalo) & (master_tree_df['level'] == thislevel)]

                    if not do_strictly or len(this_basic_df) == len(this_tree_df):
                        this_basic_df = copy_columns(this_basic_df, this_tree_df, cols_to_grab_from_tree)
                    else:
                        print(f'Concatenating N/A instead of tree stuff because df lengths do not match: {len(this_basic_df)} and {len(this_tree_df)}, set do_strictly = False to bypass this check.')
                        this_basic_df = fill_na(this_basic_df, cols_to_grab_from_tree)
                else:
                    print(f'Could not find {tree_infilename}, so concatenating N/A instead of tree stuff.')
                    this_basic_df = fill_na(this_basic_df, cols_to_grab_from_tree)

                combined_all_df = pd.concat([combined_all_df, this_basic_df])  # appending to master df
                this_basic_df.to_hdf(full_outfilename, key='combined_all_df', mode='a', append=True, format='table')  # also saving to file, just

        else:
            print(f'Could not find {this_basic_infilename}, so skipping this halo/level.')

    # ----------------------adding columns to combined basic df-------------------------
    combined_basic_df['HIfraction'] = np.log10(combined_basic_df['HImasses'] / (combined_basic_df['clumpmasses']))
    combined_basic_df['pressureratio'] = np.log10(combined_basic_df['pressures'] / combined_basic_df['envpressures'])
    combined_basic_df['velocityratio'] = np.log10(combined_basic_df['radialvelocities'] / combined_basic_df['envvels'])

    # ----------------------saving combined basic df-------------------------
    combined_basic_df.to_hdf(output_path + 'basicclumpstable_master_done.h5', key='combined_basic_df', mode='w')

    # ----------------------saving combined metfrac df-------------------------
    combined_metfrac_df.to_hdf(output_path + 'fullclumpstable_master_withmet.h5', key='combined_metfrac_df', mode='w')

    # ----------------------saving combined sustainability df-------------------------
    if do_survivability:
        combined_survivability_df.to_hdf(output_path + 'fullclumpstable_master_survivability.h5', key='combined_metfrac_df', mode='w')

    # ----------------------adding columns to combined all df-------------------------
    if make_full:
        combined_all_df['HIfraction'] = np.log10(combined_all_df['HImasses'] / (combined_all_df['clumpmasses']))
        combined_all_df['pressureratio'] = np.log10(combined_all_df['pressures'] / combined_all_df['envpressures'])
        combined_all_df['densityratio'] = np.log10(combined_all_df['avgdensity'] / combined_all_df['envdensities'])
        combined_all_df['tempratio'] = np.log10(combined_all_df['avgtemp'] / combined_all_df['envtemps'])
        combined_all_df['velocityratio'] = np.log10(combined_all_df['radialvelocities'] / combined_all_df['envvels'])

        # ----------------------saving combined all df-------------------------
        combined_all_df.to_hdf(output_path + 'fullclumpstable_test_done.h5', key='combined_all_df', mode='w')

print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
