'''
    Filename: write_clump_catalog.py
    Notes: Writes out combined clump catalog from several individual catalogs
    Author : RA, AA
    Created: 06-09-24
    Example: run write_clump_catalog.py
'''
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import os

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def fill_na(df, columns, fill_with='N/A'):
    '''
    Function to fill given columns of a given dataframe with fill_na
    Returns dataframe
    '''
    for column in columns: df[column] = fill_with
    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    makefull = False
    surv = False
    do_strictly = False

    input_path = '/Users/acharyya/Downloads/tables/'
    output_path = input_path + 'outputs/'
    os.makedirs(output_path, exist_ok=True)
    
    full_outfilename = output_path + 'fullclumpstable_master2.h5'
    basic_outfilename = output_path + 'basicclumpstable_master2.h5'
    metfrac_outfilename = output_path + 'metfracclumpstable_master2.h5'
    survivability_outfilename = output_path + 'survivabilityclumpstable_master2.h5'
    
    fulldf = pd.DataFrame()
    basicdf = pd.DataFrame()
    basic_metfrac_df = pd.DataFrame()
    basic_survivability_df = pd.DataFrame()
    
    halos = ['002392', '002878', '004123', '005016', '005036', '008508']
    levels = ['1', '2', '3', '4']

    if makefull: fulldf.to_hdf(full_outfilename, key='fulldf', mode='w', format='table')
    basicdf.to_hdf(basic_outfilename, key='basicdf', mode='w', format='table')
    basic_metfrac_df.to_hdf(metfrac_outfilename, key='basic_metfrac_df', mode='w', format='table')
    if surv: basic_survivability_df.to_hdf(survivability_outfilename, key='basic_survivability_df', mode='w', format='table')

    # read in the treefile here

    for index in range(len(halos) * len(levels)):
        thishalo = halos[int(index / len(levels))]
        thislevel = levels[index % len(levels)]
        print(f'\nStarting level {thislevel} of halo {thishalo}..')

        # -------------do halo stuff------------------------------------
        halo_filename = input_path + f'shell{thislevel}/halo_{thishalo}_nref11c_nref9f_RD0027_RD0027_box1_clump_measurements.fits'
        if os.path.exists(halo_filename):
            print(f'Found halo file {halo_filename}..')
            halo_level_df = Table.read(thisfile, format='fits').to_pandas()

            halo_level_df['halo'] = thishalo
            halo_level_df['level'] = thislevel

            basicdf = pd.concat([basicdf, halo_level_df])
            halo_level_df.to_hdf(basic_filename, mode='a', append=True, key='basicdf', format='table')

            # -------------do metfrac stuff------------------------------------
            metfrac_filename = input_path + f'metfracs/halo_{thishalo}_nref11c_nref9f_RD0027_level{thislevel}_clump_mets.fits'
            if os.path.exists(metfrac_filename):
                print(f'Found halo file {metfrac_filename}..')
                metfrac_df = Table.read(metfracfilename, format='fits').to_pandas()

                if not do_strictly or len(halo_level_df) == len(metfrac_df):
                    halo_level_df['envmets'] = metfrac_df['envmets']
                    halo_level_df['metalratio'] = halo_level_df['metallicity'] - halo_level_df['envmets']
                else:
                    print(f'NOT concatenating metfrac stuff because df lengths do not match: {len(halo_level_df)} and {len(metfrac_df)}, set do_strictly = False to bypass this check.')
                    halo_df = fill_na(halo_df, ['envmets', 'metalratio'])
            else:
                print(f'Could not find {metfrac_filename}, so not concatenating metfrac stuff.')
                halo_df = fill_na(halo_df, ['envmets', 'metalratio'])

            basic_metfrac_df = pd.concat([basic_metfrac_df, halo_level_df])
            halo_level_df.to_hdf(metfrac_filename, key='basic_metfrac_df', mode='a', append=True, format='table')

            # -------------do survivability stuff if needed------------------------------------
            if surv:
                #basic_survivability_df = pd.concat([basic_metfrac_df, halo_level_df])
                survivabilityfilename = input_path + f'survivabilityfiles/halo_{thishalo}_nref11c_nref9f_RD0027_RD0027_level{thislevel}_clump_survivability.fits'
                if os.path.exists(survivabilityfilename):
                    survdat = Table.read(survivabilityfilename, format='fits')
                    survdf = survdat.to_pd()
                    halo_level_df['Tcool_edge'] = survdf['Tcool_edge']
                    halo_level_df['Tshear'] = survdf['Tshear']
                    halo_level_df['survivability'] = survdf['survivability']

                else:
                    print('NOPE')
                    halo_level_df['Tcool_edge'] = 'N/A'
                    halo_level_df['Tshear'] = 'N/A'
                    halo_level_df['survivability'] = 'N/A'
                basic_survivability_df = pd.concat([basic_survivability_df, halo_level_df])
                halo_level_df.to_hdf(survivability_filename, key='basic_survivability_df', mode='a')

        else:
            print(f'Could not find {halofilename}, so skipping this halo/level.')
            continue

        #### adding tree files for environment differences:
        if makefull:
            if os.path.exists(halofilename) and os.path.exists(treefilename):
                thisfile = halofilename
                dat = Table.read(thisfile, format='fits')
                halo_level_df = dat.to_pd()
                for leveli in levellist:
                    if leveli in thisfile: halo_level_df['level'] = leveli
                for haloi in halos:
                    if haloi in thisfile: halo_level_df['halo'] = haloi

                t = yt.load(treefilename)
                sanitycheckmasses = []
                avgdensity = []
                avgtemp = []
                cid = []

                for c in t.leaves:
                    sanitycheckmasses.append(c['grid', 'cell_mass'].sum().in_units("Msun"))
                    avgdensity.append(c['grid', 'density'].mean())
                    avgtemp.append(c['grid', 'temperature'].mean())
                    cid.append(c['clump', 'clump_id'])

                halo_level_df['sanitycheckmasses'] = np.array(sanitycheckmasses)
                halo_level_df['avgdensity'] = np.array(avgdensity)
                halo_level_df['avgtemp'] = np.array(avgtemp)
                halo_level_df['cid'] = np.array(cid)

                fulldf = pd.concat([fulldf, halo_level_df])
                halo_level_df.to_hdf(full_filename, key='fulldf', mode='a')

            else:
                print('Could not find combination of ' + halofilename + ' and ' + treefilename)
                continue

    if makefull: fulldf['HIfraction'] = np.log10(fulldf['HImasses'] / (fulldf['clumpmasses']))
    if makefull: fulldf['pressureratio'] = np.log10(fulldf['pressures'] / fulldf['envpressures'])
    if makefull: fulldf['densityratio'] = np.log10(fulldf['avgdensity'] / fulldf['envdensities'])
    if makefull: fulldf['tempratio'] = np.log10(fulldf['avgtemp'] / fulldf['envtemps'])
    if makefull: fulldf['velocityratio'] = np.log10(fulldf['radialvelocities'] / fulldf['envvels'])

    basicdf['HIfraction'] = np.log10(basicdf['HImasses'] / (basicdf['clumpmasses']))
    basicdf['pressureratio'] = np.log10(basicdf['pressures'] / basicdf['envpressures'])
    # basicdf['densityratio'] = np.log10(basicdf['avgdensity']/basicdf['envdensities'])
    # basicdf['tempratio'] = np.log10(basicdf['avgtemp']/basicdf['envtemps'])
    basicdf['velocityratio'] = np.log10(basicdf['radialvelocities'] / basicdf['envvels'])

    if makefull: fulldf.to_hdf('/Users/ramonaaugustin/WORK/FOGGIE/Outputs/clumptables/fullclumpstable_test_done.h5', key='fulldf', mode='w')
    basicdf.to_hdf('/Users/ramonaaugustin/WORK/FOGGIE/Outputs/clumptables/basicclumpstable_master_done.h5', key='basicdf', mode='w')

    basic_metfrac_df.to_hdf('/Users/ramonaaugustin/WORK/FOGGIE/Outputs/clumptables/fullclumpstable_master_withmet.h5', key='basic_metfrac_df', mode='w')
    if surv: basic_survivability_df.to_hdf('/Users/ramonaaugustin/WORK/FOGGIE/Outputs/clumptables/fullclumpstable_master_survivability.h5', key='basic_metfrac_df', mode='w')

print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
