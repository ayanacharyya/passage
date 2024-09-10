'''
    Filename: display_images.py
    Notes: Displays existing images by reading them from saved pngs, based on an input array of IDs (which is used to determine the image filenames to be read in)
    Author : Ayan
    Created: 15-07-24
    Example: run display_images.py
    Afterwards, to make animations, one can do something like: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par061/diagnostics_and_extractions/ --rootname Par061_*_diagnostics_and_extractions.png --delay 0.1
'''

import pandas as pd
from PIL import Image
from datetime import datetime, timedelta
from pathlib import Path
import shutil
import os, subprocess, re

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    '''
    df=pd.read_table('/Users/acharyya/Work/astro/passage/passage_data/Par050/id_list_spatially_resolved_galaxies_Par50.txt', comment='#', names=['id'])
    image_path = Path('/Users/acharyya/Work/astro/passage/passage_output/Par050/')

    quantity = 'full' # 'line' #
    new_dir = image_path / f'Malkan_allIDs_{quantity}'
    new_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for this_id in df['id'].values:
        try:
            this_id = f'{int(this_id.strip()):05d}'
            file = image_path / f'{this_id}/Par050_{this_id}.{quantity}.png'

            #im=Image.open(file)
            #im.show()
            shutil.copy(file, new_dir / '.')

            count += 1
            print(this_id, count)
        except (FileNotFoundError, ValueError) as e:
            print('Skipping', this_id)
            pass
    '''

    output_dir = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/')
    df = pd.read_csv(output_dir / 'allpar_venn_EW,PA_df.txt')

    df = df[df['field']=='Par028']
    df = df[df['Hb_EW']/df['OIII_EW'] > 0.05]

    suffix = 'diagnostics_and_extractions'
    new_dir = output_dir / f'allpar_venn_EW,PA_{suffix}'
    shutil.rmtree(new_dir)
    new_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for index in range(len(df)):
        thisfield = df.iloc[index]['field']
        thisid = df.iloc[index]['objid']
        print(f'Doing {thisfield}:{thisid}, which is {index + 1} out of {len(df)}..')
        try:
            inpath = output_dir / thisfield / suffix
            file = inpath / f'{thisfield}_{thisid:05d}_{suffix}.png'
            shutil.copy(file, new_dir / '.')
            count += 1
        except (FileNotFoundError, ValueError) as e:
            print(f'Skipping {thisid} due to {e}')
            pass

    dummy = subprocess.run(['python', '/Users/acharyya/Work/astro/ayan_codes/animate_png.py', '--inpath', f'{new_dir}', '--rootname', f'Par*_*_{suffix}.png', '--delay', '0.1'])

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
