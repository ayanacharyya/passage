'''
    Filename: display_images.py
    Notes: Displays existing images by reading them from saved pngs, based on an input array of IDs (which is used to determine the image filenames to be read in)
    Author : Ayan
    Created: 15-07-24
    Example: run display_images.py
'''

import pandas as pd
from PIL import Image
from datetime import datetime, timedelta
from pathlib import Path
import shutil

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

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

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
