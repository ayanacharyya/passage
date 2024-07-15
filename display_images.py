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
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    df=pd.read_table('/Users/acharyya/Work/astro/passage/passage_data/Par050/id_list_spatially_resolved_galaxies_Par50.txt', comment='#', names=['id'])
    image_path = '/Users/acharyya/Work/astro/passage/passage_output/Par050/'

    count = 0
    for this_id in df['id'].values:
        try:
            this_id = f'{int(this_id.strip()):05d}'
            im=Image.open(image_path + f'{this_id}/Par050_{this_id}.line.png')
            im.show()
            count += 1
            print(this_id, count)
        except (FileNotFoundError, ValueError) as e:
            pass

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
