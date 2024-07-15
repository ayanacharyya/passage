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
    count = 0
    for id in df['id'].values:
        try:
            im=Image.open(f'/Users/acharyya/Work/astro/passage/passage_output/Par050/{id}/Par050_{id}.line.png')
            im.show()
            count+=1
            print(id, count)
        except:
            pass

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
