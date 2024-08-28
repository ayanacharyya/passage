'''
    Filename: make_combined_animation_from_gdrive.py
    Notes: Batch script to download reduced PASSAGE data from Goggle Drive, then run make_diagnostic_maps.py and combine_diagnostics_and_extractions.py
    and animate_png.py to finally produce an animation of combined data reduction adn diagnostic frames
    Author : Ayan
    Created: 28-08-24
    Example: run make_combined_animation_from_gdrive.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_combined_animation_from_gdrive.py
'''

from header import *
from util import *

from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from google_auth_oauthlib.flow import InstalledAppFlow
import io

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def query_google_drive_folder(folder_id):
    '''
    Queries google drive folder and returns file/folder list
    Based on top answer in https://stackoverflow.com/questions/76485003/how-do-i-download-all-files-from-a-google-drive-folder-with-more-than-50-files
    '''
    # Define the scopes
    SCOPES = ['https://www.googleapis.com/auth/drive.readonly']

    # Obtain your Google credentials
    def get_credentials():
        flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
        creds = flow.run_local_server(port=0)
        return creds

    # Build the downloader
    creds = get_credentials()
    drive_downloader = build('drive', 'v3', credentials=creds)

    query = f"'{folder_id}' in parents"  # this works  ref https://stackoverflow.com/q/73119251/248616

    results = drive_downloader.files().list(q=query, pageSize=1000).execute()
    items = results.get('files', [])

    print(f'Found total {len(items)} files/folders in the folder.')

    return items, drive_downloader

# --------------------------------------------------------------------------------------------------------------------
def download_folder_from_google_drive(folder_id, destination_folder):
    '''
    Downloads the given file IDs from google drive
    Based on top answer in https://stackoverflow.com/questions/76485003/how-do-i-download-all-files-from-a-google-drive-folder-with-more-than-50-files
    '''
    destination_folder.mkdir(parents=True, exist_ok=True)
    items, drive_downloader = query_google_drive_folder(folder_id)

    for item in items:
        request = drive_downloader.files().get_media(fileId=item['id'])
        f = io.FileIO(destination_folder / item['name'], 'wb')
        downloader = MediaIoBaseDownload(f, request)
        done = False
        while not done:
            status, done = downloader.next_chunk()
            print(f'Downloaded {item["name"]}: {int(status.progress() * 100)}\%', end='\r')
        print('\n')

    print(f'Downloaded {len(items)} files from the folder.')


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    field_list = passage_fields_in_cosmos
    passage_url_id = '1oE1V76aeDlR8vxZNefyzk3yp0eO1rkHw'

    # ------determining directories and global variables---------
    args.plot_radial_profiles, args.only_seg, args.snr_cut = True, True, 3 #
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    snr_text = f'_snr{args.snr_cut:.1f}' if args.snr_cut is not None else ''
    description_text1 = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'
    description_text2 = f'diagnostics_and_extractions'

    # --------query passage folder and see which fields available------------------
    items, _ = query_google_drive_folder(passage_url_id)
    field_url_dict = {item['name']:item['id'] for item in items if item['name'].startswith('Par') and item ['mimeType'].endswith('.folder')}
    field_list = field_url_dict.keys()
    print(f'Found {len(field_list)} PASSAGE fields in the folder.')

    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{field:03}'
        short_field_name = f'Par{field}'
        print(f'\nCommencing field {args.field} which is {index+1} of {len(field_list)}..')

        # ----------determining filenames etc. to check for presence----------
        products_path = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / args.field
        file_to_check_for = f'{args.field}__{description_text2}_anim.mp4'

        if os.path.exists(output_dir / file_to_check_for):
            print(f'{file_to_check_for} already present for {args.field}, so skipping this field.')

        else:
            # ------------download the files------------------
            if os.path.exists(products_path / 'spec1D'):
                print(f'Downloads already present, so moving to the next step.')
            else:
                print(f'Downloading form google drive..')
                folder_id = field_url_dict[field]
                download_folder_from_google_drive(folder_id, products_path)

            # ------------unzip the downloaded files------------------
            zipped_files = glob.glob(str(products_path) + '/*.gz')
            if len(zipped_files) == 0:
                print(f'All files already unzipped, so moving to the next step.')
            else:
                for ind, thisfile in zipped_files:
                    print(f'Unzipping {ind + 1} of {len(zipped_files)} zipped files..')
                    shutil.unpack_archive(thisfile, products_path)
                    os.remove(thisfile) #remove zipped files after unzipping

            # ------------rename the downloaded files------------------
            files_to_rename = glob.glob(str(products_path) + '/' + short_field_name + '*')
            if len(files_to_rename) == 0:
                print(f'All files/folders already have the correct nomenclature, so moving to the next step.')
            else:
                print(f'Renaming files/folders..')
                for thisfile in files_to_rename:
                    os.rename(thisfile, thisfile.replace(short_field_name, args.field))

            # ------------run make_diagnostic_maps.py------------------
            diag_results_file = output_dir / f'{args.field}_all_diag_results.txt'
            if os.path.exists(diag_results_file) and not args.clobber:
                print(f'Diagnostic results file already present, so moving to the next step.')
            else:
                print(f'Running make_diagnostic_maps.py..')
                dummy = subprocess.check_output([f'python make_diagnostic_maps.py --field {args.field} --do_all_obj --plot_radial_profiles --only_seg --snr_cut 3 --write_file'], shell=True)

            # ------------run combine_diagnostics_and_extractions.py------------------
            diagnostic_img_files = glob.glob(str(output_dir / f'{description_text1}') + f'/{args.field}_*_{description_text1}.png')
            extraction_img_files = glob.glob(str(output_dir / f'{description_text2}') + f'/{args.field}_*_{description_text2}.png')
            if len(extraction_img_files) == len(diagnostic_img_files):
                print(f'All extraction images already present, so moving to the next step.')
            else:
                print(f'Running combine_diagnostics_and_extractions.py..')
                dummy = subprocess.check_output([f'python combine_diagnostics_and_extractions.py --field {args.field} --keep --do_all_obj'], shell=True)

            # ------------make the final animation with the combined images------------------
            file_to_move = output_dir / f'{description_text2}' / f'{args.field}__{description_text2}_anim.mp4'
            if os.path.exists(file_to_move):
                print(f'Animation already present, so moving to the next step.')
            else:
                print(f'Running animate_png.py..')
                dummy = subprocess.check_output([f'python /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath {output_dir}/{description_text2} --rootname {args.field}_*_{description_text2}.png - -delay 0.1'], shell=True)

            # ------------moving the final animation------------------
            print(f'Moving the animation file..')
            dummy = shutil.move(file_to_move, file_to_check_for)


        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
