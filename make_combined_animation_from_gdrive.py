'''
    Filename: make_combined_animation_from_gdrive.py
    Notes: Batch script to download reduced PASSAGE data from Goggle Drive, then run make_diagnostic_maps.py and combine_diagnostics_and_extractions.py
    and animate_png.py to finally produce an animation of combined data reduction adn diagnostic frames
    Author : Ayan
    Created: 28-08-24
    Example: run make_combined_animation_from_gdrive.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_combined_animation_from_gdrive.py --only_download
             run make_combined_animation_from_gdrive.py --field Par47
'''

from header import *
from util import *

from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from google_auth_oauthlib.flow import InstalledAppFlow, Flow
import io

from get_field_stats import natural_keys

start_time = datetime.now()


# --------------------------------------------------------------------------------------------------------------------
def get_credentials():
    '''
    Obtain your Google credentials
    Based on top answer in https://stackoverflow.com/questions/76485003/how-do-i-download-all-files-from-a-google-drive-folder-with-more-than-50-files
    '''
    SCOPES = ['https://www.googleapis.com/auth/drive.readonly'] # Define the scopes

    flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
    creds = flow.run_local_server(port=0)

    return creds

# --------------------------------------------------------------------------------------------------------------------
def query_google_drive_folder(folder_id, credentials):
    '''
    Queries google drive folder and returns file/folder list
    Based on top answer in https://stackoverflow.com/questions/76485003/how-do-i-download-all-files-from-a-google-drive-folder-with-more-than-50-files
    '''
    # Build the downloader
    drive_downloader = build('drive', 'v3', credentials=credentials)

    query = f"'{folder_id}' in parents"  # this works  ref https://stackoverflow.com/q/73119251/248616

    results = drive_downloader.files().list(q=query, pageSize=1000).execute()
    items = results.get('files', [])

    print(f'Found total {len(items)} files/folders in the folder.')

    return items, drive_downloader

# --------------------------------------------------------------------------------------------------------------------
def download_folder_from_google_drive(folder_id, destination_folder, credentials):
    '''
    Downloads the given file IDs from google drive
    Based on top answer in https://stackoverflow.com/questions/76485003/how-do-i-download-all-files-from-a-google-drive-folder-with-more-than-50-files
    '''
    destination_folder.mkdir(parents=True, exist_ok=True)
    items, drive_downloader = query_google_drive_folder(folder_id, credentials)

    for item in items:
        if item['mimeType'].endswith('.folder'):
            print(f'Downloading folder {item["name"]} from google drive..')
            download_folder_from_google_drive(item['id'], products_path / item['name'], credentials)
        else:
            request = drive_downloader.files().get_media(fileId=item['id'])
            f = io.FileIO(destination_folder / item['name'], 'wb')
            downloader = MediaIoBaseDownload(f, request)
            done = False
            while not done:
                status, done = downloader.next_chunk()
                print(f'Downloaded {item["name"]}: {int(status.progress() * 100)}%', end='\r')
            #print('\n')

    print(f'Downloaded {len(items)} files from the folder.')


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    fields_of_interest = [f'Par{item}' for item in passage_fields_in_cosmos]
    passage_url_id = '1oE1V76aeDlR8vxZNefyzk3yp0eO1rkHw'

    # ------determining directories and global variables---------
    args.plot_radial_profiles, args.only_seg, args.snr_cut = True, True, 3 #
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    snr_text = f'_snr{args.snr_cut:.1f}' if args.snr_cut is not None else ''
    description_text1 = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'
    description_text2 = f'diagnostics_and_extractions'

    # --------query passage folder and see which fields available------------------
    all_photcat_files = [args.input_dir / f'Par{int(field[3:]):03}' / 'Products' / f'Par{int(field[3:]):03}_photcat.fits' for field in fields_of_interest]
    if np.array([os.path.exists(file) for file in all_photcat_files]).all():
        field_list = fields_of_interest
        print(f'All fields are of interest have been downloaded already, so skipping the goggle drive query step.')
    else:
        credentials = get_credentials()
        items, _ = query_google_drive_folder(passage_url_id, credentials)
        field_url_dict = {item['name']:item['id'] for item in items if item['name'].startswith('Par') and item ['mimeType'].endswith('.folder')}
        fields_in_gdrive = list(field_url_dict.keys())
        print(f'..out of which {len(fields_in_gdrive)} are PASSAGE fields...')

        field_list = list(set(fields_of_interest).intersection(fields_in_gdrive))

    field_list = np.unique(field_list + [f'Par{int(item.split("Par")[1])}' for item in args.field_arr]).tolist() # adding some additional fields to the mix, if any
    field_list.sort(key=natural_keys)
    print(f'...out of which {len(field_list)} fields are of interest.')

    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}'
        print(f'\n\nCommencing field {args.field} which is {index+1} of {len(field_list)}..')

        # ----------determining filenames etc. to check for presence----------
        products_path = args.input_dir / args.field / 'Products'
        extraction_path = args.input_dir / args.field / 'Extractions'
        output_dir = args.output_dir / args.field
        file_to_check_for = f'{args.field}__{description_text2}_anim.mp4'

        if os.path.exists(output_dir / file_to_check_for):
            print(f'{file_to_check_for} already present for {args.field}, so skipping this field.')

        else:
            # ------------download the files------------------
            test_filename = products_path / f'{args.field}_photcat.fits'
            if os.path.exists(test_filename):
                print(f'Downloads already present, so proceeding to unzipping.')
            else:
                print(f'Downloading folder {field} from google drive..')
                folder_id = field_url_dict[field]
                download_folder_from_google_drive(folder_id, products_path, credentials=credentials)

            if args.only_download:
                print('--only_download option was used, hence skipping all other steps. Remove this option to access subsequent steps')
                continue

            # ------------unzip the downloaded files------------------
            zipped_files = glob.glob(str(products_path) + '/*.gz')
            if len(zipped_files) == 0:
                print(f'All files already unzipped, so proceeding to renaming.')
            else:
                for ind, thisfile in enumerate(zipped_files):
                    print(f'Unzipping {ind + 1} of {len(zipped_files)} zipped files..')
                    shutil.unpack_archive(thisfile, products_path)
                    os.remove(thisfile) #remove zipped files after unzipping

            # ------------rename the files within the downloaded folders------------------
            subfolders = glob.glob(str(products_path) + '/*/')
            for subfolder in subfolders:
                files_to_rename = glob.glob(str(subfolder) + '/' + field + '*')
                if len(files_to_rename) == 0:
                    print(f'All files within sub-folder {subfolder.split("/")[-2]} already have the correct nomenclature, so proceeding to next subfolder.')
                else:
                    print(f'Renaming files within sub-folder {subfolder.split("/")[-2]}/..')
                    for thisfile in files_to_rename:
                        os.rename(thisfile, thisfile.replace(field, args.field))

            # ------------rename the downloaded files------------------
            files_to_rename = glob.glob(str(products_path) + '/' + field + '*')
            if len(files_to_rename) == 0:
                print(f'All other files already have the correct nomenclature, so proceeding to making diagnostic plots.')
            else:
                print(f'Renaming other files..')
                for thisfile in files_to_rename:
                    os.rename(thisfile, thisfile.replace(field, args.field))

            # ------------run make_diagnostic_maps.py------------------
            diag_results_file = output_dir / f'{args.field}_all_diag_results.txt'
            if os.path.exists(diag_results_file) and not args.clobber:
                print(f'Diagnostic results file already present, so proceeding to making combined diagnostics and extraction images.')
            elif not (os.path.exists(products_path / 'maps') and os.path.exists(products_path / 'spec1D')):
                print(f'Supposed to run make_diagnostic_maps.py but either {str(products_path / "maps")} OR spec_1D/ does not exist, therefore cannot run. Moving to the next field.')
                continue
            else:
                print(f'Running make_diagnostic_maps.py..')
                dummy = subprocess.run(['python', 'make_diagnostic_maps.py', '--field', f'{args.field}', '--do_all_obj', '--plot_radial_profiles', '--only_seg', '--snr_cut', '3', '--write_file', '--hide'])

            # ------------run combine_diagnostics_and_extractions.py------------------
            diagnostic_img_files = glob.glob(str(output_dir / f'{description_text1}') + f'/{args.field}_*_{description_text1}.png')
            extraction_img_files = glob.glob(str(output_dir / f'{description_text2}') + f'/{args.field}_*_{description_text2}.png')
            if len(extraction_img_files) == len(diagnostic_img_files):
                print(f'All combined extraction images already present, so proceeding to the next step.')
            elif not (os.path.exists(products_path / 'plots') or os.path.exists(extraction_path)):
                print(f'Supposed to run combine_diagnostics_and_extractions.py but neither {extraction_path} nor {str(products_path / "plots")} exists, therefore cannot run. Moving to the next field.')
                continue
            else:
                print(f'Running combine_diagnostics_and_extractions.py..')
                dummy = subprocess.run(['python', 'combine_diagnostics_and_extractions.py', '--field', f'{args.field}', '--do_all_obj', '--hide'])

            # ------------make the final animation with the combined images------------------
            file_to_move = output_dir / f'{description_text2}' / f'{args.field}__{description_text2}_anim.mp4'
            if os.path.exists(file_to_move):
                print(f'Animation already present, but in {file_to_move}, so proceeding to moving it.')
            else:
                print(f'Running animate_png.py..')
                dummy = subprocess.run(['python', '/Users/acharyya/Work/astro/ayan_codes/animate_png.py', '--inpath', f'{output_dir}/{description_text2}/', '--rootname', f'{args.field}_*_{description_text2}.png', '--delay', '0.1'])

            # ------------move the final animation------------------
            print(f'Moving the animation file to {file_to_check_for}..')
            dummy = shutil.move(file_to_move, output_dir / file_to_check_for)

        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'All {len(field_list)} fields done. Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
