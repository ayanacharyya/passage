import os, sys, tarfile, fnmatch, shutil
from pathlib import Path
from box_sdk_gen import BoxClient, BoxDeveloperTokenAuth
from datetime import datetime, timedelta

start_time = datetime.now()

# --- CONFIGURATION ---
BOX_FOLDER_ID = '305422038072'  # Extracted directly from your URL
BOX_SHARED_LINK = 'https://app.box.com/folder/305422038072'
DEVELOPER_TOKEN = '0UBpObUj4UlonxBFJBqeOPo6FD3Q3rHd' # Expires: 8 June 2026 at 04:03:15 GMT-7; Generate a new one here https://app.box.com/developers/console/app/2590971/configuration

FIELDS = ['Par003', 'Par017', 'Par023', 'Par024', 'Par025', 'Par026',
          'Par028', 'Par029', 'Par051', 'Par052', 'Par053']

BASE_TARGET_DIR = Path("/Volumes/Ayan_SSD/Ayan_macbook/Users/aacharyya/Work/astro/passage/passage_data/v0.6")
# ---------------------

def get_final_clean_name(raw_name):
    """Helper to predict what the file name will look like after processing."""
    if '_fixz' in raw_name:
        return raw_name.replace('_fixz', '')
    return raw_name

def clean_unpacked_directories(target_dir):
    """
    Scans the directory for any subfolders containing '_fixz' 
    and renames them bottom-up so path changes don't break the loop.
    """
    # top_down=False ensures we rename inner subdirectories before outer parent directories
    for root, dirs, files in os.walk(target_dir, topdown=False):
        for dir_name in dirs:
            if '_fixz' in dir_name:
                old_dir_path = Path(root) / dir_name
                clean_dir_name = dir_name.replace('_fixz', '')
                new_dir_path = Path(root) / clean_dir_name
                
                # Check if the clean directory name already exists before moving
                if not new_dir_path.exists():
                    old_dir_path.rename(new_dir_path)
                    print(f"    Renamed Folder: {dir_name} -> {clean_dir_name}")
                else:
                    print(f"    Folder {clean_dir_name} already exists; merged content manually if needed.")

# ---------------------
if __name__ == '__main__':
    # Authenticate via the modern token manager
    auth = BoxDeveloperTokenAuth(token=DEVELOPER_TOKEN)
    client = BoxClient(auth=auth)
    
    # Pass the shared link header so Box authenticates our access to this ID
    shared_headers = {"Box-Shared-Link": BOX_SHARED_LINK}
    
    print(f"Connecting directly to Box Folder ID: {BOX_FOLDER_ID}...")
    download_count = 0
    skipped_count = 0

    try:
        # 1. Directly gather all top-level elements inside the folder ID
        root_items_page = client.folders.get_folder_items(folder_id=BOX_FOLDER_ID, extra_headers=shared_headers)
    except Exception as e:
        print(f"Connection failed. Is your Developer Token expired?\nDetails: {e}")
        sys.exit()

    print("Scanning fields...")
    
    for item in root_items_page.entries:
        if item.type == 'folder' and item.name in FIELDS:
            field_name = item.name
            print(f"\n📁 Found target field directory: {field_name}")
            
            # Form your destination local directories
            products_dir = BASE_TARGET_DIR / field_name / "Products"
            products_dir.mkdir(parents=True, exist_ok=True)
            
            # 2. Open up the specific field sub-folder on Box
            sub_items_page = client.folders.get_folder_items(folder_id=item.id, extra_headers=shared_headers)
            
            for sub_item in sub_items_page.entries:

                if sub_item.type != 'file':
                    continue
                
                # Check match criteria: _fixz.tar.gz OR _photcat.fits OR *drz*.fits
                name_lower = sub_item.name.lower()
                is_tarball = '.tar.gz' in name_lower or name_lower.endswith('.tgz')
                is_photcat = name_lower.endswith('_photcat.fits')
                is_drz = fnmatch.fnmatch(name_lower, '*drz*.fits')
                
                if not (('fixz' in name_lower) or is_photcat or is_drz):
                    continue # Skip unrelated files in the directory

                print(f"  Found matching file: {sub_item.name}")
                temp_file_path = products_dir / sub_item.name

                # Predict final destination path to determine resume/skip
                final_local_name = get_final_clean_name(sub_item.name)

                if is_tarball:
                    # Tarballs extract contents. Let's guess the internal name by stripping '.tar.gz'
                    # e.g., Par003_maps_fixz.tar.gz -> Par003_maps_fixz -> Par003_maps
                    base_archive_name = sub_item.name.split('.tar.gz')[0].split('.tgz')[0]
                    predicted_unpacked_name = get_final_clean_name(base_archive_name)
                    
                    # We check if a directory or file starting with that clean string exists
                    already_exists = any(products_dir.glob(f"{predicted_unpacked_name}*"))
                else:
                    # Flat FITS images check
                    already_exists = (products_dir / final_local_name).exists()
                
                if already_exists:
                    print(f'{final_local_name} already exists, so skipping this download..')
                    skipped_count += 1
                    continue # Seamless resume skip!

                print(f"    Downloading to local drive...")   
                # 4. Stream down the file directly using the modern downloads engine
                file_stream = client.downloads.download_file(file_id=sub_item.id, extra_headers=shared_headers)

                with open(temp_file_path, 'wb') as f:
                    # shutil automatically optimizes the network buffer size (typically 16KB-64KB chunks)
                    # and streams it directly to disk without getting stuck in an idle loop
                    shutil.copyfileobj(file_stream, f)
                download_count += 1

                # Handle Archive extraction
                if is_tarball:
                    print("    Unpacking tar.gz archive contents...")
                    try:
                        with tarfile.open(temp_file_path, 'r:gz') as tar_ref:
                            tar_ref.extractall(path=products_dir)
                            extracted_files = tar_ref.getnames()
                        
                        temp_file_path.unlink() # Delete original tarball wrapper
                        
                        for ext_file in extracted_files:
                            pure_filename = Path(ext_file).name
                            local_ext_path = products_dir / pure_filename
                            
                            if local_ext_path.is_file() and '_fixz' in local_ext_path.name:
                                clean_name = local_ext_path.name.replace('_fixz', '')
                                local_ext_path.rename(local_ext_path.with_name(clean_name))
                                print(f"    Unpacked & Cleaned: {clean_name}")
                        # Clean up any nested directories extracted from the tarball
                        clean_unpacked_directories(products_dir)
                    except Exception as extract_err:
                        print(f"    Error unpacking file: {extract_err}")
                
                # Handle straight flat files (FITS maps/catalogs)
                else:
                    if '_fixz' in temp_file_path.name:
                        clean_name = temp_file_path.name.replace('_fixz', '')
                        temp_file_path.rename(temp_file_path.with_name(clean_name))
                        print(f"    Cleaned filename: {clean_name}")
                    else:
                        print(f"    Saved: {temp_file_path.name}")

    print(f"\n Processing complete in {timedelta(seconds=(datetime.now() - start_time).seconds)}!")
    print(f" Total files skipped (already completed): {skipped_count}")
    print(f" Total new files successfully fetched: {download_count}")
