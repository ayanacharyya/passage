'''
    Filename: download_ngdeep_pngs.py
    Notes: Download all .png files from a password-protected NGDEEP/Grizli web directory
           and organize them into subfolders by product type (full/, sed/, stack/, 1d/, etc.)
    Author : Gemini, modified by Ayan    
    Created: 23-11-26
    Example: run download_ngdeep_pngs.py --recursive --types full sed
             run download_ngdeep_pngs.py --url https://web.corral.tacc.utexas.edu/ngdeepinternal/NIRISSData/Grizli_Extractions/epoch2/v1.0/ --output_dir /Users/acharyya/Work/astro/passage/passage_data/v1.0/NGDEEP/plots_v1.0
             run download_ngdeep_pngs.py
'''

import os
import re
import getpass
import argparse
from pathlib import Path
from datetime import datetime, timedelta
from urllib.parse import urljoin, unquote, urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from requests.auth import HTTPBasicAuth, HTTPDigestAuth
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

try:
    from bs4 import BeautifulSoup
    HAS_BS4 = True
except ImportError:
    HAS_BS4 = False

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description='Download and categorize Grizli PNG files from NGDEEP portal.')
    parser.add_argument('--url', type=str, default='https://web.corral.tacc.utexas.edu/ngdeepinternal/NIRISSData/Grizli_Extractions/epoch2/v1.2/Extractions_wrefit/', help='Base URL of the remote directory')
    parser.add_argument('--out_dir', type=Path, default='/Users/acharyya/Work/astro/passage/passage_data/v1.0/NGDEEP/plots_v1.2', help='Target local directory to store subfolders')
    parser.add_argument('--username', type=str, default='ngdeepHUDF', help='Portal username (or set NGDEEP_USER)')
    parser.add_argument('--password', type=str, default='ngdeep4theW!n', help='Portal password (or set NGDEEP_PASS)')
    parser.add_argument('--types', nargs='+', default=None, help='Optional subset of PNG types to download (e.g., full sed stack)')
    parser.add_argument('--recursive', action='store_true', help='Also crawl subdirectories under the base URL')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel download threads (default: 8)')
    parser.add_argument('--overwrite', action='store_true', help='Re-download files even if they already exist locally')
    return parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
def build_session(username, password, test_url=None):
    '''
    Create a requests.Session with automatic Basic/Digest auth detection and retries
    '''
    session = requests.Session()
    session.headers.update({'User-Agent': 'Mozilla/5.0'})
    retry = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=['GET', 'HEAD'],
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
    session.mount('https://', adapter)
    session.mount('http://', adapter)

    session.auth = HTTPBasicAuth(username, password)
    if test_url is not None:
        probe = session.get(test_url, timeout=30)
        www_auth = probe.headers.get('WWW-Authenticate', '')
        if probe.status_code == 401 and 'digest' in www_auth.lower():
            session.auth = HTTPDigestAuth(username, password)
        elif probe.status_code == 401:
            print(f'Auth failed (WWW-Authenticate: {www_auth}). Check username/password quoting.')

    return session

# --------------------------------------------------------------------------------------------------------------------
def extract_links(html_text, base_url):
    '''
    Parse all href targets from an Apache/Nginx directory index page
    '''
    if HAS_BS4:
        soup = BeautifulSoup(html_text, 'html.parser')
        raw_hrefs = [a.get('href') for a in soup.find_all('a', href=True)]
    else:
        raw_hrefs = re.findall(r'href=["\']([^"\']+)["\']', html_text, flags=re.IGNORECASE)

    full_urls = []
    for href in raw_hrefs:
        # Skip Apache column-sorting query links, anchors, and parent directory links
        if not href or href.startswith(('?', '#', '/')) or href in ('../', './'):
            continue
        resolved = urljoin(base_url, href)
        # Ensure we stay strictly inside base_url
        if resolved.startswith(base_url) and resolved != base_url:
            full_urls.append(resolved)

    return list(dict.fromkeys(full_urls))

# --------------------------------------------------------------------------------------------------------------------
def collect_png_urls(session, base_url, recursive=False):
    '''
    Collect all .png file URLs from base_url (and subdirectories if recursive=True)
    '''
    if not base_url.endswith('/'):
        base_url += '/'

    png_urls = []
    visited_dirs = set()
    dirs_to_visit = [base_url]

    while dirs_to_visit:
        current_dir = dirs_to_visit.pop(0)
        if current_dir in visited_dirs:
            continue
        visited_dirs.add(current_dir)

        print(f'Scanning {current_dir} ...')
        resp = session.get(current_dir, timeout=30)
        resp.raise_for_status()

        for link in extract_links(resp.text, current_dir):
            path_part = urlparse(link).path
            if path_part.lower().endswith('.png'):
                png_urls.append(link)
            elif recursive and path_part.endswith('/'):
                dirs_to_visit.append(link)

    return list(dict.fromkeys(png_urls))

# --------------------------------------------------------------------------------------------------------------------
def get_subfolder_name(filename):
    '''
    Infer subfolder name from Grizli PNG naming convention:
      e.g., ngdeep-v1.0_00123.full.png  -> 'full'
            ngdeep-v1.0_00123.sed.png   -> 'sed'
            ngdeep-v1.0_00123.1d.png    -> '1d'
            ngdeep-v1.0_00123_sed.png   -> 'sed'
            ngdeep-v1.0_00123.png       -> 'other'
    '''
    stem = Path(filename).stem  # removes trailing .png
    # 1. Standard dot-separated Grizli product tag (e.g., .full, .sed, .stack, .1d, .2d, .line, .beams)
    match = re.search(r'\.((?:1d|2d|[a-zA-Z][a-zA-Z0-9-]{0,15}))$', stem)
    if match:
        return match.group(1).lower()

    # 2. Fallback for underscore-separated product tag (e.g., _sed, _full)
    match = re.search(r'_((?:1d|2d|[a-zA-Z][a-zA-Z0-9-]{0,15}))$', stem)
    if match:
        return match.group(1).lower()

    return 'other'

# --------------------------------------------------------------------------------------------------------------------
def download_single_file(session, file_url, target_path, overwrite=False):
    '''
    Stream a single PNG file to target_path; skips if already present unless overwrite=True
    '''
    if target_path.exists() and target_path.stat().st_size > 0 and not overwrite:
        return 'skipped', target_path

    target_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target_path.with_suffix(target_path.suffix + '.part')

    with session.get(file_url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(tmp_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=65536):
                if chunk:
                    f.write(chunk)

    tmp_path.replace(target_path)
    return 'downloaded', target_path

# --------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    start_time = datetime.now()
    args = parse_args()
    args.out_dir.mkdir(exist_ok=True, parents=True)

    username = args.username or input('Portal username: ').strip()
    password = args.password or getpass.getpass('Portal password: ')

    session = build_session(username, password, test_url=args.url)

    # ----------scrape directory listing for PNG files----------
    png_urls = collect_png_urls(session, args.url, recursive=args.recursive)
    print(f'\nFound {len(png_urls)} total .png files.')

    if len(png_urls) == 0:
        print('No .png files found at top level. Try re-running with --recursive if files are inside subfolders.')
        raise SystemExit(0)

    # ----------map each URL to its destination subfolder----------
    allowed_types = {t.lower().lstrip('.') for t in args.types} if args.types else None
    tasks = []
    category_counts = {}

    for url in png_urls:
        filename = unquote(Path(urlparse(url).path).name)
        subfolder = get_subfolder_name(filename)
        if allowed_types is not None and subfolder not in allowed_types:
            continue
        dest_path = args.out_dir / subfolder / filename
        tasks.append((url, dest_path))
        category_counts[subfolder] = category_counts.get(subfolder, 0) + 1

    print('Planned file breakdown by subfolder:')
    for sub, count in sorted(category_counts.items()):
        print(f'  {sub}/ : {count} files')

    # ----------download files in parallel----------
    print(f'\nStarting download of {len(tasks)} files into {args.out_dir} using {args.workers} workers...')
    downloaded_count, skipped_count, failed_count = 0, 0, 0

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_single_file, session, url, dest, args.overwrite): (url, dest)
            for url, dest in tasks
        }
        for idx, future in enumerate(as_completed(futures), start=1):
            url, dest = futures[future]
            try:
                status, _ = future.result()
                if status == 'downloaded':
                    downloaded_count += 1
                else:
                    skipped_count += 1
            except Exception as exc:
                failed_count += 1
                print(f'  [FAILED] {dest.name}: {exc}')

            if idx % 200 == 0 or idx == len(tasks):
                print(f'  Progress: {idx}/{len(tasks)} (downloaded={downloaded_count}, skipped={skipped_count}, failed={failed_count})')

    print(f'\nDone! Downloaded={downloaded_count}, Skipped={skipped_count}, Failed={failed_count} in {timedelta(seconds=(datetime.now() - start_time).seconds)}')