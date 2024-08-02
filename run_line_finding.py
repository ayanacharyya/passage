'''
    Filename: run_line_finding.py.py
    Notes: Runs Kalina's interactive line-finding algorithm (Based on Kalina's mainPASSAGE.py file)
    Author : Ayan
    Created: 02-08-24
    Example: run run_line_finding.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
'''

from header import *
from util import *

# --------to make sure wisp_analysis and mpfit can be imported------------
LF_DIR = str(HOME / 'Work/astro/passage/line-finding/')
WISP_DIR = LF_DIR + '/wisp_analysis'
sys.path.append(LF_DIR)
sys.path.append(WISP_DIR)
import wisp_analysis as wisp

# -------To open two ds9 windows---------
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')

# To run the code, first run wisp.loop_field_cwt() to create the line list for thr field you're working on
# Once that is done, that line can be commented out, unless the line list needs to be recreated (e.g. if you're changing fields).
# Next, run wisp.measure_z_interactive() to fit the spectra.

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    parno = args.field[3:]
    #output_dir = args.output_dir / args.field

    os.chdir(args.output_dir) # move to the directory where you want your outputs
    wisp.loop_field_cwt(path_to_wisp_data=str(args.input_dir) + '/', path_to_code=WISP_DIR, parno=parno)
    wisp.measure_z_interactive(path_to_wisp_data=str(args.input_dir) + '/', path_to_code=WISP_DIR, parno=parno)

