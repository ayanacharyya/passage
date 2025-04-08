'''
    Filename: compare_metallicities.py
    Notes: Compares spatially resolved metallicity map across different diagnostics for a given set of galaxies
    Author : Ayan
    Created: 02-04-25
    Example: run plots_for_zgrad_paper.py --field Par028 --id 300,1303,1634,1849,2171,2727,2867 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,R3,P25,NB --use_original_NB_grid --colorcol radius --debug_Zdiag
             run plots_for_zgrad_paper.py --field Par028 --id 300 --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --AGN_diag Ne3O2 --Zdiag O3O2,R3,NB --Zbranch low --colorcol radius --debug_Zdiag --exclude_line SII
'''
from make_diagnostic_maps import plot_2D_map, plot_radial_profile
from header import *
from util import *

plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
