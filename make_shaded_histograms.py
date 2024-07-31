'''
    Filename: make_shaded_histograms.py
    Notes: Convert given histograms into shaded KDE plots (for Beena)
    Author : Ayan
    Created: 01-08-24
    Example: run make_shaded_histograms.py
'''

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path

HOME = Path.home()
start_time = datetime.now()

try:
    import hdbscan
except:
    print('BEWARE: generating dummy data because HDBSCAN did not work.')
    pass

# --------------------------------------------------------------------------------------------------------------------
def clusterer(df):
    '''
    To computing clustering of the substructures
    Returns 2 numbers i.e., number of clusters and substructures
    '''

    try:
        # ----clustering algo from Beena------
        mcs, ms = 5, 5
        yc1 = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, cluster_selection_method='eom', cluster_selection_epsilon=30).fit(df[['XPRIME', 'YPRIME']])
        yc2 = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, cluster_selection_method='leaf').fit(df[['XPRIME', 'YPRIME']])
        clusters_AC1 = len(np.unique(yc1.labels_)) - 1 # number of EOM clusters in this age group
        clusters_AC2 = len(np.unique(yc2.labels_)) - 1 # number of leaf clusters in this age group
    except:
        clusters_AC1, clusters_AC2 = [int(item * 10) for item in np.random.random(2)] # random data generator

    return clusters_AC1, clusters_AC2

# --------------------------------------------------------------------------------------------------------------------
def make_histogram(df, bin_edges):
    '''
    To computing histogram of the number of substructures for every age bin, for a given list of bin edges
    Returns number of clusters in each bin as array, for given bin_edges
    '''

    nclusters_arr1, nclusters_arr2 = [], []

    # --------splitting the dataframe for each age bin--------
    for index in range(len(bin_edges) - 1):
        left_edge = bin_edges[index]
        right_edge = bin_edges[index + 1]
        df_split = df[df['age'].between(left_edge, right_edge)]

        nc1_thisbin, nc2_thisbin = clusterer(df_split)

        nclusters_arr1.append(nc1_thisbin)
        nclusters_arr2.append(nc2_thisbin)

    return nclusters_arr1, nclusters_arr2

# ----------------------------------------------------------------------------------
def read_catalog(filename, age_max=None):
    '''
    Reads in the star catalog
    Returns pandas dataframe with only the relevant columns
    '''

    df = pd.read_table(filename, comment='#')
    df  = df[['XPRIME', 'YPRIME', 'age']]
    df = df.sort_values(by='age') # age in Myr

    if age_max is not None: df = df[df['age'] <= age_max] # curtails dataframe to a given max age, to reduce memory requirements

    return df

# ----------------------------------------------------------------------------------
def collapse_and_plot_histogram(data, ax, color='grey', alpha=0.3, lw=2):
    '''
    Collapse a 2D array oh histograms and plot it as a dark line + shaded region, on a given axis
    Returns axis handle
    '''
    median_arr = np.median(data, axis=0)
    std_arr = np.std(data, axis=0)

    ax.plot(this_bin_centers, median_arr, color=color, alpha=1, lw=lw)
    ax.fill_between(this_bin_centers, median_arr - std_arr, median_arr + std_arr, color=color, alpha=alpha)

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    # ------------initialising user input variables-------------------
    col_arr = ['salmon', 'cornflowerblue', 'khaki']
    alpha, lw, fontsize = 0.1, 2, 10
    plot_individual_iterations, silent, keep = False, True, False
    min_shift, max_shift = 1, 10 # amounts to shift each bin edge by, in Myr

    # ------------reading data------------------
    filepath = HOME / 'Downloads/ngc4449_deproject_ages.txt'
    df = read_catalog(filepath, age_max=100) # age_max in Myr

    # ------------declaring varieties of base bin edges------------------
    fiducial_bin_edges = np.array([0, 5, 10, 15, 20, 30, 40, 50, 75, 100])

    bigger_bin_edges = np.hstack([fiducial_bin_edges[::2], [fiducial_bin_edges[-1]]])

    intermediate_bins = (fiducial_bin_edges[:-1] + fiducial_bin_edges[1:])/2
    smaller_bin_edges = np.sort(np.hstack([fiducial_bin_edges, intermediate_bins]))

    all_base_bins = [fiducial_bin_edges, bigger_bin_edges, smaller_bin_edges]

    # --------initialising figure and variables-------------------
    if not keep: plt.close('all')
    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    shifts = np.arange(min_shift, max_shift + 1)

    # ------------loop over base bins------------------
    for index, base_bins in enumerate(all_base_bins):
        print(f'\nDoing base bin {index+1} out of {len(all_base_bins)}')
        if not silent: print('Base bins are', base_bins)

        nc1_arr, nc2_arr = [], []
        # -----------loop over shifted bins------------------------
        for index2, this_shift in enumerate(shifts):
            this_bin_edges = np.hstack([base_bins[0], base_bins[1:-1] + this_shift, base_bins[-1]])
            this_bin_centers = (this_bin_edges[:-1] + this_bin_edges[1:])/2
            print(f'Shifting base bin {index+1} by {this_shift}, which is {index2 + 1} out of {len(shifts)} shifts')
            if not silent: print('Bins are', this_bin_edges)

            nc1_arr_thisbin, nc2_arr_thisbin = make_histogram(df, this_bin_edges)
            nc1_arr.append(nc1_arr_thisbin)
            nc2_arr.append(nc2_arr_thisbin)

            if plot_individual_iterations:
                axes[index].plot(this_bin_centers, nc1_arr_thisbin, c=col_arr[0], alpha=alpha)
                axes[index].plot(this_bin_centers, nc2_arr_thisbin, c=col_arr[1], alpha=alpha)

        # -----plotting histograms from variants of this base bin--------
        axes[index] = collapse_and_plot_histogram(nc1_arr, axes[index], color=col_arr[0], alpha=alpha, lw=lw)
        axes[index] = collapse_and_plot_histogram(nc2_arr, axes[index], color=col_arr[1], alpha=alpha, lw=lw)

        # -------figure labels---------------
        if index == len(all_base_bins) - 1: axes[index].set_xlabel('Age (Myr)', fontsize=fontsize)
        axes[index].set_ylabel('# counts', fontsize=fontsize)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
