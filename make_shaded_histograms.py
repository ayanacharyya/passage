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
import hdbscan
from scipy.stats import binned_statistic

HOME = Path.home()
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def clusterer(df):
    '''
    To computing clustering of the substructures
    Returns 2 numbers i.e., number of clusters and substructures
    Clustering algorithm from Beena
    '''

    mcs, ms = 5, 5
    yc1 = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, cluster_selection_method='eom', cluster_selection_epsilon=30).fit(df[['XPRIME', 'YPRIME']])
    yc2 = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, cluster_selection_method='leaf').fit(df[['XPRIME', 'YPRIME']])
    clusters_AC1 = len(np.unique(yc1.labels_)) - 1 # number of EOM clusters in this age group
    clusters_AC2 = len(np.unique(yc2.labels_)) - 1 # number of leaf clusters in this age group

    return clusters_AC1, clusters_AC2

# --------------------------------------------------------------------------------------------------------------------
def make_histogram(df, bin_edges, silent=False):
    '''
    To computing histogram of the number of substructures for every age bin, for a given list of bin edges
    Returns number of clusters in each bin as array, for given bin_edges
    '''

    nclusters_arr1, nclusters_arr2, nstars_arr = [], [], []

    # --------splitting the dataframe for each age bin--------
    for index in range(len(bin_edges) - 1):
        left_edge = bin_edges[index]
        right_edge = bin_edges[index + 1]
        df_split = df[df['age'].between(left_edge, right_edge)]

        nc1_thisbin, nc2_thisbin = clusterer(df_split)
        nstars_thisbin = len(df_split)

        nclusters_arr1.append(nc1_thisbin)
        nclusters_arr2.append(nc2_thisbin)
        nstars_arr.append(nstars_thisbin)

        if not silent: print(f'For bin [{left_edge}, {right_edge}]: {nstars_thisbin} total stars, {nc1_thisbin} clusters & {nc2_thisbin} in substructures.')

    return nclusters_arr1, nclusters_arr2, nstars_arr

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
def collapse_and_plot_histogram(data_and_bins, bin_edges, ax, color='grey', alpha=0.3, lw=2):
    '''
    Collapse a 2D array oh histograms and plot it as a dark line + shaded region, on a given axis
    Returns axis handle
    '''
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    data_and_bins = np.array(data_and_bins)

    #bins = np.digitize(data[:, 1, :], bin_edges)
    data = data_and_bins[:, 0, :].flatten()
    bins = data_and_bins[:, 1, :].flatten()
    median_arr = binned_statistic(bins, data, statistic='median', bins=bin_edges).statistic
    std_arr = binned_statistic(bins, data, statistic='std', bins=bin_edges).statistic

    ax.plot(bin_centers, median_arr, color=color, alpha=1, lw=lw)
    ax.fill_between(bin_centers, median_arr - std_arr, median_arr + std_arr, color=color, alpha=alpha)

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    # ------------initialising user input variables-------------------
    col_arr = ['salmon', 'cornflowerblue', 'khaki']
    alpha, lw, fontsize = 0.1, 2, 15
    min_shift, max_shift = 0, 10 # amounts to shift each bin edge by, in Myr
    input_file_path = HOME / 'Downloads'
    output_path = None
    age_max = 100 # Myr

    # ------------initialising user input flags-------------------
    plot_bin_edges = True
    plot_individual_iterations = True
    silent = True
    keep = False
    normalise = True

    # ------------reading data------------------
    filename = 'ngc4449_deproject_ages.txt'
    df = read_catalog(input_file_path / filename, age_max=age_max) # age_max in Myr

    # ------------declaring varieties of base bin edges------------------
    fiducial_bin_edges = np.array([0, 5, 10, 15, 20, 30, 40, 50, 75, 100])

    bigger_bin_edges = np.hstack([fiducial_bin_edges[::2], [fiducial_bin_edges[-1]]])

    intermediate_bins = (fiducial_bin_edges[:-1] + fiducial_bin_edges[1:])/2
    smaller_bin_edges = np.sort(np.hstack([fiducial_bin_edges, intermediate_bins]))

    all_base_bins = [bigger_bin_edges, fiducial_bin_edges, smaller_bin_edges]

    # --------initialising figure and variables-------------------
    if not keep: plt.close('all')
    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    shifts = np.arange(min_shift, max_shift + 1)

    # ------------loop over base bins------------------
    for index, base_bin_edges in enumerate(all_base_bins):
        print(f'\nDoing base-bin {base_bin_edges}, which is {index+1} out of {len(all_base_bins)}')

        nc1_and_bins_arr, nc2_and_bins_arr, nstars_and_bins_arr = [], [], []
        # -----------loop over shifted bins------------------------
        for index2, this_shift in enumerate(shifts):
            #this_bin_edges = np.hstack([base_bin_edges[0], base_bin_edges[1:-1] + this_shift, base_bin_edges[-1]]) # shift all but first and last bin edge, i.e. first bin gets larger and last bin gets smaller
            this_bin_edges = base_bin_edges + this_shift # shift ALL bins, i.e. all bin sizes stay the same, but some data is lost from the first bin as it moves up
            this_bin_centers = (this_bin_edges[:-1] + this_bin_edges[1:])/2
            print(f'Shifting base-bin-{index+1} by {this_shift}, which is {index2 + 1} out of {len(shifts)} shifts')

            nc1_arr_thisbin, nc2_arr_thisbin, nstars_arr_thisbin = make_histogram(df, this_bin_edges, silent=silent)

            # ------normalising by the max--------
            if normalise:
                nc1_arr_thisbin /= np.max(nc1_arr_thisbin)
                nc2_arr_thisbin /= np.max(nc2_arr_thisbin)
                nstars_arr_thisbin /= np.max(nstars_arr_thisbin)

            # ------appending each shifted iteration--------
            nc1_and_bins_arr.append([nc1_arr_thisbin, this_bin_centers])
            nc2_and_bins_arr.append([nc2_arr_thisbin, this_bin_centers])
            nstars_and_bins_arr.append([nstars_arr_thisbin, this_bin_centers])

            # -----plotting histograms of individual iterations--------
            if plot_individual_iterations:
                axes[index].plot(this_bin_centers, nc1_arr_thisbin, c=col_arr[0], alpha=alpha)
                axes[index].plot(this_bin_centers, nc2_arr_thisbin, c=col_arr[1], alpha=alpha)
                axes[index].plot(this_bin_centers, nstars_arr_thisbin, c=col_arr[2], alpha=alpha)

            # -----decorating axes tick labels etc--------
            axes[index].tick_params(axis='both', which='major', labelsize=fontsize)

        # -------plotting bin edges-----------------
        if plot_bin_edges:
            for edge in base_bin_edges: axes[index].axvline(edge, ls='dashed', lw=0.5, c='grey')

        # -----plotting histograms from variants of this base bin--------
        axes[index] = collapse_and_plot_histogram(nc1_and_bins_arr, base_bin_edges, axes[index], color=col_arr[0], alpha=alpha, lw=lw)
        axes[index] = collapse_and_plot_histogram(nc2_and_bins_arr, base_bin_edges, axes[index], color=col_arr[1], alpha=alpha, lw=lw)
        axes[index] = collapse_and_plot_histogram(nstars_and_bins_arr, base_bin_edges, axes[index], color=col_arr[2], alpha=alpha, lw=lw)

        # -------annotate figure and save---------------
        axes[0].text(0.99, 0.98, 'Clusters', color=col_arr[0], ha='right', va='top', fontsize=fontsize, transform=axes[0].transAxes)
        axes[0].text(0.99, 0.85, 'Substructures', color=col_arr[1], ha='right', va='top', fontsize=fontsize, transform=axes[0].transAxes)
        axes[0].text(0.99, 0.70, 'Total stars', color=col_arr[2], ha='right', va='top', fontsize=fontsize, transform=axes[0].transAxes)

        if index == len(all_base_bins) - 1: axes[index].set_xlabel('Age (Myr)', fontsize=fontsize)
        axes[index].set_ylabel('Normalised #', fontsize=fontsize)

        axes[index].set_xlim(0, age_max * 1.05)
        if normalise: axes[index].set_ylim(0, 1.1)

    fig.subplots_adjust(left=0.11, bottom=0.1, top=0.98, right=0.98, hspace=0.05, wspace=0.0)

    if output_path is None: output_path = input_file_path
    figname = output_path / f'{filename.split("_")[0]}_clustering_histograms.png'
    fig.savefig(figname)
    print(f'Saved figure at {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
