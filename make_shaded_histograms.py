'''
    Filename: make_shaded_histograms.py
    Notes: Compute histograms by shifting bins from a given star catalog, and plot shaded histograms (for Beena)
    Author : Ayan
    Created: 01-08-24
    How to use:
    Detail example: run make_shaded_histograms.py --normalise --plot_bin_edges --plot_individual_iterations --keep --verbose --input_dir /Users/acharyya/Downloads/ --min_shift 0 --max_shift 10 --fontsize 15 --alpha 0.1 --age_max 100 --color_arr salmon,cornflowerblue,khaki
    Short example: run make_shaded_histograms.py --normalise --plot_bin_edges --plot_individual_iterations

'''

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path
import hdbscan
from scipy.stats import binned_statistic
import argparse

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
def make_histogram(df, bin_edges, verbose=False):
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

        if verbose: print(f'For bin [{left_edge}, {right_edge}]: {nstars_thisbin} total stars, {nc1_thisbin} clusters & {nc2_thisbin} in substructures.')

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
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces clustering histograms from stellar catalog.')

    parser.add_argument('--input_dir', metavar='input_dir', type=str, action='store', default=str(HOME / 'Downloads'), help='Where does the input file reside?')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, action='store', default=None, help='Where do you want to store the outputs?')
    parser.add_argument('--input_file', metavar='input_file', type=str, action='store', default='ngc4449_deproject_ages.txt', help='What is the filename of the input catalog?')

    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=15, help='fontsize of plot labels, etc.; default is 15')
    parser.add_argument('--alpha', metavar='alpha', type=float, action='store', default=0.1, help='alpha for the individual histograms, if plotted; default is 0.1')
    parser.add_argument('--color_arr', metavar='color_arr', type=str, action='store', default='salmon,cornflowerblue,khaki', help='Colors to be used for cluster, substructures and total stars--in order, separated by comma; Default is salom, cornflowerblue and khaki')

    parser.add_argument('--min_shift', metavar='min_shift', type=int, action='store', default=0, help='Minimum amount of shift (in Myr) of the histogram bin edges; default is 0')
    parser.add_argument('--max_shift', metavar='max_shift', type=int, action='store', default=10, help='Maximum amount of shift (in Myr) of the histogram bin edges; default is 10')
    parser.add_argument('--age_max', metavar='age_max', type=float, action='store', default=100, help='Maximum age (in Myr) which the analysis is restricted to; default is 100')

    parser.add_argument('--plot_bin_edges', dest='plot_bin_edges', action='store_true', default=False, help='Plot bin edges as vertical lines? Default is no.')
    parser.add_argument('--plot_individual_iterations', dest='plot_individual_iterations', action='store_true', default=False, help='Plot histograms for individual shifted bin iterations? Default is no.')
    parser.add_argument('--normalise', dest='normalise', action='store_true', default=False, help='Normalise the histograms so that peak =1? Default is no.')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=False, help='Suppress some generic print statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    args.input_dir = Path(args.input_dir)
    if args.output_dir is None: args.output_dir = args.input_dir
    else: args.output_dir = Path(args.output_dir)
    args.color_arr = [item for item in args.color_arr.split(',')]

    return args

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    # ------------initialising user input variables-------------------
    lw = 2
    args = parse_args()

    # ------------reading data------------------
    df = read_catalog(args.input_dir / args.input_file, age_max=args.age_max) # age_max in Myr

    # ------------declaring varieties of base bin edges------------------
    fiducial_bin_edges = np.array([0, 5, 10, 15, 20, 30, 40, 50, 75, 100])

    bigger_bin_edges = np.hstack([fiducial_bin_edges[::2], [fiducial_bin_edges[-1]]])

    intermediate_bins = (fiducial_bin_edges[:-1] + fiducial_bin_edges[1:])/2
    smaller_bin_edges = np.sort(np.hstack([fiducial_bin_edges, intermediate_bins]))

    all_base_bins = [bigger_bin_edges, fiducial_bin_edges, smaller_bin_edges]

    # --------initialising figure and variables-------------------
    if not args.keep: plt.close('all')
    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    shifts = np.arange(args.min_shift, args.max_shift + 1)

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

            nc1_arr_thisbin, nc2_arr_thisbin, nstars_arr_thisbin = make_histogram(df, this_bin_edges, verbose=args.verbose)

            # ------normalising by the max--------
            if args.normalise:
                nc1_arr_thisbin /= np.max(nc1_arr_thisbin)
                nc2_arr_thisbin /= np.max(nc2_arr_thisbin)
                nstars_arr_thisbin /= np.max(nstars_arr_thisbin)

            # ------appending each shifted iteration--------
            nc1_and_bins_arr.append([nc1_arr_thisbin, this_bin_centers])
            nc2_and_bins_arr.append([nc2_arr_thisbin, this_bin_centers])
            nstars_and_bins_arr.append([nstars_arr_thisbin, this_bin_centers])

            # -----plotting histograms of individual iterations--------
            if args.plot_individual_iterations:
                axes[index].plot(this_bin_centers, nc1_arr_thisbin, c=args.color_arr[0], alpha=args.alpha)
                axes[index].plot(this_bin_centers, nc2_arr_thisbin, c=args.color_arr[1], alpha=args.alpha)
                axes[index].plot(this_bin_centers, nstars_arr_thisbin, c=args.color_arr[2], alpha=args.alpha)

            # -----decorating axes tick labels etc--------
            axes[index].tick_params(axis='both', which='major', labelsize=args.fontsize)

        # -------plotting bin edges-----------------
        if args.plot_bin_edges:
            for edge in base_bin_edges: axes[index].axvline(edge, ls='dashed', lw=0.5, c='grey')

        # -----plotting histograms from variants of this base bin--------
        axes[index] = collapse_and_plot_histogram(nc1_and_bins_arr, base_bin_edges, axes[index], color=args.color_arr[0], alpha=args.alpha, lw=lw)
        axes[index] = collapse_and_plot_histogram(nc2_and_bins_arr, base_bin_edges, axes[index], color=args.color_arr[1], alpha=args.alpha, lw=lw)
        axes[index] = collapse_and_plot_histogram(nstars_and_bins_arr, base_bin_edges, axes[index], color=args.color_arr[2], alpha=args.alpha, lw=lw)

        # -------annotate figure and save---------------
        axes[0].text(0.99, 0.98, 'Clusters', color=args.color_arr[0], ha='right', va='top', fontsize=args.fontsize, transform=axes[0].transAxes)
        axes[0].text(0.99, 0.85, 'Substructures', color=args.color_arr[1], ha='right', va='top', fontsize=args.fontsize, transform=axes[0].transAxes)
        axes[0].text(0.99, 0.70, 'Total stars', color=args.color_arr[2], ha='right', va='top', fontsize=args.fontsize, transform=axes[0].transAxes)

        if index == len(all_base_bins) - 1: axes[index].set_xlabel('Age (Myr)', fontsize=args.fontsize)

        axes[index].set_xlim(0, args.age_max * 1.05)
        if args.normalise:
            axes[index].set_ylim(0, 1.1)
            axes[index].set_ylabel('Normalised #', fontsize=args.fontsize)
        else:
            axes[index].set_yscale('log')
            axes[index].set_ylabel('Counts', fontsize=args.fontsize)

    fig.subplots_adjust(left=0.11, bottom=0.1, top=0.98, right=0.98, hspace=0.05, wspace=0.0)

    figname = args.output_dir / f'{args.input_file.split("_")[0]}_clustering_histograms.png'
    fig.savefig(figname)
    print(f'Saved figure at {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
