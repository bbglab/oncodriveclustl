"""
Contains functions to generate a cluster plot
"""
import os
import math
from collections import defaultdict

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec


def concat_regions_to_plot(element, regions, strand):
    """
    Generate info to map regions into a concatenated sequence to plot

    Args:
        element (str): element of analysis
        regions (IntervalTree): intervaltree of genomic regions for an element
        strand (str): positive or negative strand

    Returns:
        concat_regions_d (dict): dictionary of genomic region start (keys) and corresponding concatenated sequence
            index of start and end
        plot_regions_xcoords (list): list of indexes to plot corresponding to start/end of genomic sequences.
            Indexes are reversed for negative strand elements
    """
    concat_regions_d = defaultdict(dict)
    plot_regions_xcoords = []

    start = 0
    for region in sorted(regions):
        length = region.end - region.begin
        end = start + length - 1
        concat_regions_d[element][region.begin] = (start, end)
        plot_regions_xcoords.append(start)
        start = end + 1
    plot_regions_xcoords.append((end + 1))

    element_length = plot_regions_xcoords[-1]

    if strand == '-':
        plot_regions_xcoords = list(map(lambda x: element_length - x, plot_regions_xcoords))

    return concat_regions_d, plot_regions_xcoords


def concat_smooth(smooth_tree, strand):
    """
    Concatenate smoothing arrays for each genomic region to plot

    Args:
        smooth_tree (IntervalTree): intervaltree of smoothing in genomic regions for an element
        strand (str): positive or negative strand

    Returns:
        smooth (list): smoothing scores per position in an element. List is reversed for negative strand elements
    """

    smooth = []
    for region in sorted(smooth_tree):
        smooth.extend(region.data)
    if strand == '-':
        smooth = smooth[::-1]

    return smooth


def mutations_index(element, mutations, concat_regions_d, strand, length):
    """
    Map mutations to concatenated sequence to plot

    Args:
        element (str): element under analysis
        mutations (list): list of mutations for an element
        concat_regions_d (dict): dictionary of genomic region start (keys) and corresponding concatenated sequence
            index of start and end
        strand (str): positive or negative strand
        length (int): length of the genomic element

    Returns:
        mutations_xcoords_number (dict): dictionary of mutated coordinates, keys are positions, values are number of
            mutations. Reversed coordinates if strand is negative

    """
    mutations_xcoords_number = defaultdict(int)

    for mutation in mutations:
        index = mutation.position - mutation.region[0] + concat_regions_d[element][mutation.region[0]][0]
        if strand == '-':
            index = length - index
        mutations_xcoords_number[index] += 1

    return mutations_xcoords_number


def cluster_coords(element, clusters_tree, strand, length, concat_regions_d):
    """
    Generate info to plot clusters

    Args:
        element (str): element under analysis
        clusters_tree (IntervalTree): intervaltree of observed clusters for an element
        strand (str): positive or negative strand
        length (int): length of the genomic element
        concat_regions_d (dict): dictionary of genomic region start (keys) and corresponding concatenated sequence
            index of start and end

    Returns:
        plot_cluster_xcoords (list): list of indexes to plot observed clusters.
            Indexes are reversed for negative strand elements
    """
    plot_cluster_xcoords = []
    for region in clusters_tree:
        for cluster, info in region.data.items():
            for cluster_margin in ['left_m', 'right_m']:
                index = info[cluster_margin][0] + concat_regions_d[element][region[0]][0]
                plot_cluster_xcoords.append(index)
    if strand == '-':
        plot_cluster_xcoords = list(map(lambda x: length - x, plot_cluster_xcoords))

    return plot_cluster_xcoords


def cluster_coords_sig(element, clusters_tree, strand, length, concat_regions_d, pvalue_threshold):
    """
    Generate info to plot clusters

    Args:
        element (str): element under analysis
        clusters_tree (IntervalTree): intervaltree of observed clusters for an element
        strand (str): positive or negative strand
        length (int): length of the genomic element
        concat_regions_d (dict): dictionary of genomic region start (keys) and corresponding concatenated sequence
            index of start and end

    Returns:
        plot_cluster_xcoords (list): list of indexes to plot observed clusters.
            Indexes are reversed for negative strand elements
    """
    plot_cluster_xcoords = []
    for region in clusters_tree:
        for cluster, info in region.data.items():
            if info['p'] < pvalue_threshold:
                for cluster_margin in ['left_m', 'right_m']:
                    index = info[cluster_margin][0] + concat_regions_d[element][region[0]][0]
                    plot_cluster_xcoords.append(index)
    if strand == '-':
        plot_cluster_xcoords = list(map(lambda x: length - x, plot_cluster_xcoords))

    return plot_cluster_xcoords




def clusters_plot(element,
                  plot_regions_xcoords,
                  smooth,
                  mutations_xcoords_number,
                  plot_cluster_xcoords,
                  output,
                  fig_height=4,
                  fig_width=8,
                  title=None
                  ):
    """

    Args:
        element (str): element under analysis
        plot_regions_xcoords (list): list of indexes to plot corresponding to start/end of genomic sequences.
            Indexes are reversed for negative strand elements
        smooth (list): smoothing scores per position in an element. List is reversed for negative strand elements
        mutations_xcoords_number (dict): dictionary of mutated coordinates, keys are positions,
            values are number of mutations.
            Reversed coordinates if strand is negative
        plot_cluster_xcoords (list): list of indexes to plot observed clusters.
            Indexes are reversed for negative strand elements
        output (str): path to output
        fig_height (int): figure height
        fig_width (int): figure width
        title (str): title

    Returns:
        None

    """

    fig, ax = plt.subplots(nrows=10, ncols=1)
    fig.set_size_inches(fig_width, fig_height)
    gs = gridspec.GridSpec(2, 1, height_ratios=[9, 1])
    gs.update(wspace=0.025, hspace=0.05)
    plt.rc('axes', linewidth=1.25)

    # Grid 1
    ax0 = plt.subplot(gs[0])
    ax0_2 = ax0.twinx()
    ax0.get_xaxis().set_visible(False)

    if title:
        ax0.set_title(title, fontsize=16)
    else:
        ax0.set_title(element.split('//')[0], fontsize=16)

    # Regions
    ax0.axhline(y=0, color='lightgrey', linestyle='--', linewidth=2)
    for index, region_xcoord in enumerate(plot_regions_xcoords):
        index += 1
        if not index % 2 and region_xcoord != plot_regions_xcoords[-1]:
            ax0.axvspan(xmin=region_xcoord, xmax=plot_regions_xcoords[index], color='grey', alpha=0.1)
        ax0.axvline(x=region_xcoord, color='grey', linestyle='--', linewidth=2, alpha=0.25)

    # Smoothing
    ax0_2.plot(smooth, color='#0571b0', linewidth=2, solid_capstyle='butt', alpha=0.75)
    ax0_2.get_yaxis().set_ticks([])

    # Mutations
    for xcoord, ycoord in mutations_xcoords_number.items():
        x = (xcoord, xcoord)
        y = (0, ycoord)
        ax0.plot(x, y, color='black', linewidth=2, alpha=0.75)
    ax0.set_ylabel('# mutations', color='black', fontsize=16)
    max_mutations = math.ceil(max(mutations_xcoords_number.values()))
    if max_mutations >= 30:
        yint = range(0, max_mutations + 1, 10)
    elif 30 > max_mutations >= 5:
        yint = range(0, max_mutations + 1, 5)
    else:
        yint = range(0, max_mutations + 1, 1)
    ax0.get_yaxis().set_ticks(yint)
    keys = list(map(lambda x: int(x), mutations_xcoords_number.keys()))
    values = list(map(int, mutations_xcoords_number.values()))
    ax0.plot(keys, values, linewidth=0, marker='o', markersize=4, color='black', alpha=1)
    ax0.set_axisbelow(True)

    # Grid 2
    ax1 = plt.subplot(gs[1], sharex=ax0)
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel('Clusters', color='black', fontsize=10)
    ax1.get_yaxis().set_ticks([])

    # Regions
    ax1.axhline(y=0, color='lightgrey', linestyle='--', linewidth=1)
    for index, region_xcoord in enumerate(plot_regions_xcoords):
        index += 1
        if not index % 2 and region_xcoord != plot_regions_xcoords[-1]:
            ax1.axvspan(xmin=region_xcoord, xmax=plot_regions_xcoords[index], color='grey', alpha=0.1)
        ax1.axvline(x=region_xcoord, color='grey', linestyle='--', linewidth=2, alpha=0.25)

    for index, cluster_xcoord in enumerate(plot_cluster_xcoords):
        index += 1
        if index % 2: # and cluster_xcoord != plot_cluster_xcoords[-1]:
            ax1.axvspan(xmin=cluster_xcoord, xmax=plot_cluster_xcoords[index],
                        color='#0571b0', alpha=0.75, linewidth=1.5)

    plt.savefig(output, bbox_inches='tight')


def make_clustplot(elements_results, clusters_results, global_info_results, directory):
    """

    Args:
        elements_results (dict): keys are elements, values are elements results
        clusters_results (dict): keys are elements, values are cluster results
        global_info_results (dict): keys are elements, values are global information about the element
        directory (str): path to output

    Returns:
        info (str): message

    """
    info_cluster_plots = []

    # Parse
    for element, data in elements_results.items():
        _, _, _, _, score, _, pvalue, _ = data
        observed_clusters, smooth_tree, probabilities = clusters_results[element]
        regions, _, strand, mutations, length, _ = global_info_results[element]

        if type(observed_clusters) != float:

            # Preprocess data
            concat_regions_d, plot_regions_xcoords = concat_regions_to_plot(element, regions, strand)
            smooth = concat_smooth(smooth_tree, strand)
            mutations_xcoords_number = mutations_index(element, mutations, concat_regions_d, strand, length)
            plot_cluster_xcoords = cluster_coords(element, observed_clusters, strand, length, concat_regions_d)
            output = os.path.join(directory, '{}_plot.png'.format(element.split('//')[0]))
            # Plot
            clusters_plot(element,
                          plot_regions_xcoords,
                          smooth,
                          mutations_xcoords_number,
                          plot_cluster_xcoords,
                          output,
                          )
            info = 'Cluster plot for {} generated at: {}'.format(element.split('//')[0], output)
        else:
            info = 'No clusters found in {}. No cluster plot is generated for this element'.format(element.split('//')[0])
        info_cluster_plots.append(info)

    return info_cluster_plots


def make_clustplot_sig(elements_results, clusters_results, global_info_results, directory):
    """

    Args:
        elements_results (dict): keys are elements, values are elements results
        clusters_results (dict): keys are elements, values are cluster results
        global_info_results (dict): keys are elements, values are global information about the element
        directory (str): path to output

    Returns:
        info (str): message

    """
    info_cluster_plots = []

    # Parse
    for element, data in elements_results.items():
        _, _, _, _, score, _, pvalue, _ = data
        observed_clusters, smooth_tree, probabilities = clusters_results[element]
        regions, _, strand, mutations, length, _ = global_info_results[element]

        if type(observed_clusters) != float:

            # Preprocess data
            concat_regions_d, plot_regions_xcoords = concat_regions_to_plot(element, regions, strand)
            smooth = concat_smooth(smooth_tree, strand)
            mutations_xcoords_number = mutations_index(element, mutations, concat_regions_d, strand, length)
            plot_cluster_xcoords = cluster_coords_sig(element, observed_clusters, strand, length, concat_regions_d, 0.05)
            if len(plot_cluster_xcoords) > 0:
                output = os.path.join(directory, '{}_plot.png'.format(element.split('//')[0]))
                # Plot
                clusters_plot(element,
                            plot_regions_xcoords,
                            smooth,
                            mutations_xcoords_number,
                            plot_cluster_xcoords,
                            output,
                            )
                info = 'Cluster plot for {} generated at: {}'.format(element.split('//')[0], output)
            else:
                info = 'No significant clusters found in {}. No cluster plot is generated for this element'.format(element.split('//')[0])
        else:
            info = 'No clusters found in {}. No cluster plot is generated for this element'.format(element.split('//')[0])
        info_cluster_plots.append(info)

    return info_cluster_plots

