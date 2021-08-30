"""
Contains functions to perform the clustering analysis
"""
import math as m
from collections import defaultdict

import numpy as np
from scipy.signal import argrelextrema
from intervaltree import IntervalTree


def find_locals(smooth_tree, concat_regions_d):
    """
    Find local maximum and minimum of a smoothing curve

    Args:
        smooth_tree (IntervalTree): Interval are genomic regions, data np.array of smoothing score by position.
        concat_regions_d (dict): keys are start genomic regions, values are positions (index) relative to the start

    Returns:
        index_tree (IntervalTree): data are lists of tuples containing information from local maximum and minimum
    """
    index_tree = IntervalTree()
    reverse_cds_d = IntervalTree()

    # Reverse dictionary into IntervalTree
    if concat_regions_d:
        for genomic, cds in concat_regions_d.items():
            reverse_cds_d.addi(cds[0], cds[1] + 1, genomic)  # + 1 because end not included

    for interval in smooth_tree:
        indexes_info = []
        smooth = interval.data
        # Greater or equal
        max_eq = argrelextrema(smooth, np.greater_equal, order=1, mode='clip')[0].tolist()
        # Smaller or equal
        min_eq = argrelextrema(smooth, np.less_equal, order=1, mode='clip')[0].tolist()
        # Find maximums and minimums
        indexes = list(set(max_eq).symmetric_difference(set(min_eq)))
        # Add information to index
        for index in sorted(indexes):
            score_index = smooth[index]
            local = 0 if index in min_eq else 1  # 0 if minimum, 1 if maximum
            if not concat_regions_d:
                genomic = interval.begin + index
            else:
                for info in reverse_cds_d[index]:
                    genomic = info.data + index - info[0]
            indexes_info.append((local, index, genomic, score_index))

        # # Double check that there are no contiguous maximum positions
        # maximums = [i for i in indexes_info if i[0] == 1]
        # for i, maximum in enumerate(maximums.copy()[:-1]):
        #     left_max = maximum
        #     right_max = maximums[i + 1]
        #     if left_max[1] + 1 == right_max[1]:
        #         indexes_info.remove(left_max)

        # Add to new interval tree
        index_tree.addi(interval[0], interval[1], indexes_info)

    return index_tree


def find_clusters(index_tree):
    """
    Define a root cluster for each smoothing maximum

    Args:
        index_tree (IntervalTree): data are lists of tuples of 4 elements (min or max, cds region, genomic position,
                       smoothing score).

    Returns:
        clusters_tree (IntervalTree): data are dict of dict
    """

    clusters_tree = IntervalTree()

    for interval in index_tree:
        clusters = defaultdict(dict)
        j = 0
        indexes = interval.data

        # Iterate through all maximum and generate a cluster per maximum
        generator_maxs = (i for i in indexes if i[0] == 1)
        for maximum in generator_maxs:
            i = indexes.index(maximum)
            # Add maximum
            clusters[j]['max'] = (maximum[1], maximum[2], maximum[3])
            # Add margins
            # if maximum not in first nor last position
            if i != 0 and i != len(indexes) - 1:
                # if no contiguous left max
                if indexes[i - 1][0] != 1:
                    clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
                else:
                    clusters[j]['left_m'] = (maximum[1], maximum[2], maximum[3])
                # if no contiguous right max
                if indexes[i + 1][0] != 1:
                    clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])
                else:
                    clusters[j]['right_m'] = (maximum[1], maximum[2], maximum[3])
            # if first position
            elif i == 0:
                clusters[j]['left_m'] = (maximum[1], maximum[2], maximum[3])
                clusters[j]['right_m'] = (indexes[i + 1][1], indexes[i + 1][2], indexes[i + 1][3])
            # if last position
            else:
                clusters[j]['left_m'] = (indexes[i - 1][1], indexes[i - 1][2], indexes[i - 1][3])
                clusters[j]['right_m'] = (maximum[1], maximum[2], maximum[3])
            j += 1
        clusters_tree.addi(interval[0], interval[1], clusters)

    return clusters_tree


def merge(clusters_tree, window):
    """
    Iterate through clusters and merge them if they are closer than a given length.

    Args:
        clusters_tree (IntervalTree): genomic regions are intervals, data are clusters (dict of dict)
        window (int): clustering window

    Returns:
        merged_clusters_tree (IntervalTree): genomic regions are intervals, data are merged clusters (dict of dict)
    """
    merged_clusters_tree = IntervalTree()
    missed_clusters = defaultdict(dict)

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = interval.data.copy()
        n_clusters = len(clusters.keys())
        # Get the maximums to iterate
        maxs = []
        for c, v in clusters.items():
            maxs.append(v['max'][0])
        maxs_set = set(maxs)

        # Start iterations
        iterate = 1
        while iterate != 0:  # Iterate until no clusters updates occur
            stop = 0
            for x in range(n_clusters):
                # When x is a key in clusters
                if x in clusters.keys():
                    maximum = clusters[x]['max']
                    left_margin = clusters[x]['left_m']
                    right_margin = clusters[x]['right_m']
                    # Define the interval of search
                    search_r = set(range(right_margin[0], right_margin[0] + window + 1))    # +1, range()
                    # If there is a maximum in the search region
                    if search_r.intersection(maxs_set):
                        # Analyze only the closest cluster
                        intersect_clusters = sorted(list(search_r.intersection(maxs_set)))
                        if maxs.index(intersect_clusters[0]) == x:
                            if len(intersect_clusters) > 1:
                                intersect_cluster = maxs.index(intersect_clusters[1])
                            else:
                                intersect_cluster = maxs.index(intersect_clusters[0])
                        else:
                            intersect_cluster = maxs.index(intersect_clusters[0])
                        stop = 1
                        # When the testing max is greater than the intersected max,
                        # expand the right border and delete intersected cluster from clusters
                        if maximum[2] > clusters[intersect_cluster]['max'][2]:
                            clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                            del clusters[intersect_cluster]
                            maxs_set.remove(maxs[intersect_cluster])
                        # When the testing max is smaller than the intersected max,
                        # expand the left border of intersected max and remove the testing cluster (i)
                        elif maximum[2] < clusters[intersect_cluster]['max'][2]:
                            clusters[intersect_cluster]['left_m'] = left_margin
                            del clusters[x]
                            maxs_set.remove(maxs[x])
                        # When the testing max is equal than the intersected max
                        else:
                            # If contiguous maximums' positions, choose the first maximum
                            if maximum[0] == (clusters[intersect_cluster]['max'][0] - 1):
                                # Expand the right border and delete intersected cluster from clusters
                                clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                                del clusters[intersect_cluster]
                                maxs_set.remove(maxs[intersect_cluster])
                            # If not contiguous maximums' positions
                            elif maximum[0] < (clusters[intersect_cluster]['max'][0]):
                                # If contiguous clusters (right and left margings are equal), merge
                                if clusters[x]['right_m'] == clusters[intersect_cluster]['left_m']:
                                    # Expand the right border and delete intersected cluster from clusters
                                    clusters[x]['right_m'] = clusters[intersect_cluster]['right_m']
                                    del clusters[intersect_cluster]
                                    maxs_set.remove(maxs[intersect_cluster])
                                # Do not merge
                                else:
                                    missed_clusters[x] = clusters[x]
                                    del clusters[x]
                                    maxs_set.remove(maxs[x])
                            # If same cluster
                            else:
                                # Do not iterate again
                                maxs_set.remove(maxs[x])
            if stop == 0:
                iterate = 0

        # Merge clusters with overlapping margins
        # Iterate through all sorted clusters except the last one
        sorted_clusters_keys = sorted(clusters.keys())
        removed_clusters = set()
        if len(sorted_clusters_keys) >= 2:
            for i, cluster in enumerate(sorted_clusters_keys[:-1]):
                if cluster not in removed_clusters:
                    next_cluster = sorted_clusters_keys[i + 1]
                    # If there is overlap between right and left borders (they are equal)
                    # merge clusters
                    try:
                        a = clusters[cluster]['right_m']
                    except KeyError as e:
                        print('KeyError')
                        print(cluster)
                        print(clusters[cluster])

                    if clusters[cluster]['right_m'] == clusters[next_cluster]['left_m']:
                        clusters[cluster]['right_m'] = clusters[next_cluster]['right_m']
                        # maximum is the highest maximum, or, if they are egual, the first maximum
                        maximum = sorted(
                            [clusters[cluster]['max'], clusters[next_cluster]['max']], key=lambda x: x[2], reverse=True)[0]
                        clusters[cluster]['max'] = maximum
                        del clusters[next_cluster]
                        removed_clusters.add(next_cluster)

        for k, v in missed_clusters.items():
            clusters[k] = v

        merged_clusters_tree.addi(interval[0], interval[1], clusters)
        missed_clusters = defaultdict(dict)

    return merged_clusters_tree


def mapmut_and_filter(clusters_tree, mutations_in, cluster_mutations_cutoff):
    """
    Get the number of mutations within a cluster, remove those clusters below cutoff mutations

    Args:
        clusters_tree (IntervalTree): genomic regions are intervals, data are merged clusters (dict of dict)
        mutations_in (list): list of mutations fitting in regions
        cluster_mutations_cutoff (int): number of cluster mutations cutoff

    Returns:
        filter_clusters_tree (IntervalTree): genomic regions are intervals, data are filtered clusters (dict of dict)
    """
    filter_clusters_tree = IntervalTree()

    # Iterate through all regions
    for interval in clusters_tree:
        clusters = interval.data.copy()
        for cluster, values in interval.data.items():
            left = values['left_m'][1]
            right = values['right_m'][1]
            # Search mutations
            cluster_muts = [i for i in mutations_in if left <= i.position <= right]
            cluster_samples = set()
            for mut in cluster_muts:
                sample = mut.sample
                cluster_samples.add(sample)
            if len(cluster_muts) >= cluster_mutations_cutoff:
                clusters[cluster]['mutations'] = cluster_muts
                clusters[cluster]['samples'] = cluster_samples
                clusters[cluster]['fra_uniq_samples'] = len(cluster_samples)/len(cluster_muts)
            else:
                del clusters[cluster]
        filter_clusters_tree.addi(interval[0], interval[1], clusters)

    return filter_clusters_tree


def trim(clusters_tree, concat_regions_d):
    """
    Trim artificial clusters margins generated by KDE smooth. Trim margins to be mutated positions. Clusters left and
    right margins are updated.

    Args:
        clusters_tree (IntervalTree): genomic regions are intervals, data are filtered clusters (dict of dict)
        concat_regions_d (dict): keys are start genomic regions, values are positions (index) relative to the start

    Returns:
        trim_clusters_tree (IntervalTree): genomic regions are intervals, data are trimmed clusters (dict of dict)

    """
    trim_clusters_tree = IntervalTree()

    for interval in clusters_tree:
        clusters = interval.data.copy()
        for cluster, values in clusters.items():
            # Find new margins
            mutation_left = min(values['mutations'], key=lambda x: x.position)
            mutation_right = max(values['mutations'], key=lambda x: x.position)
            # Get index of mutation in cds
            if concat_regions_d:
                left_correction = concat_regions_d[mutation_left.region[0]].start
                right_correction = concat_regions_d[mutation_right.region[0]].start
            else:
                left_correction = 0
                right_correction = 0
            mutation_left_index = (mutation_left.position - mutation_left.region[0]) + left_correction
            mutation_right_index = (mutation_right.position - mutation_right.region[0]) + right_correction
            # Update clusters coordinates
            values['left_m'] = (mutation_left_index, mutation_left.position, np.nan)
            values['right_m'] = (mutation_right_index, mutation_right.position, np.nan)

        trim_clusters_tree.addi(interval[0], interval[1], clusters)

    return trim_clusters_tree


def score(clusters_tree, regions, mutations_element):
    """
    Score clusters with fraction of mutations formula and number of cluster's mutations

    Args:
        clusters_tree( IntervalTree): genomic regions are intervals, data are trimmed clusters (dict of dict)
        regions (IntervalTree): IntervalTree where intervals are genomic positions of an element
        mutations_element (int): number of mutations in the element

    Returns:
        score_clusters_tree (IntervalTree): genomic regions are intervals, data are scored clusters (dict of dict)
    """
    score_clusters_tree = IntervalTree()
    root = m.sqrt(2)

    for interval in clusters_tree:
        clusters = interval.data.copy()
        for cluster, values in clusters.items():
            score_ = 0
            mutated_positions_d = defaultdict(int)
            # Get number of mutations on each mutated position
            for mutation in values['mutations']:
                mutated_positions_d[mutation.position] += 1
            # Map mutated position and smoothing maximum to region
            for position, count in mutated_positions_d.items():
                map_mut_pos = set()
                map_smo_max = set()
                if regions[position]:
                    for i in regions[position]:
                        map_mut_pos = i
                    for i in regions[values['max'][1]]:
                        map_smo_max = i
                    # Calculate distance of position to smoothing maximum
                    if map_mut_pos[0] == map_smo_max[0]:
                        distance_to_max = abs(position - values['max'][1])
                    elif map_mut_pos[0] < map_smo_max[0]:
                        distance_to_max = (map_mut_pos[1] - position) + (values['max'][1] - map_smo_max[0])
                    else:
                        distance_to_max = (map_smo_max[1] - values['max'][1]) + (position - map_mut_pos[0])
                    # Calculate fraction of mutations
                    numerator = (count / mutations_element) * 100
                    # Calculate cluster score
                    denominator = m.pow(root, distance_to_max)
                    score_ += (numerator / denominator)
            # Update
            clusters[cluster]['score'] = score_ * len(values['mutations'])
        score_clusters_tree.addi(interval[0], interval[1], clusters)

    return score_clusters_tree
