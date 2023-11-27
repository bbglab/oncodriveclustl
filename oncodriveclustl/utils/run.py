"""
Contains the main class of the method
"""
import os.path
from concurrent.futures import ProcessPoolExecutor as Pool
import math
from collections import defaultdict
from collections import namedtuple

import bgreference as bgr
import daiquiri
from intervaltree import IntervalTree
import numpy as np
import pickle
from tqdm import tqdm

from oncodriveclustl.utils import smoothing as smo
from oncodriveclustl.utils import clustering as clu
from oncodriveclustl.utils import score
from oncodriveclustl.utils import analyticalpval as ap
from oncodriveclustl.utils.mutability import Mutabilities


# Logger
logger = daiquiri.getLogger()
Mutation = namedtuple('Mutation', 'position, region, alt, sample, group')
Cds = namedtuple('Cds', 'start, end')


class Experiment:
    """Class to analyze genomic elements of a cancer dataset"""

    def __init__(self,
                 regions_d,
                 concat_regions_d,
                 chromosomes_d,
                 strands_d,
                 mutations_d,
                 samples_d,
                 genome,
                 groups_d,
                 path_pickle,
                 mutability,
                 mutability_config,
                 element_mutations_cutoff,
                 cluster_mutations_cutoff,
                 smooth_window,
                 cluster_window,
                 kmer,
                 n_simulations,
                 simulation_mode,
                 simulation_window,
                 cores,
                 is_plot,
                 seed):

        """Initialize the Experiment class

        Args:
            regions_d (dict): dict, dictionary of IntervalTrees containing genomic regions from all analyzed elements
            concat_regions_d (dict): dictionary of dictionaries with positions of genomic regions relative to their
                start if analysis in concatenate mode
            chromosomes_d (dict): dictionary of elements (keys) and chromosomes (values)
            strands_d (dict): dictionary of elements (keys) and strands (values)
            mutations_d (dict): dictionary of elements (keys) and list of mutations formatted as namedtuple (values)
            samples_d (dict): dictionary of samples (keys) and number of mutations per sample (values)
            genome (str): genome to use
            groups_d (dict): dictionary of elements (keys) and groups (values)
            path_pickle (str): path to signature file (cache)
            element_mutations_cutoff (int): minimum number of mutations per genomic element to undertake analysis
            cluster_mutations_cutoff (int): minimum number of mutations to define a cluster
            smooth_window (int): Tukey kernel smoothing window length
            cluster_window (int): clustering window length
            kmer (int): context nucleotides to calculate the mutational probabilities (trinucleotides or
                pentanucleotides)
            n_simulations (int): number of simulations
            simulation_mode (str): simulation mode
            simulation_window (int): window length to simulate mutations
            cores (int): number of CPUs to use
            is_plot (bool): True generates extra data to build the cluster plot of an element
            seed (int): seed to randomize

        Returns:
            None

        """

        global Mutation
        global Cds

        self.regions_d = regions_d
        self.concat_regions_d = concat_regions_d
        self.chromosomes_d = chromosomes_d
        self.strands_d = strands_d
        self.mutations_d = mutations_d
        self.samples_d = samples_d
        self.genome = genome
        self.path_pickle = path_pickle
        self.groups_d = groups_d
        self.mutability = mutability
        self.mutability_config = mutability_config
        self.element_mutations_cutoff = element_mutations_cutoff
        self.cluster_mutations_cutoff = cluster_mutations_cutoff
        self.smooth_window = smooth_window + (1 - smooth_window % 2)
        # Calculate tukey filter
        self.tukey_filter = self.tukey(self.smooth_window)
        self.cluster_window = cluster_window
        self.kmer = kmer
        self.n_simulations = n_simulations
        self.simulation_mode = simulation_mode
        self.simulation_window = simulation_window + (1 - simulation_window % 2)
        self.cores = cores
        self.is_plot = is_plot
        self.main_seed = seed


        logger.warning(path_pickle)

        # Read CGC
        if self.genome in ['hg38', 'hg19']:
            # TODO: should we update this?
            with open(os.path.join(os.path.dirname(__file__), '../data/CGC_all_Oct15_10_29_09_2018.tsv'), 'r') as fd:
                next(fd)
                self.cgc_genes = set([line.split('\t')[0] for line in fd])
        else:
            self.cgc_genes = set()

    @staticmethod
    def tukey(window):
        """Tukey smoothing function generates tukey_filter for smoothing

        Args:
            window (int): smoothing window

        Returns:
            filter_ (np.array): tukey filter. Positions sum to 1

        """
        half_window = window // 2
        tukey = lambda x: np.maximum(1 - x ** 2, 0) ** 2
        filter_ = tukey(np.arange(-half_window, half_window + 1) / (half_window + 1))
        filter_ = filter_ / sum(filter_)
        return filter_

    @staticmethod
    def normalize(element, probs):
        """
        Given an array of probabilities, normalize them to 1

        Args:
            element (str): genomic element of analysis
            probs (list): list of probabilities

        Returns:
            normalized_probabilities (list): normalized probabilities

        """

        if sum(probs) != 0:
            prob_factor = 1 / sum(probs)
            normalized_probabilities = [prob_factor * p for p in probs]
        else:
            logger.error('No mutational probabilities derived from signatures in element {}'.format(element))
            # TODO: revise this error
            # Whenever this is raised there is an error downstream
            # Traceback (most recent call last):
            # File "/home/fcalvet/bin/mambaforge/envs/oncodriveclustl/lib/python3.10/concurrent/futures/process.py", line 246, in _process_worker
            #     r = call_item.fn(*call_item.args, **call_item.kwargs)
            # File "/home/fcalvet/bin/mambaforge/envs/oncodriveclustl/lib/python3.10/concurrent/futures/process.py", line 205, in _process_chunk
            #     return [fn(*args) for args in chunk]
            # File "/home/fcalvet/bin/mambaforge/envs/oncodriveclustl/lib/python3.10/concurrent/futures/process.py", line 205, in <listcomp>
            #     return [fn(*args) for args in chunk]
            # File "/home/fcalvet/projects/oncodriveclustl_dev/oncodriveclustl/oncodriveclustl/utils/run.py", line 400, in simulate_and_analysis
            #     simulations = np.random.choice(range(start_index, end_index), size=n_sim,
            # File "numpy/random/mtrand.pyx", line 953, in numpy.random.mtrand.RandomState.choice
            # TypeError: object of type 'bool' has no len()
            normalized_probabilities = False

        return normalized_probabilities

    def mut_probabilities(self, element):
        """
        Generate mutational probabilities per position of an element using the sequence context observed mutational
        probabilities calculated from the input cohort/s.

        Args:
            element (str): element to calculate pre-smoothing

        Returns:
            probs_tree (IntervalTree): IntervalTree of genomic regions. Length == 3*(genomic + simulation window)
            skip (bool): if True skip further analysis

        """
        nucleot = {'A', 'C', 'G', 'T'}
        probs_tree = defaultdict(IntervalTree)
        delta = 1 if self.kmer == 3 else 2
        half_window = self.simulation_window // 2
        skip = False

        # Check signatures dictionaries per group
        signatures_d = defaultdict()
        for group in self.groups_d[element]:
            # logger.info(f"{group}")
            if os.path.isfile(self.path_pickle):
                signature = pickle.load(open(self.path_pickle, "rb"))
                try:
                    signatures_d[group] = signature[group]
                except KeyError:
                    raise Exception('Signatures for group {} are missing in signatures dictionary\n'
                                    'Please check signatures file {}'.format(group, self.path_pickle))
            else:
                skip = True

        if not skip:
            # Iterate through genomic regions to get their sequences
            sequence = ''
            for interval in self.regions_d[element]:
                # logger.info(f"{element} {interval}")
                # logger.info(element)
                probabilities = defaultdict(list)
                expected_length = interval[1] - interval[0] + half_window*2
                start = interval[0] - half_window - delta
                size = interval[1] - interval[0] + half_window*2 + delta*2
                try:
                    sequence = bgr.refseq(self.genome, self.chromosomes_d[element], start, size)
                except ValueError as e:
                    logger.error(e, element, start, size, interval[0], interval[1])

                if sequence:
                    if self.mutability:
                        self.mutability_dict = Mutabilities(element, self.chromosomes_d[element], start, start + size, self.mutability_config).mutabilities_by_pos
                    # Search kmer probabilities
                    for n in range(delta, len(sequence)-delta):  # start to end
                        ref_kmer = sequence[n - delta: n + delta + 1]
                        prob = defaultdict(list)
                        if ref_kmer.count('N') == 0:
                            ref = ref_kmer[self.kmer//2]
                            alts = sorted(list(nucleot.difference({ref})))
                            if self.mutability:
                                # calculate mutational prob to any other kmer
                                # sort alternates to keep track                                
                                for alt in alts: # this is to iterate over the different possible mutations
                                    for group, signature in signatures_d.items():
                                        # logger.info(f"{start + n}\t{ref}\t{alt}\t{self.mutability_dict[start + n].get((ref, alt), 0)}")

                                        # TODO: pseudocount added to ensure that no position included in
                                        #  the regions under study has mutation probability equal to 0
                                        prob[group].append(self.mutability_dict[start + n].get((ref, alt), 0) + 0.5) 

                            else:
                                # calculate mutational prob to any other kmer
                                # sort alternates to keep track
                                for alt in alts: # this is to iterate over the different possible mutations
                                    for group, signature in signatures_d.items():
                                        prob[group].append(signature.get('{}>{}'.format(ref_kmer, alt), 0))
                        else:
                            logger.warning('Mutational probabilities for position {0} could not be calculated. '
                                        'Reverting {0}>ALT probabilities to 0'.format(n))
                            for group, signature in signatures_d.items():
                                prob[group].extend([0, 0, 0])
                        # Extend position info
                        for group in signatures_d.keys():
                            probabilities[group].extend(prob[group])
                    
                    # logger.info(probabilities)
                    # logger.info(f'{element}\t{start}\t{len(probabilities["mutations"])}')
                    # logger.info("end\n\n")

                    # Check and add
                    for group in signatures_d.keys():
                        if sum(probabilities[group]) != 0 and len(probabilities[group]) == 3 * expected_length:
                            probs_tree[group].addi(interval[0], interval[1], probabilities[group])
                        elif sum(probabilities[group]) == 0:
                            logger.critical('Context based mutational probabilities in {} '
                                            'region {}-{} equal to 0\n'.format(element, interval[0], interval[1]))
                            skip = True
                            break
                        elif len(probabilities[group]) != 3 * expected_length:
                            logger.warning('{} probabilities list length is different than expected'.format(element))
                            skip = True
                            break
                    if skip:
                        break
                else:
                    skip = True
                    break
        if skip:
            logger.critical('Context based mutational probabilities could not be calculated for {0}\n'
                            '{0} analysis is skipped'.format(element))
            
        # logger.info(f"{element} end")

        return probs_tree, skip

    def analysis(self, element, mutations, analysis_mode='sim', is_plot=False):
        """
        Run clustering analysis for observed or simulated mutations of an element

        Args:
            element (str): element of analysis.
            mutations (list): list of mutations of an element
            analysis_mode (str): observed or simulated; default simulated
            is_plot (bool): True returns smoothing array

        Returns:
            score_clusters_tree (IntervalTree): IntervalTree with scored clusters
            element_score (float): element score

        """
        if self.concat_regions_d:
            concat_regions_d = self.concat_regions_d[element]
        else:
            concat_regions_d = {}

        smooth_tree, mutations_in = smo.smooth_nucleotide(
            self.regions_d[element],
            concat_regions_d,
            mutations,
            self.tukey_filter,
            self.simulation_window
        )
        index_tree = clu.find_locals(smooth_tree, concat_regions_d)
        raw_clusters_tree = clu.find_clusters(index_tree)
        merge_clusters_tree = clu.merge(raw_clusters_tree, self.cluster_window)
        filter_clusters_tree = clu.mapmut_and_filter(merge_clusters_tree, mutations_in, self.cluster_mutations_cutoff)
        trim_clusters_tree = clu.trim(filter_clusters_tree, concat_regions_d)
        score_clusters_tree = clu.score(trim_clusters_tree, self.regions_d[element], len(mutations_in))
        logger.debug('Clusters scores calculated')
        element_score = score.element_score(score_clusters_tree, analysis_mode)
        logger.debug('Element score calculated')

        if is_plot:
            return smooth_tree
        else:
            return score_clusters_tree, element_score

    def simulate_and_analysis(self, item):
        """
        Simulate mutations and analyze simulations

        Args:
            item (tuple): element of analysis data containing element (str), probs_tree (IntervalTree),
                conseq_tree (IntervalTree), n_sim (int)
        Returns:
            element (str): element of analysis
            sim_scores_chunk (list): simulated element's results
            sim_cluster_chunk (list): simulated cluster's results
        """
        element, probs_tree, n_sim, seed = item
        sim_scores_chunk = []
        sim_cluster_chunk = []
        df_simulated_mutations = []
        half_window = (self.simulation_window - 1) // 2
        nucleot = {'A', 'C', 'G', 'T'}
        np.random.seed(seed)

        # Simulate mutations
        for mutation in self.mutations_d[element]:

            # Get coordinates of randomization window
            expected_hotspot_begin = mutation.position - half_window
            expected_hotspot_end = mutation.position + half_window

            if self.simulation_mode == 'region_restricted':
                """
                Region restricted mode samples simulated mutations in a window of length l that fits in the genomic 
                element. 
                
                First, it checks that the genomic region where the mutation is going to be simulated is longer or equal
                than l.            
                
                If this is true, it calculates the expected start and end positions of the simulation window. 
                If one of them falls outside the genomic element, the window of length l is displaced to fit in the 
                genomic element. If both expected start and end positions fall outside the genomic region, the 
                simulation window is trimmed and simulations are performed inside the genomic region. 
                
                If the genomic region is smaller than l, the simulation window becomes the genomic region. This means 
                that the simulation window is trimmed and simulations are performed between the end and the start of 
                the genomic region. 
                """

                if (mutation.region[1] - mutation.region[0]) >= self.simulation_window:
                    # Check if hospot outside region
                    check_5 = expected_hotspot_begin < mutation.region[0]
                    check_3 = expected_hotspot_end > (mutation.region[1] - 1)

                    if check_5 and check_3:
                        hotspot_begin = mutation.region[0]
                        hotspot_end = mutation.region[1] - 1  # regions end +1 in tree
                    elif check_5:
                        hotspot_begin = mutation.region[0]
                        hotspot_end = mutation.region[0] + self.simulation_window - 1  # window //2 per side
                    elif check_3:
                        hotspot_end = mutation.region[1] - 1  # regions end +1 in tree
                        hotspot_begin = hotspot_end - self.simulation_window + 1  # window //2 per side
                    else:
                        hotspot_begin = expected_hotspot_begin
                        hotspot_end = expected_hotspot_end
                else:
                    hotspot_begin = mutation.region[0]
                    hotspot_end = mutation.region[1] - 1  # regions end +1 in tree
            else:
                """
                Simulations are `mutation centered`, they are centered in the mutated position and can fall outside 
                the genomic region. 
                """
                hotspot_begin = expected_hotspot_begin
                hotspot_end = expected_hotspot_end

            # Map to index
            # 3* accounts for alternates in the array of probabilities
            # half_window added in probabilities array
            start_index = 3*(hotspot_begin - (mutation.region[0] - half_window))
            end_index = 3*(hotspot_end - (mutation.region[0] - half_window) + 1)  # +1, range and slice
            for interval in probs_tree[mutation.group][mutation.region[0]]:  # unique iteration
                # if element in ["VHL//VHL","PBRM1//PBRM1","MTOR//MTOR"] and sum(interval.data[start_index:end_index]) == 0 or \
                #                 ( mutation.region == (52563281, 52563496) ):
                #     logger.info(f"{mutation.group}\t{mutation.region}\t{start_index}-{end_index}\t{element}\t{interval.data[start_index:end_index]}")
                #     logger.info(f"{interval}")
                simulations = np.random.choice(range(start_index, end_index), size=n_sim,
                                                p=self.normalize(element, interval.data[start_index:end_index]))
                # Add info per simulated mutation
                list_simulations_per_mutation = []
                for count, index in enumerate(simulations):
                    position = mutation.region[0] - half_window + index // 3
                    ref_nucleotide = bgr.refseq(self.genome, self.chromosomes_d[element], position, 1)
                    # Calculate sorted alternates and obtain simulated alternated from index
                    if round(index / 3, 1) == (0.7 + index // 3):
                        alternate_index = 2
                    else:
                        alternate_index = 1 if round(index / 3, 1) == (0.3 + index // 3) else 0

                    alternate = sorted(list(nucleot.difference({ref_nucleotide})))[alternate_index]
                    # Simulated mutation
                    list_simulations_per_mutation.append(
                        Mutation(position, mutation.region, alternate, mutation.sample, mutation.group)
                    )
                df_simulated_mutations.append(list_simulations_per_mutation)

        # Start analysis
        logger.debug('Start analyzing simulations')
        for simulated_mutations in zip(*df_simulated_mutations):
            cutoff_clusters, element_score = self.analysis(element, list(simulated_mutations))
            sim_scores_chunk.append(element_score)
            for interval in cutoff_clusters:
                clusters = interval.data.copy()
                for cluster, values in clusters.items():
                    sim_cluster_chunk.append(values['score'])

        return element, sim_scores_chunk, sim_cluster_chunk

    @staticmethod
    def empirical_pvalue(observed, simulations):
        """
        Calculate empirical p-value using pseudocount (1)

        Args:
            observed (float): observed score
            simulations (list): simulated score

        Returns:
            float, p-value
        """
        return (len([x for x in simulations if x >= observed]) + 1) / len(simulations)

    def length(self, element):
        """Calculate length of an element (sum of input regions). If clustering in protein sequence, length in aa

        Args:
            element (str): element to analyze
        Returns:
            length (int): length of an element (bp or aa)
        """
        length = 0

        for region in self.regions_d[element]:
            length += region[1] - region[0]

        return length

    def post_process(self, items):
        """
        Post processing of results

        Args:
            items (list): list of tuples containing results for a chunk of elements

        Returns:
            results (list): list of results for a chunk of elements
        """

        pseudo_pvalue = 1.1102230246251566e-19
        results = []
        for item in items:
            element, obs_clusters, obs_score, sim_clusters_scores, sim_element_scores, seed = item
            np.random.seed(seed)
            mut_in_clusters = 0

            # Get GCG boolean
            cgc = element.split('//')[0] in self.cgc_genes

            # Calculate length
            element_length = self.length(element)

            # If analyzed element and element has observed clusters
            if type(obs_clusters) != float:
                if type(sim_element_scores) != float:

                    n_clusters_sim = len(sim_clusters_scores)

                    # If no simulated clusters, all simulated element scores are 0
                    if sum(sim_element_scores) == 0:

                        # Calculate p-value for element
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_element_scores)
                        obs_pvalue = obj.calculate(obs_score)

                        # Add p-value and number of mutations for each cluster
                        n_clusters = 0
                        for interval in obs_clusters:
                            for cluster, values in interval.data.items():
                                interval.data[cluster]['p'] = pseudo_pvalue
                                n_clusters += 1
                                mut_in_clusters += len(values['mutations'])

                        empirical_pvalue = self.empirical_pvalue(obs_score, sim_element_scores)
                        analytical_pvalue = top_cluster_pvalue = obs_pvalue

                    # Simulated clusters
                    else:
                        # Element score empirical p-value
                        empirical_pvalue = self.empirical_pvalue(obs_score, sim_element_scores)

                        # Element score analytical p-value
                        sim_scores_array_1000 = np.random.choice(sim_element_scores, size=1000, replace=False)
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_scores_array_1000)
                        analytical_pvalue = obj.calculate(obs_score)

                        # Cluster analytical p-values
                        obs_clusters_score = []

                        # Check how many simulated clusters are and add 0 if this number < number of simulations
                        # 0 means no simulated cluster
                        if n_clusters_sim < self.n_simulations:
                            missing_clusters = self.n_simulations - len(sim_clusters_scores)
                            sim_clusters_scores = sim_clusters_scores + [0] * missing_clusters

                        # Random choice 1000 simulated cluster scores
                        sim_clusters_scores = np.random.choice(sim_clusters_scores, size=1000, replace=False)

                        # Fit distribution
                        obj = ap.AnalyticalPvalue()
                        obj.calculate_bandwidth(sim_clusters_scores)
                        for interval in obs_clusters:
                            for cluster, values in interval.data.items():
                                cluster_p_value = obj.calculate(values['score'])
                                interval.data[cluster]['p'] = cluster_p_value
                                obs_clusters_score.append((values['score'], cluster_p_value))
                                mut_in_clusters += len(values['mutations'])

                        n_clusters = len(obs_clusters_score)

                        # Element top cluster analytical p-value
                        top_cluster_pvalue = sorted(obs_clusters_score, key=lambda k: (k[0], -k[1]), reverse=True)[0][1]
                        logger.debug('P-values calculated')

                else:
                    n_clusters = obs_score = mut_in_clusters = 0
                    n_clusters_sim = float('nan')
                    empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')
            else:
                n_clusters = n_clusters_sim = obs_score = mut_in_clusters = float('nan')
                empirical_pvalue = analytical_pvalue = top_cluster_pvalue = float('nan')

            if not self.is_plot:
                mutations = []
                smooth_tree = {}
                probs_tree = {}
            else:
                # Plot
                mutations = self.mutations_d[element]
                smooth_tree = self.analysis(element, self.mutations_d[element], analysis_mode='obs', is_plot=True)
                probs_tree, _ = self.mut_probabilities(element)

            results.append((
                element,
                (len(self.mutations_d[element]), mut_in_clusters, n_clusters, n_clusters_sim, obs_score,
                 empirical_pvalue, analytical_pvalue, top_cluster_pvalue),
                (obs_clusters, smooth_tree, probs_tree),
                (self.regions_d[element], self.chromosomes_d[element], self.strands_d[element],
                 mutations, element_length, cgc)
            ))

        return results

    def run(self):
        """
        Run clustering analysis of genomic elements

        Returns:
            elements_results (dict): dict of dict, results of elements analysis
            clusters_results (dict): dict of dict, results of clusters analysis
            global_info_results (dict): dict of dict, extra information needed to write/plot results
        """
        elements_results = defaultdict()
        clusters_results = defaultdict()
        global_info_results = defaultdict()
        np.random.seed(self.main_seed)
        cores_minus_one = self.cores - 1 if self.cores != 1 else self.cores

        # Filter elements >= cutoff mutations
        analyzed_elements = []
        belowcutoff_elements = []
        for elem in sorted(self.regions_d.keys()):
            if len(self.mutations_d[elem]) >= self.element_mutations_cutoff:
                analyzed_elements.append(elem)
            else:
                belowcutoff_elements.append(elem)
        logger.info('Calculating results {} element{}...'.format(len(analyzed_elements),
                    's' if len(analyzed_elements) > 1 else ''))

        # Performance tunning
        pf_num_elements = 100
        pf_num_simulations = 100000

        # Run analysis on each element
        analyzed_elements = list(chunkizator(analyzed_elements, pf_num_elements))
        chunk_seed_list = np.random.randint(2**32 - 1, size=len(analyzed_elements))
        for i, elements in enumerate(analyzed_elements, start=1):
            observed_clusters_d = {}
            observed_scores_d = {}
            simulations = []
            simulated_elements = []
            noprobabilities_elements = []
            nocluster_elements = []
            logger.info("Iteration {} of {}".format(i, len(analyzed_elements)))
            second_seed = chunk_seed_list[(i - 1)]

            # TODO
            # see if we can parallelize this loop
            # the option of adding mutabilities slows down oncodriveclustl
            # since it is sequentially loading the mutabilities of each element
            for element in elements:
                # Calculate observed results
                observed_clusters_d[element], observed_scores_d[element] = self.analysis(
                                                                                element,
                                                                                self.mutations_d[element],
                                                                                analysis_mode='obs'
                                                                                )
                logger.debug('Observed mutations analyzed')
                # Calculate simulated mutations results if element has observed clusters at least in one region
                check = 0
                for interval in observed_clusters_d[element]:
                    if interval.data:
                        check = 1
                    if check == 1:
                        break
                # If there is at least one region containing clusters
                if check != 0:
                    # Check if probabilities can be calculated
                    probs_tree, skip = self.mut_probabilities(element)
                    if not skip:
                        element_size = ((len(self.mutations_d[element]) * self.n_simulations) // pf_num_simulations) + 1
                        chunk_size = (self.n_simulations // element_size) + 1
                        simulations += [(element, probs_tree, n_sim, second_seed) for n_sim in partitions_list(
                                        self.n_simulations,
                                        chunk_size,
                                        )]
                        simulated_elements.append(element)
                    # Skip analysis if problems with mutational probabilities
                    else:
                        noprobabilities_elements.append(element)
                else:
                    nocluster_elements.append(element)

            with Pool(max_workers=self.cores) as executor:
                # Process simulations
                sim_scores_list = defaultdict(list)
                sim_clusters_list = defaultdict(list)

                # Analyze elements with observed clusters
                for element, sim_scores_chunk, sim_cluster_chunk in tqdm(
                        executor.map(self.simulate_and_analysis, simulations), total=len(simulations),
                        desc="simulations".rjust(30)):
                    sim_scores_list[element] += sim_scores_chunk
                    sim_clusters_list[element] += sim_cluster_chunk

                # Add information of elements without clusters
                for element in nocluster_elements:
                    sim_scores_list[element] = sim_clusters_list[element] = float('nan')

                # Post process
                post_item_simulated = [(e,
                                        observed_clusters_d[e],
                                        observed_scores_d[e],
                                        sim_clusters_list[e],
                                        sim_scores_list[e],
                                        second_seed) for e in simulated_elements + nocluster_elements]
                post_item_nan = [(e, float('nan'), float('nan'), float('nan'), float('nan'), None) for
                                 e in noprobabilities_elements + belowcutoff_elements]
                total_items_split = list(
                    chunkizator(post_item_simulated, int(math.ceil(len(post_item_simulated) / cores_minus_one)))
                )
                total_items_split.append(post_item_nan)
                for results in tqdm(executor.map(
                        self.post_process, total_items_split), total=self.cores, desc="post processing".rjust(30)):
                    for element, elem_res, clust_res, info in results:
                        elements_results[element] = elem_res
                        clusters_results[element] = clust_res
                        global_info_results[element] = info

        return elements_results, clusters_results, global_info_results


def partitions_list(total_size, chunk_size):
    """
    Create a list of values less or equal to chunk_size that sum total_size

    Args:
        total_size (int): total size
        chunk_size (int): chunk size

    Returns:
        partitions (list): list of integers

    """
    partitions = [chunk_size for _ in range(total_size // chunk_size)]

    res = total_size % chunk_size
    if res != 0:
        partitions += [res]

    return partitions


def chunkizator(iterable, size=1000):
    """
    Creates chunks from an iterable

    Args:
        iterable (list): total list of elements
        size (int): number of elements in the chunk

    Returns:
        chunk (list): list of elements

    """
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk
