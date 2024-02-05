"""
Contains the command line parsing
"""

# Import modules
import os
import logging

import bgsignature as bgsign
import click
import daiquiri
import pickle
import json

from oncodriveclustl.utils import cluster_plot as cplot
from oncodriveclustl.utils import exceptions as excep
from oncodriveclustl.utils import parsing as pars
from oncodriveclustl.utils import postprocessing as postp
from oncodriveclustl.utils import qq_plot as qplot
from oncodriveclustl.utils import run as exp
from oncodriveclustl.utils.mutability import init_mutabilities_module

# Global variables
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations')
@click.option('-r', '--regions-file', default=None, required=True, type=click.Path(exists=True),
              help='File with the genomic regions to analyze')
@click.option('-o', '--output-directory', default=None, required=True,
              help='Output directory to be created')
@click.option('-sig', '--input-signature', default=None, required=False, type=click.Path(exists=True),
              help='File containing input context based mutational probabilities (signature)')
@click.option('-ef', '--elements-file', default=None, type=click.Path(exists=True),
              help='File with the symbols of the elements to analyze')
@click.option('-e', '--elements', default=None, multiple=True, type=click.STRING,
              help='Symbol of the element(s) to analyze')
@click.option('-g', '--genome', default='hg19',
              type=click.Choice(['hg38', 'hg19', 'mm10', 'c3h', 'car', 'cast', 'f344']),
              help='Genome to use')
@click.option('-emut', '--element-mutations', type=click.INT, default=2,
              help='Cutoff of element mutations. Default is 2')
@click.option('-cmut', '--cluster-mutations', type=click.INT, default=2,
              help='Cutoff of cluster mutations. Default is 2')
@click.option('-sw', '--smooth-window', type=click.IntRange(3, 101), default=11,
              help='Smoothing window. Default is 11')
@click.option('-cw', '--cluster-window', type=click.IntRange(3, 101), default=11,
              help='Cluster window. Default is 11')
@click.option('-kmer', '--kmer', default='3', help='K-mer nucleotide context',
              type=click.Choice(['3', '5']))
@click.option('-n', '--n-simulations', type=click.INT, default=1000,
              help='number of simulations. Default is 1000')
@click.option('-sim', '--simulation-mode', default='mutation_centered', help='Simulation mode',
              type=click.Choice(['mutation_centered', 'region_restricted']))
@click.option('-simw', '--simulation-window', type=click.IntRange(19, 101), default=31,
              help='Simulation window. Default is 31')
@click.option('-sigcalc', '--signature-calculation', default='frequencies',
              help='Signature calculation: mutation frequencies (default) or k-mer mutation counts normalized by k-mer '
                   'region counts',
              type=click.Choice(['frequencies', 'region_normalized']))
@click.option('-siggroup', '--signature-group', default=None,
              help='Header of the column to group signatures calculation',
              type=click.Choice(['SIGNATURE', 'SAMPLE', 'CANCER_TYPE']))
@click.option('-mutab', '--mutabilities-path',default=None, required=False, type=click.Path(exists=True),
              help='Path to the mutabilities configuration file.')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it will use all the available cores.')
@click.option('--seed', type=click.INT, default=None, help='Seed to use in the simulations')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
@click.option('--concatenate', is_flag=True, help='Calculate clustering on concatenated genomic regions (e.g., exons '
                                                  'in coding sequences)')
@click.option('--clustplot', is_flag=True, help='Generate a needle plot with clusters for an element')
@click.option('--clustplotsig', is_flag=True, help='Generate a needle plot with the significant clusters for an element')
@click.option('--qqplot', is_flag=True, help='Generate a quantile-quantile (Q-Q) plot for a dataset')
@click.option('--gzip', is_flag=True, help='Gzip compress files')
def main(input_file,
         regions_file,
         output_directory,
         input_signature,
         elements_file,
         elements,
         genome,
         element_mutations,
         cluster_mutations,
         smooth_window,
         cluster_window,
         kmer,
         n_simulations,
         simulation_mode,
         simulation_window,
         signature_calculation,
         signature_group,
         mutabilities_path,
         cores,
         seed,
         log_level,
         concatenate,
         clustplot,
         clustplotsig,
         qqplot,
         gzip
         ):
    """
    OncodriveCLUSTL is a sequence based clustering method to identify cancer drivers across the genome

    Args:
        input_file (str): path to mutations file
        regions_file (str): path to input genomic coordinates file
        output_directory(str): path to output directory. Output files will be generated in it.
        input_signature (str): path to file containing input context based mutational probabilities.
            By default (when no input signatures), OncodriveCLUSTL will calculate them from the mutations input file.
        elements_file (str): path to file containing one element per row (optional) to analyzed the listed elements.
            By default, OncodriveCLUSTL analyzes all genomic elements contained in `regions_file`.
        elements (str): genomic element symbol (optional). The analysis will be performed only on the specified GEs.
        genome (str): genome to use: 'hg38', 'hg19', 'mm10', 'c3h', 'car', 'cast' and 'f344'
        element_mutations (int): minimum number of mutations per genomic element to undertake analysis
        cluster_mutations (int): minimum number of mutations to define a cluster
        smooth_window (int): Tukey kernel smoothing window length
        cluster_window (int): clustering window length
        kmer (int): context nucleotides to calculate the mutational probabilities (trinucleotides or pentanucleotides)
        n_simulations (int): number of simulations
        simulation_mode (str): simulation mode
        simulation_window (int): window length to simulate mutations
        signature_calculation (str): signature calculation, mutation frequencies (default) or mutation counts
            normalized by k-mer region counts
        signature_group (str): header of the column to group signatures. One signature will be computed for each group
        mutabilities_path (str): path to file containing the configuration of the mutabilities.
        cores (int): number of CPUs to use
        seed (int): seed
        log_level (str): verbosity of the logger
        concatenate (bool): flag to calculate clustering on collapsed genomic regions (e.g., coding regions in a gene)
        clustplot (bool): flag to generate a needle plot with clusters for an element
        qqplot (bool): flat to generate a quantile-quantile (QQ) plot for a dataset
        gzip (bool): flag to generate GZIP compressed output files

    Returns:
        None

    """

    global logger

    # Get output directory
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Get cache directory
    path_cache = os.path.join(output_directory, 'cache')
    os.makedirs(path_cache, exist_ok=True)

    # Get output files name
    if elements_file is not None:
        output_file = elements_file.split('/')[-1]
    else:
        output_file = 'results'
    if gzip:
        elements_output_file = 'elements_{}.txt.gz'.format(output_file)
        clusters_output_file = 'clusters_{}.tsv.gz'.format(output_file)
    else:
        elements_output_file = 'elements_{}.txt'.format(output_file)
        clusters_output_file = 'clusters_{}.tsv'.format(output_file)

    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(
            filename=output_file + '.log',
            directory=output_directory
        )
    ))
    logger = daiquiri.getLogger()

    # suppress log messages from some libraries
    daiquiri.getLogger('bgdata').setLevel(logging.WARNING)
    daiquiri.getLogger('bgsignature').setLevel(logging.WARNING)

    logger.info('OncodriveCLUSTL')
    logger.info('\n'.join([
        '',
        'input_file: {}'.format(input_file),
        'regions_file: {}'.format(regions_file),
        'input_signature: {}'.format(input_signature),
        'output_directory: {}'.format(output_directory),
        'genome: {}'.format(genome),
        'element_mutations: {}'.format(element_mutations),
        'cluster_mutations: {}'.format(cluster_mutations),
        'concatenate: {}'.format(concatenate),
        'smooth_window: {}'.format(smooth_window),
        'cluster_window: {}'.format(cluster_window),
        'k-mer: {}'.format(kmer),
        'simulation_mode: {}'.format(simulation_mode),
        'simulation_window: {}'.format(simulation_window),
        'n_simulations: {}'.format(n_simulations),
        'signature_calculation: {}'.format(signature_calculation),
        'signature_group: {}'.format(signature_group),
        'mutabilities_path: {}'.format(mutabilities_path),
        'cores: {}'.format(cores),
        'gzip: {}'.format(gzip),
        'seed: {}'.format(seed)
    ]))
    logger.info('Initializing OncodriveCLUSTL...')

    # Check parameters
    if simulation_window == 31 and smooth_window == 11 and cluster_window == 11:
        logger.warning(
            '\nRunning with default simulating, smoothing and clustering OncodriveCLUSTL parameters. '
            'Default parameters may not be optimal for your data.\n'
            'Please, read Supplementary Methods to perform model selection for your data.'
        )

    if mutabilities_path is not None:
        logger.info(f"Computing missense mut probabilities using mutabilities...")
        mutab_config = json.load(open(mutabilities_path))
        init_mutabilities_module(mutab_config)
        mutability = True
    else:
        mutab_config = None
        mutability = False

    if not input_signature and signature_calculation == 'frequencies':
        logger.warning(
            '\nSignatures will be calculated as mutation frequencies: '
            '# mutated ref>alt k-mer counts / # total substitutions\n'
            'Please, read Supplementary Methods to perform a more accurate signatures calculation'
        )

    if signature_calculation == 'region_normalized':
        logger.warning(
            '\nMutation k-mer counts will be normalized by k-mer region counts in {}\n'
            'Only mutations inside regions will contribute for the signature calculation'.format(regions_file)
        )

    if n_simulations < 1000:
        raise excep.UserInputError('Invalid number of simulations: please choose an integer greater than 1000')

    if clustplot:
        if len(elements) > 10:
            raise excep.UserInputError('Needle plots can only be generated for a maximum of 10 elements')

    # Create a list of elements to analyze
    if elements is not None:
        elements = set(elements)
    if elements_file is not None:
        elements |= set([line.strip().split()[0] for line in open(elements_file, 'r')])
    if elements is None and elements_file is None:
        elements = set([])
    if elements:
        logger.info(
            'Input element{}: {}'.format('s' if len(elements) > 1 else '', len(elements))
        )
        logger.info(', '.join(elements))

    # Parse regions and dataset mutations
    logger.info('Parsing genomic regions and mutations...')
    regions_d, concat_regions_d, chromosomes_d, strands_d, mutations_d, samples_d, groups_d = pars.parse(
        regions_file,
        elements,
        input_file,
        concatenate,
        signature_group
    )
    # Summary
    mut = 0
    elem = 0
    element_mutations_cutoff = False
    for k, v in mutations_d.items():
        mut += len(v)
        elem += 1
        if not element_mutations_cutoff:
            if len(v) >= element_mutations:
                element_mutations_cutoff = True
    logger.info('Validated elements in genomic regions: {}'.format(len(regions_d.keys())))
    logger.info('Validated elements with mutations: {}'.format(elem))
    logger.info('Total substitution mutations: {}'.format(mut))
    if not element_mutations_cutoff:
        raise excep.UserInputError('No element found with enough mutations to perform analysis')

    # Signature
    file_prefix = input_file.split('/')[-1].split('.')[0]
    path_pickle = os.path.join(path_cache, '{}_kmer_{}.pickle'.format(file_prefix, kmer))

    if not input_signature:
        """
        Calculate signatures
        
        By default, all substutions are taken into account to calculate the relative frequencies for each
        k-mer ref>alt. 
        
        Alternatively, when specified through '--signature-calculation region_normalized', k-mer mutation counts 
        can be normalized by the k-mer counts in the regions under analysis listed in 'regions_file'. In this case, 
        only substitutions that fall inside the regions will contribute to the signature calculation.         
        
        For both options, k-mers are not collapsed (192 channels) and do not include N (unknown reference nucleotides).
        """
        logger.info('Computing signature{}...'.format('s for each group' if signature_group else ''))

        if signature_calculation == 'region_normalized':
            normalize_regions_file, signature_calc_function = regions_file, bgsign.normalize
        else:
            normalize_regions_file, signature_calc_function = None, bgsign.relative_frequency

        signatures_dict = signature_calc_function(
            mutations_file=input_file,
            regions_file=normalize_regions_file,
            kmer_size=int(kmer),
            genome_build=genome,
            collapse=None,
            group=signature_group,
            cores=cores
        )
        # Reformat dictionary
        if not signature_group:
            signatures_to_pickle = {file_prefix: signatures_dict}
        else:
            signatures_to_pickle = signatures_dict
        # Save to cache
        with open(path_pickle, 'wb') as fd:
            pickle.dump(signatures_to_pickle, fd, protocol=2)
        logger.info('Signature{} computed'.format('s' if signature_group else ''))
    else:
        try:
            load_sign = bgsign.file.load(file=input_signature)
        except UnicodeDecodeError:
            raise excep.UserInputError('Error in input signatures file {}\n'
                                       'Please, check signatures file format (JSON)'.format(input_signature))
        # Check format and save pickle to cache
        keys = set(load_sign.keys())
        if not signature_group:
            # Expects 'file_prefix' to be a key in signatures dictionary of dictionaries
            if file_prefix in keys:
                with open(path_pickle, 'wb') as fd:
                    pickle.dump(load_sign, fd, protocol=2)
            # When dictionary has k-mer (ex. AAA>T) as keys, accepts dictionary and adds extra key for 'file_prefix'
            elif '>' in list(keys)[0]:
                with open(path_pickle, 'wb') as fd:
                    pickle.dump({file_prefix: load_sign}, fd, protocol=2)
            # Error, 'file_prefix' and k-mers not found as keys, check format
            else:
                raise excep.UserInputError('Incorrect format for input signature file {}'.format(input_signature))
        else:
            if not '>' in list(keys)[0]:
                with open(path_pickle, 'wb') as fd:
                    pickle.dump(load_sign, fd, protocol=2)
            # n signature dictionaries are expected (n = groups)
            else:
                raise excep.UserInputError('Groups are missing in signature dictionary at {}'.format(input_signature))

        logger.info('Input signature{} ready'.format('s' if signature_group else ''))

    # Initialize Experiment class variables and run
    elements_results, clusters_results, global_info_results = exp.Experiment(
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
                                                                            mutab_config,
                                                                            element_mutations,
                                                                            cluster_mutations,
                                                                            smooth_window,
                                                                            cluster_window,
                                                                            int(kmer),
                                                                            n_simulations,
                                                                            simulation_mode,
                                                                            simulation_window,
                                                                            cores,
                                                                            clustplot,
                                                                            seed
                                                                            ).run()

    # Write elements results (returns list of ranked elements)
    sorted_list_elements = postp.write_element_results(
        genome=genome,
        results=(elements_results, global_info_results),
        directory=output_directory,
        file=elements_output_file,
        is_gzip=gzip
    )
    logger.info('Elements results calculated')

    # Write clusters results
    postp.write_cluster_results(
        genome=genome,
        results=(clusters_results, global_info_results),
        directory=output_directory,
        file=clusters_output_file,
        sorter=sorted_list_elements,
        is_gzip=gzip
    )
    logger.info('Clusters results calculated')

    # Cluster plot
    if clustplotsig:
        info_cluster_plots = cplot.make_clustplot_sig(
            elements_results,
            clusters_results,
            global_info_results,
            directory=output_directory
        )
        for message in info_cluster_plots:
            logger.info(message)

    elif clustplot:
        info_cluster_plots = cplot.make_clustplot(
            elements_results,
            clusters_results,
            global_info_results,
            directory=output_directory
        )
        for message in info_cluster_plots:
            logger.info(message)

    # Quantile-quantile plot
    if qqplot:
        if len(elements) < 30:
            logger.warning('Q-Q plot will be generated for less than 30 elements')
        input_qqplot_file = os.path.join(output_directory, elements_output_file)
        output_qqplot_file = os.path.join(output_directory, 'quantile_quantile_plot.png')
        message = qplot.make_qqplot(
            file=input_qqplot_file,
            output=output_qqplot_file,
            col_values='P_ANALYTICAL',
            top=10
        )
        if message:
            logger.warning(message)
        else:
            logger.info('Q-Q plot plot generated at : {}'.format(output_qqplot_file))
    logger.info('Finished')


if __name__ == '__main__':
    main()
