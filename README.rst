.. _readme:

OncodriveCLUSTL
================

OncodriveCLUSTL is a sequence-based clustering method to identify significant clustering
signals in nucleotide sequence.

One of the main goals of cancer research is the identification of the genomic elements that drive tumorigenesis. OncodriveCLUSTL is a new nucleotide sequence-based clustering algorithm to detect significant clustering signals across genomic regions. OncodriveCLUSTL is based on a local background model derived from the nucleotide context mutational probabilities of the cohort under study. Our method is able to identify well-known cancer drivers in coding regions and it can be applied to non-coding regions and non-human data.

.. _readme license:

License
-------

OncodriveCLUSTL is available to the general public subject to certain conditions described in its `LICENSE <LICENSE>`_.


.. _readme install:

Installation
------------

OncodriveCLUSTL depends on Python 3.5 and some external libraries. We recommend to install it using the `Anaconda Python distribution <https://www.anaconda.com/download/>`_::

        $ conda install -c bbglab oncodriveclustl


OncodriveCLUSTL can also be installed using pip::

        $ pip install oncodriveclustl

You can obtain the latest code from the repository and install it with pip::

        $ git clone git@bitbucket.org:bbglab/oncodriveclustl.git
        $ cd oncodriveclustl
        $ conda install -c anaconda pycurl
        $ pip install .

.. note::

    The first time that you run OncodriveCLUSTL with a given reference genome, it will download it from our servers. By default the
    downloaded datasets go to ``~/.bgdata``. If you want to move these datasets to another folder you have to define the
    system environment variable BGDATA_LOCAL with an export command.

The following command will show you the help::

        $ oncodriveclustl --help


You can also download a singularity image
from the `downloads page <https://bitbucket.org/bbglab/oncodriveclustl/downloads/>`_.
Running the image invokes the oncodriveclustl command.

.. _readme inputdata:

Input data
---------------
OncodriveCLUSTL only requires two main inputs, the mutations file and the annotations file:

- Mutations file. TSV file containing SNVs (substitutions) mapped to a reference genome (e.g., human hg19 or mouse c3h). If other mutation types are present (insertions, deletions, double base substitutions, etc), they will be filtered out during the analysis. This file must contain, at least, the following 5 columns with header:

    1. CHROMOSOME: 1, 2,..., X, Y
    2. POSITION: Integer indicating the position of the mutation
    3. REF: Reference nucleotide
    4. ALT: Alternate nucleotide
    5. SAMPLE: Identifier of the sample

    Additional columns are:

    6. CANCER_TYPE: Type of tumor. When specified, OncodriveCLUSTL will calculate one mutational profile for each cancer type and mutations will be randomized accordingly.
    7. SIGNATURE: User-defined group to compute k-mer nucleotide mutational probabilities. When specified, OncodriveCLUSTL will calculate one mutational profile for each group and will randomize each mutation accordingly.

.. note::
    OncodriveCLUSTL assumes all SNVs are mapped to the positive strand.

.. warning::
    When using the '--signature-group' option, please check that the number of mutations per group is sufficient for an accurate signatures calculation.


- Annotations file. TSV file containing the coordinates of genomic elements (GEs). This file must contain, at least, the following 5 columns with header:

    1. CHROMOSOME: 1, 2,..., X, Y
    2. START: Starting position of the genomic region
    3. END: Final position of the genomic region
    4. ELEMENT: Identifier of the GE
    5. SYMBOL: Symbol of the GE. OncodriveCLUSTL will analyze GEs as SYMBOL + ELEMENT.

    Additional columns are:

    6. STRAND: Strand of the GE coordinates ("+" or "-").

.. warning::
    Coordinates of a given GE cannot overlap.

You can check the input formats in the files provided in the example.

If you have a VCF file or directory of VCF files containing somatic mutations, you can run our VCF parser to obtain a tabular file compatible with OncodriveCLUSTL input format::

       $ parse_vcf -i [INPUT_DIRECTORY] -o [OUTPUT_FILE]

Please, check 'parsers/vcf.py' module for more details.

If you would like to run OncodriveCLUSTL using a per-calculated signature or mutational profile, you need to provide a dictionary containing the reference k-mer to alternate mutational probabilities in JSON format::

        {
            "my_dataset": {
                "GCA>G": 0.02424271083094251,
                "AGC>A": 0.023005887103025254,
                "ACG>T": 0.037613802858829135,
                "CGA>C": 0.10691031051670515,
                "GAC>G": 0.017846071811001615,
                "TTC>A": 0.024003748061871697,
                "CTT>G": 0.024149863672267024,
                "GGA>T": 0.011178562948734577,
                "AGG>C": 0.010654720767868876,
                "GGG>C": 0.012031686292218055,
                "CAA>T": 0.014478959792844522,
                "TGA>A": 0.01255651801972085,
                "GGA>A": 0.011178562948734577,
                "CGA>A": 0.03563677017223505,
                "TCC>T": 0.011158347971568658,
                "GCC>A": 0.010952316565906438,
                ...
            }
        }

OncodriveCLUSTL requires non-collapsed k-mer probabilities (192 for tri-nucleotides, 3072 for penta-nucleotides).

.. _readme outputdata:

Output data
---------------
OncodriveCLUSTL generates three output files:

- Elements results file ('elements_results.txt'). TSV file containing results of the analyzed elements:

    1. SYMBOL: GE symbol
    #. ENSID: GE ID
    #. CGC: True if GE in the COSMIC Cancer Gene Census (CGC) list (Sondka et al., 2018)
    #. CHROMOSOME: 1, 2,..., X, Y
    #. STRAND: Strand of the GE ("+" or "-")
    #. LENGTH: length (bp) of the GE
    #. TOTAL_MUT: total substitutions observed in the GE
    #. CLUSTERED_MUT: number of substitutions in a cluster
    #. CLUSTERS: number of clusters
    #. SIM_CLUSTERS: number of simulated clusters
    #. SCORE: GE score
    #. P_EMPIRICAL: empirical p-value of the GE
    #. Q_EMPIRICAL: empirical q-value of the GE
    #. P_ANALYTICAL: analytical p-value of the GE
    #. Q_ANALYTICAL: analytical q-value of the GE
    #. P_TOPCLUSTER: analytical p-value of the cluster with highest cluster score
    #. Q_TOPCLUSTER: analytical q-value of the cluster with highest cluster score


- Clusters results file ('clusters_results.tsv'). TSV file containing results of the clusters observed in the analyzed elements:

    1. RANK: Position of the GE in the list of
    #. SYMBOL: GE symbol
    #. ENSID: GE ID
    #. CGC: True if GE in the CGC list
    #. CHROMOSOME: 1, 2,..., X, Y
    #. STRAND: Strand of the GE ("+" or "-")
    #. COORDINATES: genomic coordinates of the cluster. It can be 'coord1,coord2' for clusters inside a single region or 'coord1,coord2;coord3,coord4' for those spanning regions (--concatenate flag)
    #. MAX_COORD: genomic position with the highest smoothing score inside the cluster
    #. WIDTH: cluster's width (pb)
    #. N_MUT: number of substitutions in the cluster
    #. N_SAMPLES: number of samples with a mutation in the cluster
    #. FRA_UNIQ_SAMPLES: proportion of unique samples mutated in the cluster out of the total of mutations in the cluster
    #. SCORE: cluster score
    #. P: analytical p-value of the cluster

- Log file ('results.log'). TXT file containing OncodriveCLUSTL's run information

.. _readme usage:

Usage
---------------
OncodriveCLUSTL is meant to be used through the command line.

- '-i', '--input-file': File containing mutations (required)
- '-r', '--regions-file': GZIP compressed file with the genomic regions to analyze (required)
- '-o', '--output-directory': Output directory to be created (required)
- '-sig', '--input-signature': File containing input context based mutational probabilities
- '-ef', '--elements-file': File with the symbol of a set GEs to analyze, one per row
- '-e', '--elements': Symbol of the GE(s) to analyze
- '-g', '--genome': Genome to use. Default is hg19
- '-emut', '--element-mutations': Cutoff of element mutations. Default is 2
- '-cmut', '--cluster-mutations': Cutoff of cluster mutations. Default is 2
- '-sw', '--smooth-window': Smoothing window. Default is 11
- '-cw', '--cluster-window': Cluster window. Default is 11
- '-kmer', '--kmer': K-mer nucleotide context (3 or 5) to calculate mutational probabilities. Default is 3
- '-n', '--n-simulations': Number of simulations. Default is 1000
- '-sim', '--simulation-mode': Simulation mode. Default is 'mutation_centered'
- '-simw', '--simulation-window': Simulation window. Default is 31
- '-sigcalc', '--signature-calculation': calculation of mutational probabilities as mutation frequencies ('frequencies') or k-mer mutation counts normalized by k-mer region counts ('region_normalized'). Default is frequencies
- '-siggroup', '--signature-group': Header of the column to group signatures calculation ('SIGNATURE', 'SAMPLE', 'CANCER_TYPE'). One mutational profile will be calculated for each group.
- '-c', '--cores': Number of cores to use in the computation. By default it uses all the available cores
- '--seed': seed to use in the simulations
- '--log-level': Verbosity of the logger. Default is 'info'
- '--concatenate': Calculate clustering on concatenated genomic regions (e.g., exons in coding sequences)
- '--clustplot': Needle plot with clusters for a GE
- '--qqplot': Quantile-quantile (Q-Q) plot for a dataset
- '--gzip': Gzip compress files

.. note::
    When using simulation mode 'mutation_centered', simulation windows can be simulated outside the GE.

.. note::
    When using '--signature-calculation region_normalized', k-mer mutation counts will be normalized by k-mer nucleotide counts in the genomic regions provided as input ('--regions-file').

.. _readme example:

Run the example
---------------

After installing OncodriveCLUSTL, you can run an example of TCGA pancreatic adenocarcinomas (Ellrott et al. 2018) for coding regions (Mularoni et al., 2016) using 1000 simulations.
First you need to download the example folder. Then you run OncodriveCLUSTL with default mode and parameters as::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example

The results will be saved in a folder named ``output_example``.

You can compute a more sophisticated analysis using non-default parameters and generate a quantile-quantile plot by typing::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example -sw 15 -cw 15 -simw 35 -sim region_restricted --concatenate --qqplot

If you want to run a specific GE and generate a plot its observed clusters, you can type::

        $ oncodriveclustl -i ~/example/PAAD.tsv.gz -r ~/example/cds.hg19.regions.gz -o ~/example/output_example -sw 15 -cw 15 -simw 35 -sim region_restricted --concatenate --clustplot -e KRAS


