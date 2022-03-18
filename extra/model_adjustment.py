"""
Contains functions to adjust OncodriveCLUSTL models
"""
import os

import numpy as np
import pandas as pd
from scipy import stats


"""
Functions contained in this module can be used to adjust OncodriveCLUSTL's models based on a set of tested parameters, 
including: smoothing, clustering, sampling/simulation and kmer-nucleotides. Given a set of 'elements_results.txt' 
output files for a dataset, ran using different combinations of parameters, the steps are as follows: 

1) Calculate the KS statistic and the logarithmic weighted enrichment for each combination of 
parameters (configuration):  
1.1) 'calculate_ks': computes the KS statistic
1.2) 'evaluate_enrichment_method': computes the CGC enrichment. It requires a list of CGC genes. 

2) Results from step 1 (e.g., DataFrame format) are evaluated using 'calculate_best_config' and the best configuration
 is retrieved (details in Supplementary Methods). 

"""


def calculate_ks(file, pvalue, random_set=1000, alpha=0.1):
    """
    Calculate Kolmogorov-Smirnov statistic (KS) for raw-pvalues fitness to the uniform distribution (expected).
    KS is returned as positive if the difference between observed and expected is greater than 0 in subset of
    p-values defined by alpha cutoff.

    Args:
        file (str): path to file with elements results
        pvalue (str): p-value column name in dataframe to calculate KS
        random_set (int): random set of genes to calculate KS
        alpha (float): cutoff to calculate KS sign

    Returns:
        ks_statistic (float): KS statistic
        (float): ks_pvalue
        (int): number of p-values used to calculate KS
    """

    df = pd.read_csv(file, sep='\t', header=0)
    df = df[np.isfinite(pd.to_numeric(df[pvalue]))].copy()
    df.sort_values(by=[pvalue, 'SCORE', 'CGC'], ascending=[True, False, False], inplace=True)

    if len(df) != 0:
        significant_pvalues = df.loc[df[pvalue] > alpha][pvalue]

        # Subset p-values > 0.1
        if len(significant_pvalues) > 1000:
            pvalues = np.random.choice(significant_pvalues, size=random_set, replace=False)
        else:
            pvalues = significant_pvalues

        # Calculate KS
        ks = stats.kstest(pvalues, 'uniform')

        # Calculate KS sign
        observed = len(df.loc[df[pvalue] >= alpha])
        expected = (1 - alpha) * len(df)
        if observed < expected:
            ks_statistic = ks.statistic
        else:
            ks_statistic = - ks.statistic

        return ks_statistic, ks.pvalue, len(pvalues)
    else:
        return float('nan'), float('nan'), float('nan')


def get_weight(i, weight):
    """
    Weight contributions according to ranking

    Args:
        i (int): the ith ranking
        weight (str): the type of weighting [log, normal]

    Returns:
        (float): the weight of the ith position

    """
    if weight == "log":
        return 1.0 / np.log2(i + 2)
    if weight == "normal":
        return 1.0 / i


def calculate_percentage_cgc(ranking, cgc_genes):
    """
    Calculate the percentage of CGC genes in the input list

    Args:
        ranking (list): the input list of the ranked genes
        cgc_genes (set): set of genes included in the Cancer Gene Census (CGC)

    Returns:
        (float): percentage of cgc genes in the list

    """

    n = float(sum([1.0 if gene in cgc_genes else 0.0 for gene in ranking]))
    return n / len(ranking)


def evaluate_enrichment_method(file, cgc_genes, pvalue, weight="log", ranking_limit=40):
    """

    Args:
        file (str): path to file with elements results
        cgc_genes (set): set of genes in the CGC
        pvalue (str): p-value column name in dataframe to calculate KS
        weight (str): normalization of the weight. [log,normal] default: log
        ranking_limit (int): limit to calculate the area under the curve. Default: 40

    Returns:
        (foat): The weighted area under the curve of the CGC enrichment

    """

    df = pd.read_csv(file, sep='\t', header=0)
    df = df[np.isfinite(pd.to_numeric(df[pvalue]))].copy()
    df.sort_values(by=[pvalue, 'SCORE', 'CGC'], ascending=[True, False, False], inplace=True)

    ranking = df['SYMBOL'].tolist()

    if len(ranking) > ranking_limit:
        ranking = ranking[0:ranking_limit]
    xticks = range(len(ranking))
    area = 0.0
    for i in xticks:
        weight_i = get_weight(i, weight)
        x_i = calculate_percentage_cgc(ranking[0:i + 1], cgc_genes)
        area += x_i * weight_i

    return area


def calculate_best_config(file):
    """
    Calculate the best OncodriveCLUSTL configuration of parameters for a dataset
    Args:
        file (path): path to file containing KS statistics and CGC enrichments of a dataset using a set of
                        configurations
    Returns:
        (tuple): KS statistic, enrichment, simulation window, smoothing window, clustering window, kmer-nucleotide

    """

    df = pd.read_csv(file, sep='\t', header=0)
    df_nonempty = df[np.isfinite(df['KS'])].copy()

    if len(df_nonempty) != 0:
        # Sort by KS and get best configurations
        df_nonempty.sort_values(by=['KS'], ascending=[True], inplace=True)
        df_nonempty = df_nonempty.reset_index(drop=True)
        candidate_configs = df_nonempty.loc[df_nonempty['KS'] <= min(df_nonempty['KS'].tolist()) + 0.1].copy()

        # Sort by CGC enrichment and get top 1
        candidate_configs.sort_values(by=['ENRICH_LOG'], ascending=[False], inplace=True)
        candidate_configs = candidate_configs.reset_index(drop=True)
        candidate_config = candidate_configs.iloc[0]

        # Parameters
        ks = str(candidate_config['KS'])
        enrich = str(candidate_config['ENRICH_LOG'])
        simw = str(candidate_config['SIMULATION'])
        smo = str(candidate_config['SMOOTH'])
        clu = str(candidate_config['CLUSTER'])
        kmer = 'kmer' + str(candidate_config['KMER'])
        results = (ks, enrich, simw, smo, clu, kmer)
    else:
        results = ()

    return results
