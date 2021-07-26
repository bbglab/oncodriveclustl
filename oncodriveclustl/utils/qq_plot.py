"""
Contains functions to generate a quantile-quantile plot for a dataset
"""
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def preprocess_dataframe(file, col_values, p_values_cutoff=None, remove_cgc=False):
    """
    Preprocess dataframe for subsequent analysis.

     Args:
        file (str): path to file with results
        col_values (str): p-value column name in dataframe
        p_values_cutoff (float): select p-values greater than a cutoff
        remove_cgc (bool): True removes CGC genes from the dataframe

    Returns:
        df (DataFrame): dataframe with results

    """

    q_values = 'Q_EMPIRICAL' if col_values == 'P_EMPIRICAL' else 'Q_ANALYTICAL'
    df = pd.read_csv(file, sep='\t', header=0)

    # Preprocess
    df = df[np.isfinite(pd.to_numeric(df[col_values]))].copy()
    df.sort_values(by=[col_values, 'SCORE', 'CGC'], ascending=[True, False, False], inplace=True)

    if remove_cgc:
        df = df.loc[df['CGC'] == False]

    if p_values_cutoff:
        df = df.loc[df[col_values] > p_values_cutoff]

    if len(df[col_values].tolist()) > 0:
        min_pvalue = min([x for x in df[col_values].tolist() if x > 0])
        obs = list(df[col_values].map(lambda x: -np.log10(x) if x > 0 else -np.log10(min_pvalue)))
        data = pd.DataFrame({
            'SYMBOL': df['SYMBOL'],
            'SCORE': df['SCORE'],
            'CGC': df['CGC'],
            'obs': obs,
            'Q': df[q_values]
        })

        data.sort_values(by=['obs', 'SCORE', 'CGC'], ascending=[True, True, False], inplace=True)
        exp = -np.log10(np.arange(1, len(data) + 1) / len(data))
        exp.sort()
        data['exp'] = exp

        df = data

        return df

    else:
        return 'Error'


def make_qqplot(file, col_values, top=10, output=None):
    """
    Generate quantile-quantile plot (QQ-plot) from an OncodriveCLUSTLs output

     Args:
        file (str): path to file with elements results (e.g., 'elements_results.txt')
        col_values (str): p-value column name in dataframe
        top (int): number or top ranking elements to plot
        output (str): path to output file

    Returns:
        error (str):  returns error if q-q plot couldn't be computed

    """
    # Plot params
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.set_xlabel(r"Expected p-value $(-\log_{10})$", fontsize=14)
    ax.set_ylabel(r"Observed p-value $(-\log_{10})$", fontsize=14)
    ax.set_axisbelow(True)
    ax.grid(color='lightgrey', linestyle='--', linewidth=1)

    data = preprocess_dataframe(file=file, col_values=col_values)

    if type(data) != str:

        # Scatter of the p-values
        color = np.where(data['Q'] < 0.01, '#ca0020', '#0571b0')
        ax.scatter(data['exp'], data['obs'], color=color, alpha=0.5)

        # Null hypothesis
        ax.plot(data['exp'], data['exp'], color='#ca0020', linestyle='--')

        # Plot significant elements
        if len(data) > top:
            x_text = ax.get_xlim()[1] * 1.05
            y_text = max(data['obs'])
            delta = y_text / 15
            sign = data[-top:].copy()
            sign.sort_values(by='exp', ascending=False, inplace=True)
            for count, line in sign.iterrows():
                color = 'black' if line['CGC'] is True else 'black'
                weight = 'demi' if line['CGC'] is True else 'normal'
                ax.annotate(line['SYMBOL'], xy=(line['exp'], line['obs']),
                            xytext=(x_text, y_text), color=color, fontsize=12, fontweight=weight,
                            arrowprops=dict(color='grey', shrink=0.05, width=0.5, headwidth=5, alpha=0.20), )
                y_text -= delta

        if output:
            plt.savefig(output, bbox_inches='tight')
        plt.close()

        return None

    return "Q-Q plot couldn't be generated. There are no elements with p-values"
