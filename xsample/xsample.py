import cirro
import pandas as pd
import scipy.stats as stats

import anndata as ad
import os
import glob

def ligrec_pair_from_list(ligrec_list, pair=None, type='means', axis=1,
                          samples=None, plot=False):
    if type not in ['means', 'pvalues']:
        raise ValueError("type must be either 'means' or 'pvalues'")
    ligrecs = [ligrec[type] for ligrec in ligrec_list]
    comb = pd.concat(ligrecs, axis=axis, keys=samples, join='outer')

    comb.columns.names = ['sample', 'cluster_1', 'cluster_2']

    if pair is not None:
        idx = pd.IndexSlice
        res = comb.loc[:,idx[:,pair[0],pair[1]]].copy()
    else:
        res = comb.copy()

    if(plot):
        plot_hist(res, title="LigRec {} for clusters {}".format(type, pair))

    return res

def plot_hist(df_pair, title=None, save=True):
    import matplotlib.pyplot as plt
    import seaborn as sns
    x = 8
    y = df_pair.shape[0] / 4

    plt.figure(figsize=(x, y))
    sns.heatmap(df_pair, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title(title)
    if save:
        plt.savefig(title.replace(" ", "_")+".png", dpi=300, bbox_inches='tight')

def xsample_ttest(df, group1, group2):
    res = df.copy()
    test = stats.ttest_ind(res[group1], res[group2], axis=1)
    res['pval'] = test.pvalue
    res['statistic'] = test.statistic
    res.dropna(inplace=True)
    res['pval_adj'] = stats.false_discovery_control(res['pval'], method='bh')
    res.sort_values('pval_adj', inplace=True)

    return res[sorted(res.columns)]