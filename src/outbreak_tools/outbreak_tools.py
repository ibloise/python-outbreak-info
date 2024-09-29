import pandas as pd
import numpy as np
import frozendict
import requests
import gzip
import yaml
import json

from src.outbreak_tools import outbreak_clustering

def get_colors(lins, brighten, lineage_key):
    """Heuristically assign colors to lineages to convey divergence.

     :param lins: list of (group root) lineage names.
     :param brighten: boolean allowing to brighten some lineages.
     :param lineage_key: dict mapping lineage names to tree nodes.

     :return: a list of lineage colors represented as hsv tuples."""
    colors = np.searchsorted(
        sorted([lin['alias'] for lin in lineage_key.values()]),
        [lineage_key[lin]['alias'] for lin in lins] )
    colors = colors ** 2
    colors = (colors - np.min(colors)) / (np.max(colors)-np.min(colors)) * 0.75
    return [(color, 1, 0.55 + 0.25*b) for color, b in zip(colors, brighten)]

def get_riverplot_baseline(prevalences, loads, k=128):
    """Find a baseline for drawing a river plot (a shifted scaled stacked area plot) that minimizes visual shear.

     :param prevalences: pandas df of lineage prevalences over time (See lineage_cl_prevalence())
     :param loads: pandas series of viral loads or other scaling data.
     :param k: number of iterations to run.

     :return: a pandas series representing the vertical offset of the bottom edge of the river plot."""
    c = prevalences.mul(loads.interpolate(), axis=0).dropna()
    d = c.div(loads.dropna(), axis=0)
    shear = lambda O: (c.cumsum(axis=1).add(O, axis=0).rolling(window=2).apply(np.diff).mul(d)**2).sum().sum()
    Ot = -loads.dropna()/2
    for n in range(k):
        O = np.random.normal(size=Ot.shape) / (n+48) * 2
        if shear(O+Ot) < shear(Ot):
            Ot += O
            Ot -= np.mean(Ot)
    return pd.Series(Ot, c.index).reindex(prevalences.index).interpolate()

def first_date(samples, by='collection_site_id'):
    """Get the earliest date among samples for each unique value in some column.

     :param samples: pandas dataframe of samples indexed by date.
     :param by: name of target column.

     :return: a pandas series mapping unique values to dates"""
    return samples.reset_index(level=0, names='date').groupby(by)['date'].min()

def get_ww_weights(df, loaded=True):
    """Get default weights for aggregating wastewater data.

     :param df: pandas dataframe of samples to be weighted.
     :param loaded: whether to incorporate viral load data.

     :return: a pandas series of sample weights."""
    weights = df['ww_population'].fillna(1000)
    if loaded: weights *= df['normed_viral_load'].fillna(0.5)
    return weights

def const_idx(df, const, level):
    """Set one level of a multi-indexed df to a constant.

     :param df: multi-indexed pandas dataframe.
     :param const: constant value to assign to index.
     :param level: level of index to change.

     :return: the modified dataframe."""
    df = df.copy()
    df.index = df.index.set_levels([const]*len(df), level=level, verify_integrity=False)
    return df

def datebin_and_agg(df, weights=None, freq='7D', rolling=1, startdate=None, enddate=None, column='prevalence', norm=True, variance=False, log=False, trustna=1):
    """Gather and aggregate samples into signals.

     :param df: A multi-indexed pandas dataframe; df.index[0] is assumed to be a date and df.index[1] a categorical.
     :param weights: A pandas series of sample weights. `None` is appropriate for clinical df[column] and `get_ww_weights` for wastewater.
     :param freq: Length of date bins as a string.
     :param rolling: How to smooth the data; an int will be treated as a number of bins to take the rolling mean over, and an array as a kernel.
     :param startdate: Start of date bin range as YYYY-MM-DD string.
     :param enddate: End of date bin range as YYYY-MM-DD string.
     :param column: Data column to aggregate.
     :param norm: Whether to normalize so that aggregated values across all categories in a date bin sum to 1.
     :param variance: Whether to return the rolling variances along with the aggregated values.
     :param log: Whether to do the aggregation in log space (geometric vs arithmetic mean).
     :param trustna: How much weight to place on the nan=0 assumption.

     :return: A pandas dataframe of aggregated values with rows corresponding to date bins and columns corresponding to categories."""
    if startdate is None: startdate = df.index.get_level_values(0).min()
    if enddate is None: enddate = df.index.get_level_values(0).max()
    startdate = pd.to_datetime(startdate)-pd.Timedelta('1 day')
    enddate = pd.to_datetime(enddate)+pd.Timedelta('1 day')
    if freq is None: dbins = [pd.Interval(startdate, enddate)]
    else: dbins = pd.interval_range(startdate, enddate, freq=freq)
    bins = pd.IntervalIndex(pd.cut(pd.to_datetime(df.index.get_level_values(0)) + pd.Timedelta('1 hour'), dbins))
    if weights is None: weights = df.apply(lambda x: 1, axis=1)
    df, weights, bins = df[~bins.isna()], weights[~bins.isna()], bins[~bins.isna()]
    eps = 1e-8
    clog, cexp = [(lambda x:x, lambda x:x), (lambda x: np.log(x+eps), lambda x: np.exp(x))][int(log)]
    if isinstance(rolling, int): rolling = [1] * rolling
    else: rolling = np.array(list(rolling))
    rolling = rolling / np.sum(rolling)
    rollingf = lambda x: np.convolve(rolling,  np.pad(x.fillna(0), len(rolling)//2, 'edge'), mode='valid')
    bindex = pd.MultiIndex.from_arrays([bins, df.index.get_level_values(1).str.split('-like').str[0].str.split('(').str[0]])
    def binsum(x):
        x = x.to_frame().groupby(bindex).sum(min_count=1)
        x = x.set_index(pd.MultiIndex.from_tuples(x.index)).unstack(1)
        x.columns = x.columns.droplevel(0)
        return x.reindex(dbins).sort_index().apply(rollingf, axis=0)
    nanmask = (~np.isnan(df[column])).astype(int)
    nanmask = np.clip(nanmask + trustna, 0, 1)
    prevalences = binsum(weights*nanmask*clog(df[column].fillna(0)))
    if norm:
        prevalences = prevalences.apply(cexp)
        denoms = prevalences.sum(axis=1)
        prevalences = prevalences.div(denoms, axis=0)
    else:
        denoms = binsum(weights*nanmask)
        prevalences = prevalences.div(denoms)
        prevalences = prevalences.apply(cexp)
        prevalences = prevalences.where(binsum(~np.isnan(df[column])) > 0, np.nan)
    if variance:
        means = np.array(prevalences)[
            prevalences.index.get_indexer_for(bins),
            prevalences.columns.get_indexer_for(df.index.get_level_values(1))]
        variances = binsum((weights*nanmask*(clog(df[column].fillna(0)) - clog(means)))**2)
        variances = variances.div(denoms**2, **({'axis': 0} if norm else {}))
        if log: variances = variances * prevalences**2
    return (prevalences, variances) if variance else prevalences

def get_tree(url='https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/curated_reports_prep/lineages.yml'):
    """Download and parse the lineage tree (derived from the Pangolin project).

     :param url: The URL of an outbreak-info lineages.yml file.

     :return: A nested tree of frozendicts representing the phylogenetic tree."""
    response = requests.get(url)
    response = yaml.safe_load(response.content.decode("utf-8"))
    lin_names = sorted(['*'] + [lin['name'] for lin in response])
    lindex = {lin:i for i,lin in enumerate(lin_names)}
    lineage_key = dict([(lin['name'], lin) for lin in response if 'parent' in lin])
    def get_children(node, lindex):
        return tuple( frozendict.frozendict({ 'name': lineage_key[c]['name'], 'lindex': lindex[lineage_key[c]['name']],
                                              'alias': lineage_key[c]['alias'], 'parent': node['name'],
                                              'children': get_children(lineage_key[c], lindex) })
                         for c in node['children'] if c in lineage_key and lineage_key[c]['parent'] == node['name'] )
    roots = tuple( frozendict.frozendict({ 'name': lin['name'], 'lindex': lindex[lin['name']],
                                           'alias': lin['alias'], 'parent': '*', 'children': get_children(lin, lindex) 
                             }) for lin in response if not 'parent' in lin )
    return frozendict.frozendict({ 'name': '*', 'lindex': lindex['*'], 'alias': '*',
                                   'parent': '*', 'children': roots })

def write_compressed_tree(tree, file='./tree.json.gz'):
    with gzip.open(file, 'wb') as f:
        f.write(json.dumps(tree).encode('utf-8'))
def read_compressed_tree(file='./tree.json.gz'):
    with gzip.open(file, 'rb') as f:
        return frozendict.deepfreeze(json.loads(f.read()))

def cluster_df(df, clusters, tree, lineage_key=None, norm=False, include_K=False):
    """Aggregate the columns of a dataframe into some phylogenetic groups.

     :param df: A dataframe of prevalence signals. Rows are assumed to be date bins and columns are assumed to be lineages.
     :param clusters: A tuple (U,V) of sets of root nodes representing clusters (from cluster_lineages).
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :param norm: Whether to assume that values in a row should sum to one.
     :param include_K: Whether to include fixed lineages in the output.

     :return: A tuple (data,names,is_inclusive) where data is the input dataframe with aggregated and relabeled columns, names contains the names of the root lineages for each column/group, and is_inclusive indicates whether the column's root is in U or V."""
    if lineage_key is None: tree = get_lineage_key(tree)
    (U,V,K) = clusters
    if include_K:
        U = U|K
        K = set([])
    prevalences_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'      {u["alias"]}*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
    vlabels = [f'other {v["alias"]}*' + (f' ({v["name"]})' if v["name"] != v["alias"] else '') for v in V]
    legend = list(np.array(ulabels+vlabels)[order])
    clustered_prevalences = pd.DataFrame(
        { d: { label:outbreak_clustering.get_agg_prevalence(lin, a, U|V|K)
            for label, lin in zip(legend, lins) }
        for d,a in zip(dates, prevalences_dated) } ).transpose()
    if norm:
        clustered_prevalences[np.sum(clustered_prevalences, axis=1) < 0.5] = pd.NA
        clustered_prevalences['other **'] += 1 - clustered_prevalences.sum(axis=1)
        clustered_prevalences['other **'] = np.clip(clustered_prevalences['other **'], 0, 1)
    return clustered_prevalences, [lin['name'] for lin in lins], np.array([1]*len(U)+[0]*len(V))[order]