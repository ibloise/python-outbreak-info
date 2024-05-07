import pandas as pd
import re
import warnings
import numpy as np
import frozendict
from collections import OrderedDict
import yaml
import requests
import json
import gzip
    
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
     :param prevalences: pandas df of lineage prevalences over time.
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

def datebin_and_agg(df, weights=None, freq='7D', startdate=None, enddate=None, column='prevalence', norm=True):
    """Gather and aggregate samples into signals.
     :param df: A multi-indexed pandas dataframe; df.index[0] is assumed to be a date and df.index[1] a categorical.
     :param weights: A pandas series of sample weights. `None` is appropriate for clinical data and `get_ww_weights` for wastewater.
     :param freq: Length of date bins as a string.
     :param startdate: Start of date bin range as YYYY-MM-DD string.
     :param enddate: End of date bin range as YYYY-MM-DD string.
     :param column: Data column to aggregate.
     :param norm: Whether to normalize so that aggregated values across all categories in a date bin sum to 1.
     :return: A pandas dataframe of aggregated values with rows corresponding to date bins and columns corresponding to categories."""
    if startdate is None: startdate = df.index.get_level_values(0).min()
    if enddate is None: enddate = df.index.get_level_values(0).max()
    startdate = pd.to_datetime(startdate)-pd.Timedelta('1 day')
    enddate = pd.to_datetime(enddate)+pd.Timedelta('1 day')
    if freq is None: dbins = [pd.Interval(startdate, enddate)]
    else: dbins = pd.interval_range(startdate, enddate, freq=freq)
    bins = pd.IntervalIndex(pd.cut(pd.to_datetime(df.index.get_level_values(0)), dbins))
    if weights is None: weights = df.apply(lambda x: 1, axis=1)
    df, weights, bins = df[~bins.isna()], weights[~bins.isna()], bins[~bins.isna()]
    def binsum(data, bins):
        data = data.to_frame().set_index(bins)
        return data.groupby(data.index).sum()
    agged_prevalences = binsum(df[column]*weights, pd.MultiIndex.from_arrays([bins, df.index.get_level_values(1)]))
    agged_prevalences = agged_prevalences.set_index(pd.MultiIndex.from_tuples(agged_prevalences.index)).unstack(1)
    agged_prevalences.columns = agged_prevalences.columns.droplevel(0)
    denoms = binsum((df[column] if norm else 1) * weights, bins)
    denoms = denoms[denoms.columns[0]]
    denoms.index = agged_prevalences.index
    agged_prevalences = agged_prevalences.div(denoms, axis=0)
    agged_prevalences = agged_prevalences.rename(columns = { c:c.split('-like')[0] for c in agged_prevalences.columns })
    return agged_prevalences.T.groupby(agged_prevalences.columns).sum().T
    
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

def get_compressed_tree(url='https://raw.githubusercontent.com/outbreak-info/python-outbreak-info/wastewater_sprint/tree.json.gz'):
    """Download the pre-parsed lineage tree (derived from the Pangolin project).
     :param url: The URL of the json tree file.
     :return: A nested tree of frozendicts representing the phylogenetic tree."""
    response = requests.get(url)
    return frozendict.deepfreeze(json.loads(gzip.decompress(response.content)))

def get_lineage_key(tree, field='name'):
    """Create a mapping of names to tree nodes.
     :param tree: The root of the phylo tree object (a frozendict).
     :param name: The field to map on.
     :return: An OrderedDict containing the map."""
    def get_names(tree):
        return np.concatenate([[(tree[field], tree)]] + [get_names(c) for c in tree['children']])
    return OrderedDict(sorted(get_names(tree), key=lambda x: x[0]))

def cluster_lineages(prevalences, tree, lineage_key=None, n=16, alpha=0.1):
    """Cluster lineages via greedy group-splitting on the phylo tree starting from the root based on some heuristics.
     :param prevalences: A dict, pandas series or other map between lineage names and (un-normalized) prevalences.
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :param n: The target number of clusters.
     :param alpha. Heuristic control in range (0, 1); higher values avoid more low-quality groups, but can prevent convergence on some data.
     :return: A tuple (U,V) of sets of group root lineages. Groups in U contain all descendant lineages of their roots, while groups in V are exclusive of some more distal groups in U or V."""
    if lineage_key is None: lineage_key = get_lineage_key(tree)
    prevalences = OrderedDict([(k,0) for k in lineage_key.keys()]) | dict(prevalences)
    prevalences = np.array(list(prevalences.values()))
    agg_prevalences = np.zeros_like(prevalences) - 1
    def init_agg_prevalences(node):
        if agg_prevalences[node['lindex']] < 0:
            agg_prevalences[node['lindex']] = prevalences[node['lindex']]
            agg_prevalences[node['lindex']] += np.sum([init_agg_prevalences(c) for c in node['children']])
        return agg_prevalences[node['lindex']]
    init_agg_prevalences(tree)
    def update_ancestors(node, diff, W):
        if not node in W and node['parent'] != node['name']:
            parent = lineage_key[node['parent']]
            agg_prevalences[parent['lindex']] += diff
            update_ancestors(parent, diff, W)
    def contains_descendant(node, nodeset):
        return node in nodeset or \
               len(set(node['children']) & nodeset) > 0 or \
               np.sum([contains_descendant(c, nodeset) for c in node['children']]) > 0
    U,V = set([tree]), set([])
    while len(U|V) < n:
        split_node = list(U|V)[np.argmax([ agg_prevalences[w['lindex']] * 0.1 + np.max(
            [0] + [agg_prevalences[c['lindex']] for c in w['children'] if not c in U|V] ) for w in list(U|V) ])]
        C = [c for c in split_node['children'] if not c in U|V]
        add_node = C[np.argmax([agg_prevalences[c['lindex']] for c in C])]
        update_ancestors(add_node, -agg_prevalences[add_node['lindex']], U|V)
        if split_node in U:
            U = U - set([split_node])
            V = V | set([split_node])
        if contains_descendant(add_node, U|V):
                V = V | set([add_node])
        else:   U = U | set([add_node])
        if len(U) > 1 and len(list(V - set([tree]))) > 1:
            drop_node = list(V - set([tree]))[np.argmin([agg_prevalences[v['lindex']] for v in list(V - set([tree]))])]
            if agg_prevalences[drop_node['lindex']] < alpha * np.mean([agg_prevalences[u['lindex']] for u in U]):
                V = V - set([drop_node])
                update_ancestors(drop_node, agg_prevalences[drop_node['lindex']], U|V)
    return U,V

def get_agg_prevalence(root, prevalences, W=set([])):
    """Compute the total prevalence of all lineages in a group.
     :param root: A tree node representing the root of the group (a frozendict).
     :param prevalences: A map from lineage names to prevalences.
     :param W: A set of tree nodes representing the roots of other groups. All nodes descended from root that are in W or are children of those nodes are excluded.
     :return: The prevalence of the group."""
    cs = [get_agg_prevalence(c, prevalences, W) for c in root['children'] if not c in W]
    return np.clip((prevalences[root['name']] if root['name'] in prevalences else 0) + np.sum(cs), 0, None)

def cluster_df(df, clusters, tree, lineage_key=None):
    """Aggregate the columns of a dataframe into some phylogenetic groups.
     :param df: A dataframe of prevalence signals. Rows are assumed to be date bins and columns are assumed to be lineages.
     :param clusters: A tuple (U,V) of sets of root nodes representing clusters (from cluster_lineages).
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :return: A tuple (data,names,is_inclusive) where data is the input dataframe with aggregated and relabeled columns, names contains the names of the root lineages for each column/group, and is_inclusive indicates whether the column's root is in U or V."""
    if lineage_key is None: tree = get_lineage_key(tree)
    (U,V) = clusters
    viral_load = None
    prevalences_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'      {u["alias"]}*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
    vlabels = [f'other {v["alias"]}*' + (f' ({v["name"]})' if v["name"] != v["alias"] else '') for v in V]
    legend = list(np.array(ulabels+vlabels)[order])
    clustered_prevalences = pd.DataFrame(
        { d: { label:get_agg_prevalence(lin, a, U|V)
            for label, lin in zip(legend, lins) }
        for d,a in zip(dates, prevalences_dated) } ).transpose()
    clustered_prevalences[np.sum(clustered_prevalences, axis=1) < 0.5] = pd.NA
    clustered_prevalences['other **'] += 1 - clustered_prevalences.sum(axis=1)
    clustered_prevalences['other **'] = np.clip(clustered_prevalences['other **'], 0, 1)
    if viral_load is not None: clustered_prevalences = clustered_prevalences.join(viral_load)
    return clustered_prevalences, [lin['name'] for lin in lins], np.array([1]*len(U)+[0]*len(V))[order]

def get_descendants(node):
    return set(node['children']) | set.union(*[get_descendants(c) for c in node['children']]) if len(node['children']) > 0 else set([])

def gather_groups(clusters, prevalences, count_scores = tuple([0.1, 4, 4, 1] + [0] * 256)):
    """Greedily aggregate groups into meta-groups based on some heuristics.
     :param clusters: A tuple (U,V) of sets of root nodes representing clusters (from cluster_lineages).
     :param prevalences: A map from lineage names to prevalences.
     :param count_scores: Heuristic control prioritizing different meta-group sizes.
     :return: A list of lists clustering the clusters."""
    U,V = clusters[0].copy(), clusters[1].copy()
    groups = []
    while len(V) > 0:
        groupparent = list(V)[np.argmax([
                count_scores[len(get_descendants(v) & (U|V))] * get_agg_prevalence(v, prevalences)
            for v in list(V) ])]
        descendants = get_descendants(groupparent)
        groups.append(sorted([groupparent] + list(descendants & (U|V)), key=lambda x: x['alias']))
        V = V - set([groupparent]) - descendants
        U = U - descendants
    return sorted(groups, key=lambda x: x[0]['alias'], reverse=True)
