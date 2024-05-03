import pandas as pd
import re
import warnings
import numpy as np
from frozendict import frozendict
import yaml
import requests

def normalize_viral_loads_by_site(df):
    site_vars = df.groupby('collection_site_id', observed=True)['viral_load'].std(ddof=1).rename('site_var')
    df = df.merge(site_vars, how='left', on='collection_site_id')
    df['normed_viral_load'] = df['viral_load'] / df['site_var'].fillna(df['viral_load']*2)
    df['normed_viral_load'] = df['normed_viral_load'].where(np.isfinite(df['normed_viral_load']), pd.NA)
    return df.drop(columns=['site_var'])

def datebin_and_agg(df, freq='7D', startdate=None, enddate=None, loaded=True, norm_abundance=True):
    df = df.copy()
    if startdate is None: startdate = df['collection_date'].min()
    if enddate is None: enddate = df['collection_date'].max()
    startdate = pd.to_datetime(startdate)-pd.Timedelta('1 day')
    enddate = pd.to_datetime(enddate)+pd.Timedelta('1 day')
    if freq is None: dbins = [pd.Interval(startdate, enddate)]
    else: dbins = pd.interval_range(startdate, enddate, freq=freq)
    df['date_bin'] = pd.cut(pd.to_datetime(df['collection_date']), dbins)
    df = df[~df['date_bin'].isna()]
    df['weight'] = df['normed_viral_load'].fillna(0.5) * df['ww_population'] if loaded else df['ww_population']
    agg_loads = lambda x: (x['normed_viral_load'] * x['ww_population']).sum() / x['ww_population'].sum()
    agg_abundance = lambda lin: lambda x: (x['abundance'] * (x['name'] == lin) * x['weight']).sum() / \
                                          ((x['abundance'] if norm_abundance==True else 1) * x['weight']).sum()
    bins = df.groupby('date_bin', observed=True)
    agged_loads = bins.apply(agg_loads, include_groups=False)
    agged_loads = agged_loads.rename('viral_load')
    agged_loads = agged_loads.where(np.isfinite(agged_loads), pd.NA)
    agged_abundances = [ bins.apply(agg_abundance(lin), include_groups=False).rename(lin).fillna(0) for lin in df['name'].unique() ]
    agged_abundances = pd.DataFrame(agged_abundances).T
    agged_abundances = agged_abundances.rename(columns = {c:c.split('-like')[0] for c in agged_abundances.columns})
    agged_abundances = agged_abundances.T.groupby(agged_abundances.columns).sum().T
    agged_abundances = [agged_abundances[column] for column in agged_abundances.columns]
    return pd.concat( [agged_loads] + agged_abundances, axis=1).reindex(dbins)
    
def get_tree(url='https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/curated_reports_prep/lineages.yml'):
    response = requests.get(url)
    response = yaml.safe_load(response.content.decode("utf-8"))
    lin_names = sorted(['*'] + [lin['name'] for lin in response])
    lindex = {lin:i for i,lin in enumerate(lin_names)}
    lineage_key = dict([(lin['name'], lin) for lin in response if 'parent' in lin])
    def get_children(node, lindex):
        return tuple( frozendict({ 'name': lineage_key[c]['name'], 'lindex': lindex[lineage_key[c]['name']],
                                   'alias': lineage_key[c]['alias'], 'parent': node['name'],
                                   'descendants': tuple(node['children']), 'children': get_children(lineage_key[c], lindex) })
                         for c in node['children'] if c in lineage_key and lineage_key[c]['parent'] == node['name'] )
    tree = tuple( frozendict({ 'name': lin['name'], 'lindex': lindex[lin['name']],
                               'alias': lin['alias'], 'parent': '*', 'descendants': tuple(lin['children']),
                               'children': get_children(lin, lindex) }) for lin in response if not 'parent' in lin )
    tree = frozendict({ 'name': '*', 'lindex': lindex['*'], 'alias': '*',
                        'parent': '*', 'descendants':tuple(lineage_key.keys()), 'children': tree })
    def get_names(tree):
        return np.concatenate([[(tree['name'], tree)]] + [get_names(c) for c in tree['children']])
    lineage_key = dict(sorted(get_names(tree), key=lambda x: x[0]))
    return tree, lineage_key

def cluster_lineages(tree, abundances, n=16, alpha=0.1):
    (tree, lineage_key) = tree
    abundances = np.array(abundances.reindex(lineage_key.keys()).fillna(0))
    agg_abundances = np.zeros_like(abundances) - 1
    def init_agg_abundances(node):
        if agg_abundances[node['lindex']] < 0:
            agg_abundances[node['lindex']] = abundances[node['lindex']]
            agg_abundances[node['lindex']] += np.sum([init_agg_abundances(c) for c in node['children']])
        return agg_abundances[node['lindex']]
    init_agg_abundances(tree)
    def update_ancestors(node, diff, W):
        if not node in W and node['parent'] != node['name']:
            parent = lineage_key[node['parent']]
            agg_abundances[parent['lindex']] += diff
            update_ancestors(parent, diff, W)
    def contains_descendant(node, nodeset):
        return node in nodeset or \
               len(set(node['children']) & nodeset) > 0 or \
               np.sum([contains_descendant(c, nodeset) for c in node['children']]) > 0
    U,V = set([tree]), set([])
    while len(U|V) < n:
        split_node = list(U|V)[np.argmax([ agg_abundances[w['lindex']] * 0.1 + np.max(
            [0] + [agg_abundances[c['lindex']] for c in w['children'] if not c in U|V] ) for w in list(U|V) ])]
        C = [c for c in split_node['children'] if not c in U|V]
        add_node = C[np.argmax([agg_abundances[c['lindex']] for c in C])]
        update_ancestors(add_node, -agg_abundances[add_node['lindex']], U|V)
        if split_node in U:
            U = U - set([split_node])
            V = V | set([split_node])
        if contains_descendant(add_node, U|V):
                V = V | set([add_node])
        else:   U = U | set([add_node])
        if len(U) > 1 and len(list(V - set([tree]))) > 1:
            drop_node = list(V - set([tree]))[np.argmin([agg_abundances[v['lindex']] for v in list(V - set([tree]))])]
            if agg_abundances[drop_node['lindex']] < alpha * np.mean([agg_abundances[u['lindex']] for u in U]):
                V = V - set([drop_node])
                update_ancestors(drop_node, agg_abundances[drop_node['lindex']], U|V)
    return U,V

def get_agg_abundance(lin, abundances, W=set([])):
    cs = [get_agg_abundance(c, abundances, W) for c in lin['children'] if not c in W]
    return np.clip((abundances[lin['name']] if lin['name'] in abundances else 0) + np.sum(cs), 0, None)

def cluster_df(tree, clusters, df):
    (tree, lineage_key) = tree
    (U,V) = clusters
    viral_load = None
    if 'viral_load' in df.columns:
        abundances = df.drop(columns=['viral_load']).mul(df['viral_load'], axis=0)
        viral_load = df['viral_load'].to_frame()
        df = df.drop(columns=['viral_load'])
    abundances = np.array(abundances.reindex(lineage_key.keys()).fillna(0))
    abundances = df.sum(axis=0)
    abundances_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'      {u["alias"]}*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
    vlabels = [f'other {v["alias"]}*' + (f' ({v["name"]})' if v["name"] != v["alias"] else '') for v in V]
    legend = list(np.array(ulabels+vlabels)[order])
    clustered_abundances = pd.DataFrame(
        { d: { label:get_agg_abundance(lin, a, U|V)
            for label, lin in zip(legend, lins) }
        for d,a in zip(dates, abundances_dated) } ).transpose()
    clustered_abundances[np.sum(clustered_abundances, axis=1) < 0.5] = pd.NA
    clustered_abundances['other **'] += 1 - clustered_abundances.sum(axis=1)
    clustered_abundances['other **'] = np.clip(clustered_abundances['other **'], 0, 1)
    clustered_abundances = clustered_abundances.rename_axis(viral_load.index.name)
    if viral_load is not None: clustered_abundances = clustered_abundances.join(viral_load)
    return clustered_abundances, [lin['name'] for lin in lins], np.array([1]*len(U)+[0]*len(V))[order]

def get_descendants(node):
    return set(node['children']) | set.union(*[get_descendants(c) for c in node['children']]) if len(node['children']) > 0 else set([])

def gather_groups(clusters, abundances, count_scores = tuple([0.1, 4, 4, 1] + [0] * 256)):
    U,V = clusters[0].copy(), clusters[1].copy()
    groups = []
    while len(V) > 0:
        groupparent = list(V)[np.argmax([
                count_scores[len(get_descendants(v) & (U|V))] * get_agg_abundance(v, abundances)
            for v in list(V) ])]
        descendants = get_descendants(groupparent)
        groups.append(sorted([groupparent] + list(descendants & (U|V)), key=lambda x: x['alias']))
        V = V - set([groupparent]) - descendants
        U = U - descendants
    return sorted(groups, key=lambda x: x[0]['alias'], reverse=True)