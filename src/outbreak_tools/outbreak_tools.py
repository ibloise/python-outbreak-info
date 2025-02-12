import pandas as pd
import numpy as np
import frozendict
import requests
import gzip
import yaml
import json
import warnings
from outbreak_data import outbreak_data
import re

from outbreak_tools import outbreak_clustering

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
        
#--- borrowed from SEARCH wastewater surveillance dashboard#
def convert_rbg_to_tuple( rgb ):
    print(rgb)
    rgb = rgb.lstrip( "#" )
    return tuple( int( rgb[i :i + 2], 16 ) for i in (0, 2, 4) )
def convert_tuple_to_rgb( r, g, b ):
    return '#%02x%02x%02x' % (int(r), int(g), int(b))
def lighten_field( value, alpha, gamma=2.2 ):
    return pow( pow(255, gamma) * (1 - alpha) + pow( value, gamma ) * alpha, 1 / gamma)
def lighten_color( r, g, b, alpha, gamma=2.2 ):
    return lighten_field(r, alpha, gamma ), lighten_field( g, alpha, gamma ), lighten_field( b, alpha, gamma )

def cluster_df(df, clusters, tree, lineage_key=None, norm=True, cmap = None):
    """Aggregate the columns of a dataframe into some phylogenetic groups.

     :param df: A dataframe of prevalence signals. Rows are assumed to be date bins and columns are assumed to be lineages.
     :param clusters: A tuple (U,V) of sets of root nodes representing clusters (from cluster_lineages).
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :param norm: Whether to assume that values in a row should sum to one.
     :param norm: Whether to assume that values in a row should sum to one.

     :return: A tuple (data,names,is_inclusive) where data is the input dataframe with aggregated and relabeled columns, names contains the names of the root lineages for each column/group, and is_inclusive indicates whether the column's root is in U or V."""
    if lineage_key is None: tree = get_lineage_key(tree)
    (U,V,K) = clusters
    # if include_K:
    #     U = U|K
    #     K = set([])
    prevalences_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'{u["alias"]}*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
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


def cluster_df_wcolors(df, clusters, tree, lineage_key=None, norm=True, cmap = None):
    """Aggregate the columns of a dataframe into some phylogenetic groups.

     :param df: A dataframe of prevalence signals. Rows are assumed to be date bins and columns are assumed to be lineages.
     :param clusters: A tuple (U,V) of sets of root nodes representing clusters (from cluster_lineages).
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :param norm: Whether to assume that values in a row should sum to one.
     :param norm: Whether to assume that values in a row should sum to one.

     :return: A tuple (data,names,is_inclusive) where data is the input dataframe with aggregated and relabeled columns, names contains the names of the root lineages for each column/group, and is_inclusive indicates whether the column's root is in U or V."""
    if lineage_key is None: tree = get_lineage_key(tree)
    (U,V,K) = clusters
    # if include_K:
    #     U = U|K
    #     K = set([])
    prevalences_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'{u["alias"]}*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
    vlabels = [f'other {v["alias"]}*' + (f' ({v["name"]})' if v["name"] != v["alias"] else '') for v in V]

    # get related sequences'
    ualiases = [w['alias'] for w in list(U)]
    valiases = [w['alias'] for w in list(V)]

    alias_subsetted = {}
    for u in ualiases:
        alias_subsetted[u] = list(np.sort([w for w in valiases if u.startswith(w)]))
    # get overlapping classes across groups
    overlaps = set()
    for v in alias_subsetted.values():
        for v2 in alias_subsetted.values():
            if v != v2:
                newOverlap = set(v) & set(v2)
                overlaps = overlaps | newOverlap
    for alias in alias_subsetted.keys():
        alias_subsetted[alias] = [u for u in alias_subsetted[alias] if u not in overlaps]
    #iteratively pass through the overlaps to get additional groups.
    overlaps = list(np.sort(list(overlaps)))
    while len(overlaps)>1:
        candidate_list = list(np.sort([w for w in overlaps[0:len(overlaps)-1] if overlaps[-1].startswith(w)]))
        if len(candidate_list)>5:
            candidate_list = candidate_list[len(candidate_list)-5:len(candidate_list)]
        alias_subsetted[overlaps[-1]] = candidate_list
        overlaps.remove(overlaps[-1])
        for ci in candidate_list:
            overlaps.remove(ci)
    if len(overlaps) ==1:
        alias_subsetted[overlaps[0]] = []
    lincolors = {}
    delta = 0.15
    from matplotlib.colors import rgb2hex
    for j,pkey in enumerate(alias_subsetted.keys()):
        parent_color = convert_rbg_to_tuple(rgb2hex(cmap(j)))
        for j0,ckey in enumerate(alias_subsetted[pkey]):
            lincolors[ckey] = convert_tuple_to_rgb( *lighten_color( *parent_color, alpha=1.0-(delta*(j0)) )) 
        if len(alias_subsetted[pkey]) == 0:
            lincolors[pkey] = cmap(j)
        else:
            lincolors[pkey] = convert_tuple_to_rgb( *lighten_color( *parent_color, alpha=1.0-(delta*(j0+1)) ))      
    legend = list(np.array(ulabels+vlabels)[order])
    clustered_prevalences = pd.DataFrame(
        { d: { label:outbreak_clustering.get_agg_prevalence(lin, a, U|V)
            for label, lin in zip(legend, lins) }
        for d,a in zip(dates, prevalences_dated) } ).transpose()
    if norm:
        clustered_prevalences[np.sum(clustered_prevalences, axis=1) < 0.5] = pd.NA
        clustered_prevalences['other **'] += 1 - clustered_prevalences.sum(axis=1)
        clustered_prevalences['other **'] = np.clip(clustered_prevalences['other **'], 0, 1)
    return clustered_prevalences, [lin['name'] for lin in lins], np.array([1]*len(U)+[0]*len(V))[order], lincolors, (ualiases+valiases)



def id_lookup(locations, max_results = 10, table = False):
    """
    Helps find location ID for use with outbreak_data.py
    Requires integration with get_outbreak_data
    :param locations: A string or list of location names
    :param max_results: Int, of how many results to return
    :param table: If True, returns all results as pandas DataFrame
    :return: location_id
    """
    warnings.filterwarnings("ignore", message='Warning!: Data has "results" but length of data is 0')
    #setting max_colwidth for showing full-name completely in table
    pd.set_option("max_colwidth", 300)
    locIds_of_interest=[]
    locIds_not_found=[]
    locations_not_found=[]
    if isinstance(locations, str):
        locations = [locations]
    #first pass of the query tries every location name as-is & collects malformed queries
    for i in locations:
        locId_arg = "name=" + i
        results = outbreak_data._get_outbreak_data('genomics/location', locId_arg)
        if results != None:
            if (len(results) == 0):
                locIds_not_found.extend([i])
            else:
                df = pd.DataFrame(results['results'])
                #print(df.columns)
                if (df.shape[0]==1):
                    locIds_of_interest.extend([df.id.unique()[0]])
                else:
                    locIds_not_found.extend([i])
        #everything matches up
        if (len(locIds_of_interest)==len(locations)):
            return locIds_of_interest
        #any locations not found require further investigation (*name* must be a catch-all?)
        if (len(locIds_of_interest)!=len(locations)):
            not_found=[]
            for i in locIds_not_found:
                locs=''.join(['*', i, '*'])
                not_found.extend([locs])
            locations_not_found = not_found
    all_hits = None
    #using genomic endpoint to parse location names and corrects malformed queries
    for i in range(0, len(locations_not_found)):
        locId = "name=" + locations_not_found[i]
        results = outbreak_data._get_outbreak_data('genomics/location', locId)
        hits = pd.DataFrame()
        if(len(results) >= 1):
            hits = pd.DataFrame(results['results'])
        if(hits.shape[0] == 0):
            next
        else:
            #replacing code with meaning
            hits.admin_level.replace(-1, "World Bank Region", inplace=True)
            hits.admin_level.replace(0, "country", inplace=True)
            hits.admin_level.replace(1, "state/province", inplace=True)
            hits.admin_level.replace(1.5, "metropolitan area", inplace=True)
            hits.admin_level.replace(2, "county", inplace=True)
            hits['full'] = hits.label + ' ' + " (" + ' ' + hits.admin_level + ' ' + ")"
            hits = hits[:max_results]
            hits.index= pd.MultiIndex.from_arrays([[locations_not_found[i].strip('*')] * len(hits.index)
                                                      ,list(hits.index)], names=['Query', 'Index'])
            if isinstance(all_hits, pd.core.frame.DataFrame):
                all_hits = all_hits.append(hits)
            else:
                all_hits = hits.copy()
            if not table:
                # ask questions about ambiguous names (one-to-many)
                print("\n")
                display(hits['full'])
                print('Int values must be entered in comma seperated format.')
                loc_sel = input("Enter the indices of locations of interest in dataframe above: ")
                if loc_sel == '':
                    raise ValueError('Input string is empty, enter single or multiple int comma seperated value/s before submitting.')
                input_locs = list(loc_sel.split(','))
                for i in range(len(input_locs)):
                    val = re.search(r' *[0-9]+', input_locs[i])
                    if (isinstance(val, re.Match)):
                        if (val.group() != ''):
                            val = int(val.group())
                            input_locs[i] = val
                all_int = all([isinstance(x, int) for x in input_locs])
                if all_int:
                    locIds_of_interest.extend(hits.iloc[input_locs, :].id)
                else:
                    print("Input entries are all not int. Please try again.")
                    print('\n')
                    display(hits['full'])
                    loc_sel = input("Enter the indices of locations of interest in dataframe above: ")
                    if loc_sel == '':
                        raise ValueError('Input string is empty, enter single or multiple int comma seperated value/s before submitting.')
                    input_locs = list(loc_sel.split(','))
                    for i in range(len(input_locs)):
                        val = re.search(r' *[0-9]+', input_locs[i])
                        if (isinstance(val, re.Match)):
                            if (val.group() != ''):
                                val = int(val.group())
                                input_locs[i] = val
                    all_int = all([isinstance(x, int) for x in input_locs])
                    if all_int:
                        locIds_of_interest.extend(hits.iloc[input_locs, :].id)
                print('\n')
    if table:
        #necessary identification
        return all_hits.loc[:, ['id', 'label', 'admin_level']]
    return locIds_of_interest
