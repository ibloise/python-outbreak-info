import numpy as np
import frozendict
import requests
import gzip
import json
from collections import OrderedDict

def get_compressed_tree(url='https://raw.githubusercontent.com/outbreak-info/python-outbreak-info/new_docs/tree.json.gz'):
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

def cluster_lineages(prevalences, tree, lineage_key=None, n=10, K=set([]), alpha=0.15):
    """Cluster lineages via greedy group-splitting on the phylo tree starting from the root based on some heuristics.

     :param prevalences: A dict, pandas series or other map between lineage names and (un-normalized) prevalences.
     :param tree: A frozendict representing the root of the phylo tree object.
     :param lineage_key: An OrderedDict mapping names to tree nodes.
     :param n: The target number of clusters.
     :param alpha: Heuristic control in range (0, 1); higher values avoid more low-quality groups, but can prevent convergence on some data.

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
    cs = set([])
    for c in K:
        update_ancestors(c, -agg_prevalences[c['lindex']], cs)
        cs = cs | set([c])
    root = tree
    while root['parent'] != root['name']:
        root = lineage_key[root['parent']]
    while len(U|V) < n:
        add_node_candidates = [c for w in U|V for c in w['children'] if not c in U|V|K]
        score = lambda c: agg_prevalences[c['lindex']] * agg_prevalences[lineage_key[c['parent']]['lindex']]
        add_node = add_node_candidates[np.argmax([score(c) for c in add_node_candidates])]
        split_node = lineage_key[add_node['parent']]
        update_ancestors(add_node, -agg_prevalences[add_node['lindex']], U|V|K)
        if split_node in U:
            U = U - set([split_node])
            V = V | set([split_node])
        if contains_descendant(add_node, U|V):
                V = V | set([add_node])
        else:   U = U | set([add_node])
        if len(U) > 1 and len(list(V - set([tree]))) > 1:
            drop_node_candidates = list(V - set([tree]))
            drop_node = drop_node_candidates[np.argmin([agg_prevalences[d['lindex']] for d in drop_node_candidates])]
            if agg_prevalences[drop_node['lindex']] < alpha * np.mean([agg_prevalences[u['lindex']] for u in U]):
                V = V - set([drop_node])
                update_ancestors(drop_node, agg_prevalences[drop_node['lindex']], U|V|K)
    if root != tree: K = K|set([root])
    return U,V,K

def get_agg_prevalence(root, prevalences, W=set([])):
    """Compute the total prevalence of all lineages in a group.

     :param root: A tree node representing the root of the group (a frozendict).
     :param prevalences: A map from lineage names to prevalences.
     :param W: A set of tree nodes representing the roots of other groups. All nodes descended from root that are in W or are children of those nodes are excluded.

     :return: The prevalence of the group."""
    cs = [get_agg_prevalence(c, prevalences, W) for c in root['children'] if not c in W]
    return np.clip((prevalences[root['name']] if root['name'] in prevalences else 0) + np.sum(cs), 0, None)

def get_agg_prevalences(roots, prevalences):
    """Get exclusive prevalences for all groups in a list of sets of groups."""
    return [[get_agg_prevalence(c, prevalences, set().union(*roots)) for c in sorted(list(g), key=lambda x: x['alias'], reverse=True)] for g in roots]

def get_group_root_names(roots, aliases=False):
    """Get root names for all groups in a list of sets of groups."""
    return [[c['alias' if aliases else 'name'] for c in sorted(list(g), key=lambda x: x['alias'], reverse=True)] for g in roots]

def get_descendants(node):
    """Get the set of all descendants of some node."""
    return set(node['children']) | set.union(*[get_descendants(c) for c in node['children']]) if len(node['children']) > 0 else set([])

def gather_groups(clusters, prevalences, count_scores=tuple([0.1, 4, 4, 4, 0.1] + [0] * 256)):
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
