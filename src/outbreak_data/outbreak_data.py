import sys
import requests
import warnings
import pandas as pd
import numpy as np
import json

from outbreak_data import authenticate_user

default_server = 'api.outbreak.info' # or 'dev.outbreak.info'
print_reqs = False

def _list_if_str(x):
    if isinstance(x, str): x = list(x.split(","))
    return x

def _pangolin_crumbs(pango_lin, lin_prefix=True):
    query = 'lineages=None&' if lin_prefix else ''
    return query + f'q=pangolin_lineage_crumbs:*;{pango_lin};*'

def _multiquery_to_df(data):
    return pd.concat([pd.DataFrame(v).assign(query=k) for k,v in data['results'].items()], axis=0)

def _lin_or_descendants(pango_lin, descendants, lineage_key, join=',', exclude=[]):
    if descendants:
        if lineage_key and pango_lin in lineage_key: pango_lin = lineage_key[pango_lin]['alias']
        if not lineage_key: warnings.warn('without the lineage_key parameter, descendant queries on aliased lineages with aliased children (eg JN.1 and KP.1) will not be accurate.')
        query = _pangolin_crumbs(pango_lin)
        for ex in exclude:
            if lineage_key and ex in lineage_key: ex = lineage_key[ex]['alias']
            query += f' AND NOT pangolin_lineage_crumbs:*;{ex};*'
    else: query = f'pangolin_lineage={join.join(_list_if_str(pango_lin))}'
    return query

def _lboolstr(b):
    return str(bool(b)).lower()

def _get_user_authentication():
    try: token = authenticate_user.get_authentication()
    except:
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    if token == "":
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    return {'Authorization': 'Bearer ' + token}

def _get_outbreak_data(endpoint, argstring, server=None, auth=None, collect_all=False, curr_page=0):
    """Get data via GET from the outbreak.info API, which is based on ElasticSearch.
     :param endpoint: target index or service, specified as a URL.
     :param argstring: URL-formatted args and query (endpoint specific).
     :param server: Address of a server hosting the outbreak.info API to use for the request.
     :param auth: The authorization key to use for the request.
     :param collect_all: if True, use paging mechanism to retrieve data.
     :param curr_page: iterator state for paging.
     :return: A request object containing the endpoint's response."""
    if server is None: server = default_server
    if auth is None: auth = _get_user_authentication()
    if collect_all: argstring += ('&' if len(argstring) > 0 else '') + 'fetch_all=true'
    url = f'https://{server}/{endpoint}?{argstring}'
    if print_reqs: print('GET', url)
    in_req = requests.get(url, headers=auth)
    if in_req.headers.get('content-type') != 'application/json; charset=UTF-8':
        raise ValueError('Warning!: Potentially missing endpoint. Data not being returned by server.')
    if 400 <= in_req.status_code <= 499:
        raise NameError(f'Request error (client-side/Error might be endpoint): {in_req.status_code}')
    elif 500 <= in_req.status_code <= 599:
        raise NameError(f'Request error (server-side): {in_req.status_code}')
    json_data = in_req.json()
    if collect_all:
        json_data = {k: v if isinstance(v, list) else [v] for k, v in json_data.items()}
        if 'hits' in json_data.keys() or 'results' in json_data.keys():
            scroll_id = json_data['_scroll_id'][0]
            to_scroll = 'scroll_id=' + scroll_id + '&fetch_all=true&page=' + str(curr_page)
            next_page = _get_outbreak_data( endpoint, to_scroll, server=server, auth=auth,
                                        collect_all=True, curr_page=curr_page+1 )
            for k in json_data.keys(): json_data[k].extend(next_page.get(k) or [])
    return json_data

def mutation_details(mutations, **req_args):
    """Get details of one or more mutations from clinical data.

     :param mutations: Mutation or list of mutations of interest.

     :return: A pandas dataframe with one row per mutation.

     :Parameter example: {'mutations':'s:e484*'} """
    mutations = ','.join(_list_if_str(mutations))
    data = _get_outbreak_data('genomics/mutation-details', f'mutations={mutations}', collect_all=False, **req_args)
    return pd.DataFrame(data['results']).set_index('mutation')
def wildcard_mutations(search, **req_args):
    return mutation_details(search, **req_args)

def wildcard_lineage(search, **req_args):
    """Find pango lineages via wildcard search.

     :param search: A query string for the lineage name; asterisks treated as wildcards.

     :return: A pandas dataframe containing matching lineages and clinical sequence counts.

     :Parameter example: { 'search': '*.86.*' } """
    data = _get_outbreak_data('genomics/lineage', f'name={search}', collect_all=False, **req_args)
    return pd.DataFrame(data['results']).set_index('name')

def wildcard_location(search, **req_args):
    """Find location info via wildcard search.

     :param search: A query string for the location name; asterisks treated as wildcards.

     :return: A pandas dataframe containing matching locations and their metadata.

     :Parameter example: { 'search': '*awai*' } """
    data = _get_outbreak_data('genomics/location', f'name={search}', collect_all=False, **req_args)
    return pd.DataFrame(data['results']).set_index('id')
def location_details(location, **req_args):
    return wildcard_location(location, **req_args)

def cases_by_location(location, pull_smoothed=0, **req_args):
    """Get case counts over time in a location

     :param location: String or list of location IDs

     :param pull_smoothed: 0 -> unsmoothed data, 1 -> weekly smoothed data, 2 -> both.

     :return: A pandas df of case counts indexed by location and date.

     :Parameter example: { 'location': ['USA_US-HI', 'USA_US-KY'], 'pull_smoothed': 2 } """
    location = _list_if_str(location)
    if not isinstance(location, list) or len(location) == 0:
        raise ValueError('Please enter at least 1 valid location id')
    location = ' OR '.join(location)
    smooth_vals = ['confirmed_numIncrease', 'confirmed_rolling']
    smooth_vals += [', '.join(smooth_vals)]
    if isinstance(pull_smoothed, int) and pull_smoothed in [0, 1, 2]:
        pull_smoothed = smooth_vals[pull_smoothed]
    elif not pull_smoothed in smooth_vals: raise Exception("invalid parameter value for pull_smoothed!")
    args = f'q=location_id:({location})&sort=date&fields=date,admin1,{pull_smoothed}'
    data = _get_outbreak_data('covid19/query', args, auth={}, collect_all=True, **req_args)
    data = pd.DataFrame(data['hits']).drop(columns=['_score', 'admin1'], axis=1)
    data['location'] = data.apply(lambda x: x['_id'].split(x['date'])[0], axis=1)
    return data.set_index(['location', 'date'])[pull_smoothed.split(', ')].sort_index()

def most_recent_cl_data(pango_lin, mutations=None, location=None, submission=False, **req_args):
    """Get most recent date of clinical sequencing data by location.

     :param pango_lin: A string or list of pango lineages. If a list, behavior is OR.
     :param mutations: A string or list of mutations. If a list, behavior is OR.
     :param location: A location ID. If not specified, global data are returned.
     :param submission: True -> submission dates; False -> collection dates.

     :return: The most recent date in YYYY-MM-DD.

     :Parameter example: { 'pango_lin': 'ba.2.86.1', 'location': 'USA_US-HI' } """
    query = ''
    if pango_lin is not None: query += f'&pangolin_lineage={",".join(_list_if_str(pango_lin))}'
    if location is not None: query += f'&location_id={location}'
    if mutations is not None: query += f'&mutations={",".join(_list_if_str(mutations))}'
    date_type = 'submission' if submission else 'collection'
    data = _get_outbreak_data(f'genomics/most-recent-{date_type}-date-by-location', query[1:], collect_all=False, **req_args)
    return pd.to_datetime(pd.DataFrame([data['results']])['date'][0])
def collection_date(**args):
    return most_recent_cl_data(**args, submission=False)
def submission_date(**args):
    return most_recent_cl_data(**args, submission=True)

def daily_lag(location=None, **req_args):
    """Get the daily lag between collection and submission dates of clinical sequences in a location.

     :param location: A string containing a location ID. If not specified, return lag globally.

     :return: A pandas dataframe of collection-submission date pairs and clinical sequence counts.

     :Parameter example: { 'location': 'USA_US-HI' } """
    query = f'location_id={location}' if location is not None else ''
    data = _get_outbreak_data('genomics/collection-submission', query, collect_all=False, **req_args)
    return pd.DataFrame(data['results']).set_index(['date_collected', 'date_submitted'])

def sequence_counts(location=None, sub_admin=False, **req_args):
    """Get the number of clinical sequences collected for a location.

     :param location: A string containing a location ID. If not specified, the global total counts are returned.
     :param subadmin: False -> daily counts, True -> sub-location counts.

     :return: A pandas dataframe of sequence counts.

     :Parameter example 1: { 'location': 'USA_US-HI' }
     :Parameter example 2: { 'location': 'MEX', 'sub_admin': True } """
    query = f'&cumulative={_lboolstr(sub_admin)}&subadmin={_lboolstr(sub_admin)}'
    if location is not None: query += f'&location_id={location}'
    data = pd.DataFrame(_get_outbreak_data('genomics/sequence-count', query, **req_args)['results'])
    return data.set_index('location_id' if sub_admin else 'date')

def known_mutations(pango_lin=None, descendants=False, mutations=None, lineage_key=None, freq=0.8, **req_args):
    """Get information about each mutation present in a lineage or set of lineages in clinical sequences.

     :param pango_lin: A string or list of lineage names. Return mutations occuring in any of these lineages.
     :param descendants: If True, return mutations contained in pango_lin as well as any descendants (works only with single pango_lin).
     :param mutations: A string or list of mutation names. Return only mutations co-occuring with all of these mutations.
     :param freq: A frequency threshold above which to return mutations.

     :return: A pandas dataframe of mutation information.

     :Parameter example 1: { 'pango_lin': 'BA.2.86.1', 'descendants': True }
     :Parameter example 2: { 'pango_lin': ['BA.1', 'BA.2'] } """
    query = _lin_or_descendants(pango_lin, descendants, lineage_key, join=' OR ')
    if mutations is not None: query += f'&mutations={" AND ".join(_list_if_str(mutations))}'
    query += f'&frequency={freq}'
    data = _get_outbreak_data('genomics/lineage-mutations', query, collect_all=False, **req_args)
    return _multiquery_to_df(data).drop(columns=['query']).set_index('mutation')
def lineage_mutations(**kwargs):
    return known_mutations(**kwargs)

def mutation_prevalences( mutations=None, location=None, pango_lin=None, datemin=None, datemax=None, **req_args ):
    """Get the prevalence of a set of mutations given in some subset of clinical sequences.

     :param mutations: List of mutations to query for. When datemin and datemax aren't specified, unions and intersections of lineages may be specified using OR and AND respectively. For example ['{mut1} AND ({mut2} OR {mut3})', '{mut4}']
     :param location: The ID string of a location to query within.
     :param pango_lineage: The name of a pangolin lineage to query within.
     :param datemin: (Optional). String containing start of date range to query within in YYYY-MM-DD. Specifying datemin is only supported for AND queries; lineage information is not returned.
     :param datemax: (Optional). String containing end of date range to query within in YYYY-MM-DD. Specifying datemax is only supported for AND queries; lineage information is not returned.

     :return: A pandas dataframe of mutation information.

     :Parameter example: { 'mutations': ['orf1b:r1315c', 's:l24s'], 'pango_lin': 'BA.2' } """
    if datemin is not None or datemax is not None:
        if ' OR ' in str(mutations) or ' , ' in str(mutations):
            raise ValueError('When datemin or datemax is specified, only AND queries are supported.')
        return lineage_cl_prevalence(pango_lin='.', descendants=True, location=location, mutations=mutations, datemin=datemin, datemax=datemax, lineage_key=dict())
    query = f'mutations={", ".join(_list_if_str(mutations))}'
    if pango_lin is not None: query += f'&pangolin_lineage={pango_lin}'
    if location is not None: query += f'&location_id={location}'
    if datemin is not None: query += f'&min_date={datemin}'
    if datemax is not None: query += f'&max_date={datemax}'
    data = _get_outbreak_data('genomics/mutations-by-lineage', query, **req_args)
    return _multiquery_to_df(data).set_index('query')
def mutations_by_lineage(**kwargs):
    return mutation_prevalences(**kwargs)

def lineage_cl_prevalence( pango_lin, descendants=False, location=None, mutations=None,
                           datemin=None, datemax=None, cumulative=False, lineage_key=None, exclude_descendants=[], **req_args ):
    """Get the daily prevalence of a set of lineages in clinical sequencing data.

     :param pango_lin: Lineage name or list of lineage names to query for. When descendants=False, unions of lineages may be specified using OR, for example ['{lin1} OR {lin2}', '{lin3}']
     :param descendants (Optional): If True, return mutations contained in pango_lin as well as any descendants (works only with single pango_lin).
     :param location (Optional): A string containing the location ID to query within.
     :param mutations (Optional): A list of mutation names; query within the subset of sequences containing all of these.
     :param datemin (Optional): String containing start of date range to query within in YYYY-MM-DD.
     :param datemax (Optional): String containing end of date range to query within in YYYY-MM-DD.
     :param cumulative (Optional): If true returns the cumulative global prevalence since the first day of detection.
     :param lineage_key (Optional): The lineage key for dealiasing variant names

     :return: A pandas dataframe containing prevalence data.

     :Parameter example: { 'pango_lin': 'BA.2.86.1', 'descendants': True } """
    if len(exclude_descendants) > 0: descendants = True
    query = _lin_or_descendants(pango_lin, descendants, lineage_key, exclude=exclude_descendants)
    if location is not None: query += f'&location_id={location}'
    if mutations is not None: query += f'&mutations={" AND ".join(_list_if_str(mutations))}'
    query += f'&cumulative={_lboolstr(cumulative)}'
    if datemin is not None: query += f'&min_date={datemin}'
    if datemax is not None: query += f'&max_date={datemax}'
    try:
        data = _get_outbreak_data('genomics/prevalence-by-location', query, collect_all=False, **req_args)
        return pd.DataFrame(data['results']) if cumulative else _multiquery_to_df(data).set_index(['date'])
    except KeyError:
        print(f' No results for lineage "{pango_lin}" could be found for this location.')
def prevalence_by_location(pango_lin, **kwargs):
    return lineage_cl_prevalence(pango_lin, **kwargs)
def global_prevalence(pango_lin, **kwargs):
    return lineage_cl_prevalence(pango_lin, location=None, **kwargs)

def lineage_by_sub_admin(pango_lin, mutations=None, location=None, ndays=180, detected=False, **req_args):
    """Get clinical data from the most recent date with more than zero sequences for each sublocation.

     :param pangolin_lineage: A list or string of lineage names to query for. Results for each lineage are returned on separate rows. Queries such as ['{lin1} OR {lin2}', '{lin3}'] may be used to get results for the union of lin1 and lin2.
     :param mutations: A list or string of mutations to query for; only sequences with all of these mutations will match against pangolin_lineage.
     :param location_id: A string containing a location ID. If not specified, returns global data at the country level.
     :param ndays: A positive integer number of days back from the current date to calculative cumuative counts within.
     :param detected: If true return a list of locations where at least one matching sequence has been detected.

     :return: A pandas dataframe containing sequence count data.

     :Parameter example: { 'pango_lin': ['BA.1', 'BA.2'], 'location': 'USA' } """
    query = f'pangolin_lineage={",".join(_list_if_str(pango_lin))}'
    if mutations is not None: query += f'&mutations={" AND ".join(_list_if_str(mutations))}'
    if location is not None: query += f'&location_id={location}'
    if ndays is not None: query += f'&ndays={ndays}'
    query += f'&detected={_lboolstr(detected)}'
    data = _get_outbreak_data('genomics/lineage-by-sub-admin-most-recent', query, collect_all=False, **req_args)
    return _multiquery_to_df(data).set_index(['name', 'query'])

def all_lineage_prevalences( location=None, ndays=180, nday_threshold=10, other_threshold=0.05,
                             other_exclude=None, cumulative=False, **req_args ):
    """Get prevalences of lineages circulating in a location according to clinical sequencing data.

     :param location: A string containing a location ID. If not specified, global data is returned.
     :param other_threshold: Minimum prevalence threshold below which lineages will be aggregated under "other".
     :param nday_threshold: Minimum number of days in which a lineage's prevalence must be above other_threshold in order to not be aggregated.
     :param ndays: The number of days before the current date to be used as a window to accumulate lineages under "other".
     :param other_exclude: List of lineages that are not to be included under "other".
     :param cumulative: If true return the cumulative prevalence; otherwise return daily data.

     :return: A pandas dataframe containing lineage prevalences.

     :Parameter example: { 'location': 'USA_US-HI' }
     :Parameter example: { 'cumulative': True } """
    query = ''
    if location is not None: query += f'&location_id={location}'
    query += f'&ndays={ndays}&nday_threshold={nday_threshold}&other_threshold={other_threshold}&cumulative={_lboolstr(cumulative)}'
    if other_exclude is not None: query += f'&other_exclude={",".join(_list_if_str(other_exclude))}'
    data = _get_outbreak_data('genomics/prevalence-by-location-all-lineages', query[1:], **req_args)
    data = pd.DataFrame(data['results'])
    data['lineage'] = data['lineage'].str.upper()
    return data.set_index('lineage' if cumulative else ['date', 'lineage'])

def growth_rates(lineage, location='Global'):
    """Get growth rate data for a given lineage in a given location.

     :param lineage: A list or string of lineage names.
     :param location: A list or string of location IDs.

     :return: A pandas dataframe of lineage growth data.

     :Parameter example: { 'lineage': ['XBB.1.5', 'BA.2.86'], 'location': ['Global', 'USA'] } """
    query = f'q=lineage:({" OR ".join(_list_if_str(lineage))}) AND location:({" OR ".join(_list_if_str(location))})'
    data = _get_outbreak_data('growth_rate/query', query, collect_all=False)
    return pd.concat([ pd.DataFrame(d['values'])
                         .assign(lineage = d['lineage'])
                         .assign(location = d['location']) for d in data['hits'] ], axis=0).set_index(['location', 'lineage', 'date'])

def gr_significance(location='Global', n=5):
    """Get the top lineages with the most significant growth behavior in a given location.

     :param location: List or string of location IDs.
     :param n: Number of lineages to return.

     :return: A pandas dataframe of significant lineages.

     :Parameter example: { 'location': ['USA', 'Global'] } """
    query = f'q=loc:({" OR ".join(_list_if_str(location))}) AND growing:true&sort=-sig&size=5'
    data = _get_outbreak_data('significance/query', query, collect_all=False)
    return pd.DataFrame(data['hits']).set_index('lin')
    
def _ww_metadata_query( country=None, region=None, collection_site_id=None,
                        date_range=None, sra_ids=None, viral_load_at_least=None,
                        population_at_least=None, demix_success=True, variants_success=True, **kwargs ):
    query_params = []
    if country is not None:
        query_params.append(f"geo_loc_country:{country}")
    if region is not None:
        query_params.append(f"geo_loc_region:{region}")
    if collection_site_id is not None:
        query_params.append(f"collection_site_id:{collection_site_id}")
    if date_range is not None:
        query_params.append(f"collection_date:[{date_range[0]} TO {date_range[1]}]")
    if sra_ids is not None:
        sra_query = " OR ".join([f"sra_accession:{sra_id}" for sra_id in sra_ids])
        query_params.append(f"({sra_query})")
    if viral_load_at_least is not None:
        query_params.append(f"viral_load:>={viral_load_at_least}")
    if population_at_least is not None:
        query_params.append(f"ww_population:>={population_at_least}")
    if demix_success is not None:
      query_params.append(f"demix_success:{_lboolstr(demix_success)}")
    if variants_success is not None:
      query_params.append(f"variants_success:{_lboolstr(variants_success)}")
    return " AND ".join(query_params)

def _get_ww_results(data):
    try: return pd.DataFrame(data['hits'])
    except: raise KeyError("No data for query was found.")

def _normalize_viral_loads_by_site(df):
    site_vars = df.groupby('collection_site_id', observed=True)['viral_load'].std(ddof=1).rename('site_var')
    site_vars = site_vars.reindex(df['collection_site_id'])
    site_vars.index = df.index
    normed_vl = df['viral_load'] / site_vars
    return normed_vl.where(np.isfinite(normed_vl), pd.NA)

def get_wastewater_latest(**kwargs):
    """Get date of latest wastewater sample matching a given query. Same parameters as `get_wastewater_samples`.

     :return: The date of the most recent matching sample in YYYY-MM-DD.

     :Parameter example: { 'region': 'Ohio', 'server': 'dev.outbreak.info' } """
    query = _ww_metadata_query(**kwargs)
    data = _get_outbreak_data( 'wastewater_metadata/query',
        "size=1&sort=-collection_date&fields=collection_date&q=" + query,
        server=kwargs.get('server'), auth=kwargs.get('auth') )
    return _get_ww_results(data)['collection_date'][0]

def get_wastewater_samples(**kwargs):
    """Get IDs and metadata of wastewater samples matching a given query.

     :param country: String containing name of country to search within.
     :param region: String containing name of region to search within.
     :param collection_site_id: ID of collection site.
     :param date_range: Date range in the format [start_date, end_date], with dates in YYYY-MM-DD.
     :param sra_ids: List of sample IDs.
     :param viral_load_at_least: Minimum viral load threshold for matching samples.
     :param population_at_least: Minimum population threshold for matching samples.
     :param demix_success: Whether to gather only samples with valid lineage mix data.
     :param variants_success: Whether to gather only samples with valid mutation data.

     :return: A pandas dataframe containing the IDs and metadata of matching samples.

     :Parameter example: { 'region': 'Ohio', 'date_range': ['2023-06-01', '2023-12-31'], 'server': 'dev.outbreak.info' } """
    query = _ww_metadata_query(**kwargs)
    data = _get_outbreak_data( 'wastewater_metadata/query', f"q=" + query,
                              collect_all=True, server=kwargs.get('server'), auth=kwargs.get('auth'))
    df = _get_ww_results(data).drop(columns=['_score', '_id'])
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    df['normed_viral_load'] = _normalize_viral_loads_by_site(df)
    return df.set_index('collection_date')

def get_wastewater_samples_by_lineage(lineage, descendants=False, min_prevalence=0.01, **req_args):
    """Get IDs of wastewater samples containing a certain lineage.

     :param lineage: String containing the name of the target lineage.
     :param descendants: If true, include that lineage's descendants in the query.
     :param min_prevalence: The minimum prevalence necessary for a sample to be considered to contain a lineage.

     :return: A pandas series containing IDs of samples found to contain matching lineages.

     :Parameter example: { 'lineage': 'EG.5.1', 'server': 'dev.outbreak.info' } """
    namequery = f'name:{lineage}' if not descendants else f'crumbs:*;{lineage};*'
    data = _get_outbreak_data('wastewater_demix/query', f"q=prevalence:>={min_prevalence} AND {namequery}", collect_all=True, **req_args)
    data = _get_ww_results(data).drop(columns=['_score', '_id'])
    return data.set_index(pd.Index([lineage]*len(data)))

def get_wastewater_samples_by_mutation(site, alt_base=None, min_prevalence=0.01, **req_args):
    """Get IDs of wastewater samples containing a mutation at a certain site.

     :param site: Positive integer representing the base pair index of mutations of interest.
     :param alt_base: The new base at that site (from ['G', 'A', 'T', 'C']).
     :param min_prevalence: The minimum prevalence necessary for a sample to be considered to contain a mutation.

     :return: A pandas series containing IDs of samples found to contain matching mutations.

     :Parameter example: { 'site': 1003, 'alt_base': 'G', 'server': 'dev.outbreak.info' } """
    alt_base = '' if alt_base is None else ' AND alt_base:' + alt_base
    data = _get_outbreak_data('wastewater_variants/query', f"q=prevalence:>={min_prevalence} AND site:{str(site)}{alt_base}", collect_all=True, **req_args)
    data = _get_ww_results(data).drop(columns=['_score', '_id'])
    data['mutation'] = str(site) + str(alt_base)
    return data.set_index('mutation')

def _fetch_ww_data(sample_metadata, endpoint, server=None, auth=None):
    if server is None: server = default_server
    if auth is None: auth = _get_user_authentication()
    if not isinstance(sample_metadata, pd.DataFrame): sample_metadata = pd.Series(sample_metadata).rename('sra_accession').to_frame()
    data = {"q": sample_metadata['sra_accession'].unique().tolist(), "scopes": "sra_accession"}
    url = f'https://{server}/{endpoint}/?size=1000'
    if print_reqs: print('POST', url)
    response = requests.post(url, headers=auth, json=data)
    if not response.ok:
        raise RuntimeError('Request failed. Please check that the network connection and endpoint are online.')
    df = pd.DataFrame(response.json())
    if not '_score' in df.columns:
        raise RuntimeError('Empty response. Please check the query.')
    df = df.drop(columns=['_score', '_id'])
    merged_data = pd.merge(sample_metadata.reset_index(names=sample_metadata.index.names), df, on='sra_accession').set_index(sample_metadata.index.names)
    return merged_data.drop(columns='notfound', errors='ignore')

def get_wastewater_metadata(input_df, **req_args):
    """Add wastewater sample metadata to a DataFrame containing sample IDs.

     :param input_df: Pandas DataFrame containing sample IDs as a column, as from get_wastewater_samples_by_*. A list of accession IDs is also supported.

     :return: The input dataframe joined with metadata columns.

     :Parameter example: { 'input_df': ['SRR26963071', 'SRR25666039'], 'server': 'dev.outbreak.info' } """
    if isinstance(input_df, pd.DataFrame) and 'geo_loc_country' in input_df.columns:
        raise ValueError('This DataFrame already seems to have metadata information.')
    df = _fetch_ww_data(input_df, 'wastewater_metadata/query', **req_args)
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    df['normed_viral_load'] = _normalize_viral_loads_by_site(df)
    return df.set_index('collection_date', append=True).reorder_levels([1, 0])

def get_wastewater_mutations(input_df, **req_args):
    """Add wastewater mutations data to a DataFrame containing sample IDs.

     :param input_df: Pandas DataFrame containing sample IDs as a column, as from get_wastewater_samples_by_*. A list of accession IDs is also supported.

     :return: The input dataframe joined with mutation data columns.

     :Parameter example: { 'input_df': ['SRR26963071', 'SRR25666039'], 'server': 'dev.outbreak.info' } """
    if isinstance(input_df, pd.DataFrame) and 'mutation' in input_df.columns:
        raise ValueError('This DataFrame already seems to have mutation information.')
    data = _fetch_ww_data(input_df, 'wastewater_variants/query', **req_args)
    data['mutation'] = data['site'].astype(int).astype(str) + data['alt_base'].astype(str)
    return data.set_index('mutation', append=True)

def get_wastewater_lineages(input_df, **req_args):
    """Add wastewater demix results to a DataFrame containing sample IDs.

     :param input_df: Pandas DataFrame containing sample IDs as a column, as from get_wastewater_samples_by_*. A list of accession IDs is also supported.

     :return: The input dataframe joined with lineage data columns.

     :Parameter example: { 'input_df': ['SRR26963071', 'SRR25666039'], 'server': 'dev.outbreak.info' } """
    if isinstance(input_df, pd.DataFrame) and ('name' in input_df.columns or 'lineage' in input_df.columns):
        raise ValueError('This DataFrame already seems to have lineage information.')
    data = _fetch_ww_data(input_df, 'wastewater_demix/query', **req_args)
    return data.rename(columns={'name': 'lineage'}).set_index('lineage', append=True)

def get_wastewater_lin_prevalences(**kwargs):
    """Get aggregated lineage prevalences from ww. See get_wastewater_samples for parameters."""
    server = kwargs['server'] if 'server' in kwargs else None
    kwargs['demix_success'] = True
    samples=get_wastewater_samples(kwargs)
    samples=get_wastewater_lineages(samples)
    return datebin_and_agg(samples)

def infer_mutations(mutation_df, muts_of_interest):
    """Infer which samples contain mutations with zero prevalence based on coverage data.

     :param mutation_df: A multi-indexed pandas dataframe of mutations; df.index[1] is assumed to contain mutation names.
     :param muts_of_interest: The list of mutations to infer zero-prevalence samples of.

     :return: The input df sliced down to `muts_of_interest` with additional rows corresponding to zero-mutation-prevalence samples added."""
    mutation_df = mutation_df.copy().loc[pd.IndexSlice[:, muts_of_interest],:]
    mutation_df = mutation_df.set_index(mutation_df['sra_accession'], append=True).unstack(1).stack(dropna=False)
    mutation_df_b = mutation_df.reset_index(level=1, drop=True)
    mutation_df = mutation_df_b.interpolate().bfill().ffill()
    mutation_df['prevalence'] = mutation_df_b['prevalence']
    def is_covered(x):
        for i in x['coverage_intervals']:
            if i['start'] <= x['site'] and x['site'] <= i['end']: return True
        return False
    mutation_df[mutation_df.apply(is_covered, axis=1)]
    mutation_df['prevalence'] = mutation_df['prevalence'].fillna(0)
    return mutation_df

def get_wastewater_mut_prevalences(mutations, **kwargs):
    """Get prevalences of a list of mutations. See get_wastewater_samples for parameters."""
    server = kwargs['server'] if 'server' in kwargs else None
    kwargs['variants_success'] = True
    samples = get_wastewater_samples(kwargs)
    samples = get_wastewater_mutations(samples)
    samples = infer_mutations(samples, mutations)
    return datebin_and_agg(samples)