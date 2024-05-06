import sys
import requests
import warnings
import pandas as pd
import json

from . import authenticate_user

default_server = 'api.outbreak.info' # or 'dev.outbreak.info'
print_reqs = False
_dopage = 'fetch_all=true'

def _get_user_authentication():
    """Get the authorization token.
    :return token: the users authorization token"""
    try: token = authenticate_user.get_authentication()
    except:
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    if token == "":
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    return {'Authorization': 'Bearer ' + token}

def _get_outbreak_data(endpoint, argstring, server=None, auth=None, collect_all=False, curr_page=0):
    """Receives raw data using outbreak API. 
     :param endpoint: directory in server the data is stored
     :param argstring: feature arguments to provide to API call
     :param server: Server to request from
     :param auth: Auth key (defaults to acceptable state)
     :param collect_all: if True, returns all data.
     :param curr_page: iterator state for paging
     :return: A request object containing the raw data"""
    if server is None: server = default_server
    if auth is None: auth = _get_user_authentication()
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
            to_scroll = 'scroll_id=' + scroll_id + '&fetch_all=True&page=' + str(curr_page)
            next_page = _get_outbreak_data( endpoint, to_scroll, server=server, auth=auth,
                                        collect_all=True, curr_page=curr_page+1 )
            for k in json_data.keys(): json_data[k].extend(next_page.get(k) or [])
    return json_data

def _list_if_str(x):
    if isinstance(x, str):
        x = x.replace(" ", "")
        x = list(x.split(","))
    return x

def _pangolin_crumbs(pango_lin, lin_prefix=True):
    query = 'lineages=None&' if lin_prefix else ''
    return query + f'q=pangolin_lineage_crumbs:*;{pango_lin};*'

def cases_by_location(location, pull_smoothed=0, **req_args):
    """Loads data from a location if input is a string, or from multiple locations
    if location is a list of string locations. Since this API endpoint supports paging, collect_all is used to return all data.
     :param location: A string or list of strings, separate multiple locations by ","
     :param pull_smoothed: For every value >= 0, returns 1000 obs. (paging)
     :return: A pandas dataframe"""
    location = _list_if_str(location)
    if not isinstance(location, list) or len(location) == 0:
        raise ValueError('Please enter at least 1 valid location id')
    location = ' OR '.join(location)
    smooth_vals = ['confirmed_numIncrease', 'confirmed_rolling']
    smooth_vals += [', '.join(smooth_vals)]
    if isinstance(pull_smoothed, int) and pull_smoothed in [0, 1, 2]:
        pull_smoothed = smooth_vals[pull_smoothed]
    elif not pull_smoothed in smooth_vals: raise Exception("invalid parameter value for pull_smoothed!")
    args = f'q=location_id:({location})&sort=date&fields=date,admin1,{pull_smoothed}&{_dopage}'
    data = _get_outbreak_data('covid19/query', args, auth={}, collect_all=True, **req_args)
    return pd.DataFrame(data['hits']).drop(columns=['_score', 'admin1'], axis=1)

def _multiquery_to_df(data):
    return pd.concat([pd.DataFrame(v).assign(query=k) for k,v in data['results'].items()], axis=0)

def _lin_or_ancestors(pango_lin, ancestors):
    if ancestors: query = _pangolin_crumbs(pango_lin)
    else: query = f'pangolin_lineage={",".join(_list_if_str(pango_lin))}'
    return query

def lineage_mutations(pango_lin=None, ancestors=False, mutations=None, freq=0.8, **req_args):
    """Retrieves data from all mutations in a specified lineage above a frequency threshold.
       - Use 'OR' in a string to return overlapping mutations in multiple lineages: 'BA.2 OR BA.1'
     :param pango_lin: A string; loads data for all mutations in a specified PANGO lineage
     :param ancestors: If true returns data for descendant lineages of pango_lin. Include the wildcard '*' in string to return info on all related lineages.
     :param mutations: A string; loads mutation data for the specified sequence under the specified PANGO lineage 
     :param freq: A number between 0 and 1 specifying the frequency threshold above which to return mutations (default = 0.8)
     :return: A pandas dataframe"""
    query = _lin_or_ancestors(pango_lin, ancestors)
    if mutations is not None: query += f'&mutations={mutations}'
    query += f'&frequency={freq}'
    data = _get_outbreak_data('genomics/lineage-mutations', query, collect_all=False, **req_args)
    return _multiquery_to_df(data)

def sequence_counts(location=None, cumulative=None, sub_admin=None, **req_args):
    """Returns number of sequences per day by location
     :param location: (Somewhat optional). If not specified, the global total counts are returned.
     :param cumulative: (Somewhat optional). If true returns the cumulative number of sequences till date.
     :param subadmin: (Somewhat optional). If true and cumulative=true, returns the cumulative number of sequences for the immedaite lower admin level.
     :return: A pandas dataframe."""
    query = ''    
    if location is not None: query += f'location_id={location}'
    query += f'&cumulative={str(bool(cumulative)).lower()}'
    query += f'&subadmin={str(bool(sub_admin)).lower()}'
    if query.startswith('&'): query = query[1:]
    data = _get_outbreak_data('genomics/sequence-count', f'{query}', **req_args)
    return pd.DataFrame({'Values' : data['results']} if cumulative or sub_admin else data['results'])

def mutations_by_lineage( mutation=None, location=None, pango_lin=None, ancestors=False,
                          datemin=None,  datemax=None, freq=None, **req_args ):
    """Returns the prevalence of a mutation or series of mutations across specified lineages by location
     :param mutations: (Optional). List or string of mutations separated by ",". 
     :param location_id: (Optional). A string; If not specified, return most recent date globally.
     :param pangolin_lineage: (Optional). If not specfied, returns all Pango lineages containing that mutation.
     :param frequency: (Optional) Minimimum frequency threshold for the prevalence of a mutation in a lineage.
     :param datemin: (Optional). A string representing the first cutoff date for returned date. Must be in YYYY-MM-DD format and be before 'datemax'
     :param datemax: (Optional). A string representing the second cutoff date. Must be in YYY-MM-DD format and be after 'datemin'
     :return: A pandas dataframe."""
    query = ''
    if mutation is not None:
        query += f'mutations={" AND ".join(_list_if_str(mutation))}'
    if pango_lin is not None:
        query += '&' + _lin_or_ancestors(pango_lin, ancestors)
    if location is not None: query += f'&location_id={location}'
    if datemin is not None: query += f'&datemin={datemin}'
    if datemax is not None: query += f'&datemax={datemax}'
    if query.startswith('&'): query = query[1:]
    data = _get_outbreak_data('genomics/mutations-by-lineage', query, **req_args)
    return _multiquery_to_df(data)

def lineage_cl_prevalence( pango_lin, location=None, mutations=None,  datemin=None, ancestors=False, 
                           datemax=None, cumulative=None, **req_args ):
    """Returns the daily prevalence of a PANGO lineage by location.
        :param pango_lin: (Required). List of lineages separated by ,
        :param location: (Somewhat optional). Default location: USA
        :param mutations: (Somewhat optional). List of mutations separated by AND
        :param cumulative: (Somewhat optional). If true returns the cumulative global prevalence since the first day of detection.
        :param datemin: (Optional). A string representing the first cutoff date for returned date. Must be in YYYY-MM-DD format and be before 'datemax'
        :param datemax: (Optional). A string representing the second cutoff date. Must be in YYY-MM-DD format and be after 'datemin'
        :return: A pandas dataframe."""
    query = _lin_or_ancestors(pango_lin, ancestors)
    if location is not None: query += f'&location_id={location}'
    if mutations is not None: query += f'&mutations={" AND ".join(_list_if_str(mutations))}'
    query += f'&cumulative={str(bool(cumulative)).lower()}'
    if datemin is not None: query += f'&datemin={datemin}'
    if datemax is not None: query += f'&datemax={datemax}'
    data = _get_outbreak_data('genomics/prevalence-by-location', query, collect_all=False, **req_args)
    if cumulative: data = pd.DataFrame(data['results'])
    else: data = _multiquery_to_df(data)
    return data
def prevalence_by_location(pango_lin, **kwargs):
    return lineage_cl_prevalence(pango_lin, **kwargs)
def global_prevalence(pango_lin, **kwargs):
    return lineage_cl_prevalence(pango_lin, location=None, **kwargs)

def all_lineage_prevalences( location, ndays=180, nday_threshold=10, other_threshold=0.05,
                             other_exclude=None, cumulative=None, startswith=None, **req_args ):
    """Loads prevalence data from a location
     :param location: A string
     :param other_threshold (Default: 0) Minimum prevalence threshold below which lineages must be accumulated under "Other".
     :param nday_threshold (Default: 0) Minimum number of days in which the prevalence of a lineage must be below other_threshold to be accumulated under "Other".
     :param ndays (Default: 180) The number of days before the current date to be used as a window to accumulate lineages under "Other".
     :param other_exclude: Comma separated lineages that are NOT to be included under "Other" even if the conditions specified by the three thresholds above are met.
     :param cumulative: (Default: false) If true return the cumulative prevalence.:param startswith: A string; loads data for all lineages beginning with first letter(s) of name
     :return: A pandas dataframe"""     
    query = f'location_id={location}&ndays={ndays}&nday_threshold={nday_threshold}'
    query += f'&other_threshold={other_threshold}&cumulative={str(bool(cumulative)).lower()}'
    if other_exclude is not None:
        other_exclude = other_exclude.replace(" ", "")
        query += f'&other_exclude={",".join(_list_if_str(other_exclude))}'
        
    data = _get_outbreak_data('genomics/prevalence-by-location-all-lineages', query, **req_args)
    data = pd.DataFrame(data['results'])
    if startswith:
        data = data.loc[data['lineage'].str.startswith(startswith)]
    return data
    
def lineage_by_sub_admin(pango_lin, mutations=None, location=None, ndays=None, detected=False, **req_args):
    """Cumulative prevalence of a PANGO lineage by the immediate admin level of a location
        :param pangolin_lineage: (Required). A list or string. List of lineages separated by ,
        :param mutations: (Somewhat optional). A string or list of strings. Uses AND logic.
        :param location_id: (Somewhat optional). A string. If not specified, returns cumulative prevalence at the country level globally.
        :param ndays: (Somewhat optional). An integer. Specify number of days from current date to calculative cumuative counts. If not specified, there is no limit on the window.
        :param detected: (Somewhat optional). If true returns only if at least found in location
        :return: A pandas dataframe."""
    query = f'pangolin_lineage={",".join(_list_if_str(pango_lin))}'
    if mutations is not None: query += f'&mutations={" AND ".join(_list_if_str(mutations))}'
    if location is not None: query += f'&location_id={location}'
    if ndays is not None: query += f'&ndays={ndays}'
    query += f'&detected={str(bool(detected)).lower()}'
    data = _get_outbreak_data('genomics/lineage-by-sub-admin-most-recent', query, collect_all=False, **req_args)
    return _multiquery_to_df(data)

def most_recent_cl_data(pango_lin, mutations=None, location=None, submission=False, **req_args):
    """Most recent collection date by location
     :param pango_lin: A string. (Required).
     :param mutations: (Somewhat optional). A string or list of strings. Comma separated list of mutations.
     :param location: (Somewhat optional). If not specified, return most recent date globally.
     :return: A pandas dataframe."""
    query = f'pangolin_lineage={pango_lin}'
    if location is not None: query += f'&location_id={location}'
    if mutations is not None:
        query += f'&mutations={",".join(_list_if_str(mutations))}'
    date_type = 'submission' if submission else 'collection'
    data = _get_outbreak_data(f'genomics/most-recent-{date_type}-date-by-location', query, collect_all=False, **req_args)
    return pd.DataFrame([data['results']])
def collection_date(**args):
    return most_recent_cl_data(**args, submission=False)
def submission_date(**args):
    return most_recent_cl_data(**args, submission=True)

def daily_lag(location=None, **req_args):
    """Return the daily lag between collection and submission dates by location
     :param location_id: (Somewhat optional). If not specified, return lag globally.
     :return: A pandas dataframe."""
    query = f'location_id={location}' if location is not None else ''
    data = _get_outbreak_data('genomics/collection-submission', query, collect_all=False, **req_args)
    return pd.DataFrame(data['results'])

def location_details(location, **req_args):
    """Get location details using location ID.
    :param location: A string. (Required).
    :return: Some pandas dataframes."""
    data = _get_outbreak_data('genomics/location-lookup', f'id={location}', collect_all=False, **req_args)
    return pd.DataFrame([data['results']])
    
def mutation_details(mutations, **req_args):
    """ Returns details of a mutation.
     :param mutations: (Required). Comma separated list of mutations.
     :return: A pandas dataframe."""
    mutations = ','.join(_list_if_str(mutations))
    data = _get_outbreak_data('genomics/mutation-details', f'mutations={mutations}', collect_all=False, **req_args)
    return pd.DataFrame(data['results'])

def _wildcard_handler(endpoint, search, **req_args):
    data = _get_outbreak_data(endpoint, f'name={search}', collect_all=False, **req_args)
    return pd.DataFrame(data['results'])
def wildcard_lineage(search, **req_args):
    """Match lineage name using wildcards. 
    :param name: (Required). A string. Supports wildcards. Must use '*' at end of string.
    :return: A pandas dataframe."""
    return _wildcard_handler('genomics/lineage', search, **req_args)
def wildcard_location(search, **req_args):
    """Match location name using wildcards. 
    :param name: (Required). A string. Supports wildcards. Must use '*' at end of string.
    :return: A pandas dataframe."""
    return _wildcard_handler('genomics/location', search, **req_args)
def wildcard_mutations(search, **req_args):
    """Match mutation name using wildcards. 
    :param name: (Required). A string. Supports wildcards. Must use '*' at end of string.
    :return: A pandas dataframe."""
    return _wildcard_handler('genomics/mutations', search, **req_args)

def growth_rates(lineage, location='Global'):
    """Returns the growth rate for a given lineage in a given location.
     :param lineage: (Required)  A string. 
     :param location: (Required. Default: 'Global') A list or string. Separate multiple locations with ","
     :return: A pandas dataframe."""
    query = f'q=lineage:({" OR ".join(_list_if_str(lineage))})'
    query += f'AND location:({" OR ".join(_list_if_str(location))})'
    data = _get_outbreak_data('growth_rate/query', query, collect_all=False)
    return pd.concat([ pd.DataFrame(d['values'])
                         .assign(lineage = d['lineage'])
                         .assign(location = d['location']) for d in data['hits'] ], axis=0)

def gr_significance(location='Global'):
    """Returns the lineages with the most significant growth behavior in a given location.
     :param location: (Required. Default: 'Global') A list or string. Separate multiple locations with ","
     :return: A pandas dataframe."""
    query = f'location:({" OR ".join(_list_if_str(location))})'
    data = _get_outbreak_data('significance/query', query, collect_all=False)
    return pd.DataFrame(data['hits'])

_ww_metadata_endpoint = "wastewater_metadata/query"
_ww_demix_endpoint = "wastewater_demix/query"
_ww_variants_endpoint = "wastewater_variants/query"

def _ww_metadata_query( country=None, region=None, collection_site_id=None,
                        date_range=None, sra_ids=None, viral_load_at_least=None,
                        population_at_least=None, demix_success=True, **kwargs ):
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
      query_params.append(f"demix_success:{str(bool(demix_success)).lower()}")
    return " AND ".join(query_params)

def get_wastewater_latest(**kwargs):
    """Retrieve wastewater samples data based on specified filters.
    :param country: (Optional) Country name.
    :param region: (Optional) Region name.
    :param collection_site_id: (Optional) Site ID.
    :param viral_load_at_least: (Optional) Minimum viral load threshold.
    :return: __"""
    query = _ww_metadata_query(**kwargs)
    data = _get_outbreak_data( _ww_metadata_endpoint,
        "size=1&sort=-collection_date&fields=collection_date&q=" + query,
        server=kwargs.get('server'), auth=kwargs.get('auth') ) 
    try:
        return pd.DataFrame(data['hits'])['collection_date'][0]
    except:
        raise KeyError("No data for query was found.")

def _get_ww_results(data):
    try: return pd.DataFrame(data['hits']).drop(columns=['_score', '_id'])
    except: raise KeyError("No data for query was found.")

def get_wastewater_samples(**kwargs):
    """Retrieve wastewater samples data based on specified filters.
    :param country: (Optional) Country name.
    :param region: (Optional) Region name.
    :param collection_site_id: (Optional) Site ID.
    :param date_range: (Optional) Date range in the format [start_date, end_date].
    :param sra_ids: (Optional) List of SRA IDs.
    :param viral_load_at_least: (Optional) Minimum viral load threshold.
    :param population_at_least: (Optional) Minimum population threshold.
    :return: A pandas DataFrame containing wastewater samples data."""
    query = _ww_metadata_query(**kwargs)
    data = _get_outbreak_data( _ww_metadata_endpoint, f"{_dopage}&q=" + query,
                              collect_all=True, server=kwargs.get('server'), auth=kwargs.get('auth'))
    df = _get_ww_results(data)
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    df['ww_population'] = df['ww_population'].fillna(1000)
    return df

def get_wastewater_samples_by_lineage(lineage, descendants=False, min_abundance=0.01, **req_args):
    """Gets IDs of wastewater samples containing a certain lineage
    :param lineage: Target lineage
    :return: A pandas series containing IDs of samples found to contain that lineage"""
    namequery = f'name:{lineage}' if not descendants else f'crumbs:*;{lineage};*'
    data = _get_outbreak_data(_ww_demix_endpoint, f"{_dopage}&q=abundance:>={min_abundance} AND {namequery}", collect_all=True, **req_args)
    return _get_ww_results(data)

def get_wastewater_samples_by_mutation(site, alt_base=None, min_frequency=0.01, **req_args):
    """Gets IDs of wastewater samples containing a mutation a certain site
    :param site: Base pair index of mutations of interest
    :param alt_base: The new base at that site
    :return: A pandas series containing IDs of samples found to contain that lineage"""
    alt_base = '' if alt_base is None else ' AND alt_base:' + alt_base
    data = _get_outbreak_data(_ww_variants_endpoint, f"{_dopage}&q=frequency:>={min_frequency} AND site:{site}{alt_base}", collect_all=True, **req_args)
    return _get_ww_results(data)

def _fetch_ww_data(sample_metadata, endpoint, server=None, auth=None):
    """Retrieve and join variants or demix info with sample metadata
    :param sample_metadata: DataFrame containing metadata.
    :param endpoint: API endpoint to retrieve data.
    :return: A pandas DataFrame containing merged data with exploded nested data."""
    if server is None: server = default_server
    if auth is None: auth = _get_user_authentication()
    data = {"q": sample_metadata['sra_accession'].unique().tolist(), "scopes": "sra_accession"}
    url = f'https://{server}/{endpoint}'
    if print_reqs: print('POST', url)
    response = requests.post(url, headers=auth, json=data)
    df = pd.DataFrame(response.json()).drop(columns=['_score', '_id'])
    merged_data = pd.merge(sample_metadata, df, on='sra_accession')
    return merged_data.drop(columns='notfound', errors='ignore')

def get_wastewater_metadata(sample_metadata, **req_args):
    """Add wastewater mutations data to a DF of samples
    :param sample_metadata: DataFrame containing metadata.
    :return: A pandas DataFrame containing merged wastewater mutations data with metadata."""
    df = _fetch_ww_data(sample_metadata, _ww_metadata_endpoint, **req_args)
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    # df['ww_population'] = df['ww_population'].fillna(1000)
    return df

def get_wastewater_mutations(sample_metadata, **req_args):
    """Add wastewater mutations data to a DF of samples
    :param sample_metadata: DataFrame containing metadata.
    :return: A pandas DataFrame containing merged wastewater mutations data with metadata."""
    return _fetch_ww_data(sample_metadata, _ww_variants_endpoint, **req_args)

def get_wastewater_lineages(sample_metadata, **req_args):
    """Add wastewater demix results to a DF of samples
    :param sample_metadata: DataFrame containing metadata.
    :return: A pandas DataFrame containing merged wastewater lineage abundance data with metadata."""
    return _fetch_ww_data(sample_metadata, _ww_demix_endpoint, **req_args)
