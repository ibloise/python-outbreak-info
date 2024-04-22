import sys
import requests
import warnings
import pandas as pd
import json
import numpy as np
from frozendict import frozendict
from matplotlib.colors import hsv_to_rgb
import yaml

from outbreak_data import authenticate_user

server = 'dev.outbreak.info'  # or 'api.outbreak.info'
dopage = 'fetch_all=true'  # worth verifying that this works with newer ES versions as well
covid19_endpoint = 'covid19/query'
test_server = 'test.outbreak.info'

def check_user_authentication():
    """
    Get the authorization token.
    :return token: the users authorization token
    """
    try:
        token = authenticate_user.get_authentication()
    except:
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    if token == "":
        print("Issue retrieving token, please reauthenticate.")
        sys.exit(1)
    return(token)

def get_auth(auth=None):
    """
    Get authentication header for API requests.

    Returns:
    :return: Dictionary containing the authentication header.
    """
    if auth is None:
        auth = check_user_authentication()
    return {'Authorization': 'Bearer ' + auth}


def get_outbreak_data(endpoint, argstring, server=server, auth=None, collect_all=False, curr_page=0):
    """
    Receives raw data using outbreak API.

    Arguments: 
     :param endpoint: directory in server the data is stored
     :param argstring: feature arguments to provide to API call
     :param server: Server to request from
     :param auth: Auth key (defaults to acceptable state)
     :param collect_all: if True, returns all data.
     :param curr_page: iterator state for paging
     :return: A request object containing the raw data
    """

    auth = get_auth(auth)
    url = f'https://{server}/{endpoint}?{argstring}'
    print(url)

    in_req = requests.get(url, headers=auth)
    if in_req.headers.get('content-type') != 'application/json; charset=UTF-8':
        raise ValueError('Warning!: Potentially missing endpoint. Data not being returned by server.')
    if 400 <= in_req.status_code <= 499:
        raise NameError(f'Request error (client-side/Error might be endpoint): {in_req.status_code}')
    elif 500 <= in_req.status_code <= 599:
        raise NameError(f'Request error (server-side): {in_req.status_code}')
    in_json = in_req.json()
    # checking that the request contains data
    hits = 'hits' in in_json.keys()
    results = 'results' in in_json.keys()
    contains_data = hits | results
    if collect_all is False:
        if hits and (len(in_json['hits']) == 0):
            warnings.warn('Warning!: Data has "hits" but length of data is 0')
        elif results and (len(in_json['results']) == 0):
            warnings.warn('Warning!: Data has "results" but length of data is 0')
        return in_json
    elif collect_all and not contains_data:
        return
    elif collect_all and contains_data:
        # initial dict for collecting new json data
        data_json = {k: v if isinstance(v, list) else [v] for k, v in in_json.items()}
        del data_json['_scroll_id']
        # recursion during collection
        scroll_id = in_json['_scroll_id']
        fetching_page = '&fetch_all=True&page='
        page = fetching_page + str(curr_page)
        to_scroll = 'scroll_id=' + scroll_id + page
        in_req = get_outbreak_data(endpoint, to_scroll, server=server, collect_all=True, curr_page=curr_page+1)
        if not isinstance(in_req, type(None)):
            if hits and len(in_req['hits']) == 0:
                warnings.warn('Warning!: Recursion step has "hits" key but empty data value')
            elif results and len(in_req['results']) == 0:
                warnings.warn('Warning!: Recursion step has "results" key but empty data value')
            in_req = {k: v if isinstance(v, list) else [v] for k, v in in_req.items()}
        for k in data_json.keys():
            try:
                data_json[k].extend(in_req[k])
            except TypeError:
                continue
        return data_json


def cases_by_location(location, server=server, auth=None, pull_smoothed=0):
    """
    Loads data from a location if input is a string, or from multiple locations
    if location is a list of string locations. Since this API endpoint supports paging, collect_all is used to return all data.

    Arguments:
     :param location: A string or list of strings, separate multiple locations by ","
     :param pull_smoothed: For every value >= 0, returns 1000 obs. (paging)
     :return: A pandas dataframe
    """
    # location names can be further checked to verify validity // proper format
    if isinstance(location, str):  # Converts all location input strings into lists: best universal input 
        location = location.replace(" ", "")
        location = list(location.split(","))
    if not isinstance(location, list) or len(location) == 0:
        raise ValueError('Please enter at least 1 valid location id')
    if pull_smoothed == 0:
        confirmed='confirmed_numIncrease'
    elif pull_smoothed == 1:
        confirmed='confirmed_rolling'
    elif pull_smoothed == 2:
        confirmed='confirmed_rolling, confirmed_numIncrease'
    else:
        raise Exception("invalid parameter value for pull_smoothed!")
    try:
        locations = '(' + ' OR '.join(location) + ')'
        args = f'q=location_id:{locations}&sort=date&fields=date,{confirmed},admin1&{dopage}'
        raw_data = get_outbreak_data(covid19_endpoint, args, collect_all=True)
        df = pd.DataFrame(raw_data['hits'])
        refined_table=df.drop(columns=['_score', 'admin1'], axis=1)
        return refined_table
    
    except:
        for i in location:
            raise Exception('{} is not a valid location ID'.format(i))


def all_lineage_prevalences(location, ndays=180, nday_threshold=10, other_threshold=0.05, other_exclude=None, cumulative=None, server=server, auth=None, startswith=None):
    """
    Loads prevalence data from a location

    Arguments:
     :param location: A string
     :param other_threshold (Default: 0) Minimum prevalence threshold below which lineages must be accumulated under "Other".
     :param nday_threshold (Default: 0) Minimum number of days in which the prevalence of a lineage must be below other_threshold to be accumulated under "Other".
     :param ndays (Default: 180) The number of days before the current date to be used as a window to accumulate lineages under "Other".
     :param other_exclude: Comma separated lineages that are NOT to be included under "Other" even if the conditions specified by the three thresholds above are met.
     :param cumulative: (Default: false) If true return the cumulative prevalence.:param startswith: A string; loads data for all lineages beginning with first letter(s) of name
     :return: A pandas dataframe
    """
                
    query =  f'location_id={location}&ndays={ndays}&nday_threshold={nday_threshold}&other_threshold={other_threshold}'
   
    if cumulative:
        query = query + '&' + 'cumulative=true'
    if other_exclude:
        other_exclude = other_exclude.replace(" ", "")
        query = query + '&' + f'other_exclude={other_exclude}'
        
    lins = get_outbreak_data('genomics/prevalence-by-location-all-lineages', query, server=server)
    df = pd.DataFrame(lins['results'])
    if startswith:
        return df.loc[df['lineage'].str.startswith(startswith)]
    return df


### Helper function for dealing with all 'q' queries
def pangolin_crumbs(pango_lin, mutations=None,lin_prefix=True):
    if lin_prefix:
        query = 'lineages=None&'
    else:
        query = ''
    if mutations:
        query = query + f'mutations={mutations}&'
    query = query + f'q=pangolin_lineage_crumbs:*;{pango_lin};*'
    return query


def lineage_mutations(pango_lin=None, lineage_crumbs=False, mutations=None, freq=0.8, server=server, auth=None):  ###
    """Retrieves data from all mutations in a specified lineage above a frequency threshold.
       - Use 'OR' in a string to return overlapping mutations in multiple lineages: 'BA.2 OR BA.1'

          Arguments:
             :param pango_lin: A string; loads data for all mutations in a specified PANGO lineage
             :param lineage_crumbs: If true returns data for descendant lineages of pango_lin. Include the wildcard '*' in string to return info on all related lineages.
             :param mutations: A string; loads mutation data for the specified sequence under the specified PANGO lineage 
             :param freq: A number between 0 and 1 specifying the frequency threshold above which to return mutations (default = 0.8)
             :return: A pandas dataframe"""

    # Use strings, no reason to use list format anymore
    
    if lineage_crumbs:
        query = pangolin_crumbs(pango_lin)
                
    else:
        query = f'lineages={pango_lin}'
        if 'OR' in pango_lin:
          lineages = pango_lin.split('OR')
          query = "OR".join(lineages)
        if mutations:
            query = '&' + f'mutations={mutations}' + query 
        
    if freq!=0.8:
        query = query + f'&frequency={freq}'
    raw_data = get_outbreak_data('genomics/lineage-mutations', f'{query}', collect_all=False)
    key_list = raw_data['results']
    if len(key_list) == 0:
        raise TypeError('No matches for query found')
    
    key_list = raw_data['results']
    key_list = list(key_list)
    df = pd.DataFrame(raw_data['results'][key_list[0]])
       
    return df
    

def global_prevalence(pango_lin, mutations=None, cumulative=None, lineage_crumbs=False, server=server):
   
    """Returns the global daily prevalence of a PANGO lineage
       
       Arguments:
        :param pangolin_lineage: (Required).
        :param mutations: (Somewhat optional). Comma separated list of mutations.
        :param cumulative: (Somewhat optional). If true returns the cumulative global prevalence since the first day of detection.
        :return: A pandas dataframe."""

    if lineage_crumbs:
        query = pangolin_crumbs(pango_lin)   
    else:
        if mutations:
            if isinstance(mutations, list):
                mutations = ','.join(mutations)
            elif isinstance(mutations, str):
                 mutations = mutations.replace(" ", "")
           
        query = '' + pango_lin
          
    if mutations:
        query =  query + '&' + f'mutations={mutations}' 
    if cumulative:
        query = query + '&' + 'cumulative=true'
    if lineage_crumbs:
        # using a modified formulation to access the crumbs 
        raw_data = get_outbreak_data('genomics/prevalence-by-location', query, collect_all=False)
        key_list = raw_data['results']
        key_list = list(key_list)
        if cumulative:
            for i in key_list:
                if i == key_list[0]:
                    data = {'Values' : raw_data['results'][i]}
                    df = pd.DataFrame(data) 
                else:
                    newdf = {'Values' : raw_data['results'][i]}
                    df = pd.concat([data, newdf], sort=False)  
        else:
            for i in key_list:
                if i == key_list[0]:
                    df = pd.DataFrame(raw_data['results'][i])
                else:
                    newdf = pd.DataFrame(raw_data['results'][i]) # append each df
                    df = pd.concat([df, newdf], sort=False)  
    else:
        raw_data = get_outbreak_data('genomics/global-prevalence', f'pangolin_lineage={query}')
        if cumulative:
            data = {'Values' : raw_data['results']}
            df = pd.DataFrame(data) 
        else:
            df = pd.DataFrame(raw_data['results'])
    return df

def sequence_counts(location=None, cumulative=None, sub_admin=None, server=server):
    """Returns number of sequences per day by location

    Arguments:
     :param location: (Somewhat optional). If not specified, the global total counts are returned.
     :param cumulative: (Somewhat optional). If true returns the cumulative number of sequences till date.
     :param subadmin: (Somewhat optional). If true and cumulative=true, returns the cumulative number of sequences for the immedaite lower admin level.
     :return: A pandas dataframe.
    """
        
    query = ''    
    if location:
        query = query + f'location_id={location}'
    if cumulative:
        query = query +  '&' + 'cumulative=true'
    if sub_admin:
        query = query + '&' + 'subadmin=true'
    
    raw_data = get_outbreak_data('genomics/sequence-count', f'{query}')
     
    if cumulative or sub_admin:
        data = {'Values' : raw_data['results']}
        df = pd.DataFrame(data) 
    else:
        df = pd.DataFrame(raw_data['results'])
    return df

def mutations_by_lineage(mutation=None, location=None, pango_lin=None, lineage_crumbs=False, datemin=None,  datemax=None, freq=None, server=server):
    """Returns the prevalence of a mutation or series of mutations across specified lineages by location

    Arguments:
     :param mutations: (Optional). List or string of mutations separated by ",". 
     :param location_id: (Optional). A string; If not specified, return most recent date globally.
     :param pangolin_lineage: (Optional). If not specfied, returns all Pango lineages containing that mutation.
     :param frequency: (Optional) Minimimum frequency threshold for the prevalence of a mutation in a lineage.
     :param datemin: (Optional). A string representing the first cutoff date for returned date. Must be in YYYY-MM-DD format and be before 'datemax'
     :param datemax: (Optional). A string representing the second cutoff date. Must be in YYY-MM-DD format and be after 'datemin'
     :return: A pandas dataframe.
    """
    
    if mutation: 
        if isinstance(mutation, str):
             mutation = mutation.replace(" ", "")
             mutation = list(mutation.split(","))
        mutation = '' + ' AND '.join(mutation) + ''   
        query = f'mutations={mutation}'
            
    if pango_lin and lineage_crumbs:
        query = pangolin_crumbs(pango_lin)

    if location:
        query = query + f'&location_id={location}'
    if datemin and datemax:
        query = query + f'&datemin={datemin}&datemax={datemax}'

    raw_data = get_outbreak_data('genomics/mutations-by-lineage', f'{query}')
    
    key_list = raw_data['results']
    key_list = list(key_list)
    df = pd.DataFrame(raw_data['results'][key_list[0]])
       
    if isinstance(freq, float) and freq > 0 and freq < 1:
        return df.loc[df['prevalence'] >= freq]
    return df


def prevalence_by_location(pango_lin, location, mutations=None,  datemin=None, lineage_crumbs=False, 
                           datemax=None, cumulative=None, server=server):
    """Returns the daily prevalence of a PANGO lineage by location.
   
       Arguments:
        :param pango_lin: (Required). List of lineages separated by ,
        :param location_id: (Somewhat optional). Default location: USA
        :param mutations: (Somewhat optional). List of mutations separated by AND
        :param cumulative: (Somewhat optional). If true returns the cumulative global prevalence since the first day of detection.
        :param datemin: (Optional). A string representing the first cutoff date for returned date. Must be in YYYY-MM-DD format and be before 'datemax'
        :param datemax: (Optional). A string representing the second cutoff date. Must be in YYY-MM-DD format and be after 'datemin'
        :return: A pandas dataframe."""
    if lineage_crumbs:
        query = pangolin_crumbs(pango_lin)  
        query = query + '&' + f'location_id={location}' 
    else:
        if isinstance(pango_lin, str):
            pango_lin = pango_lin.replace(" ", "")
         
        elif isinstance(pango_lin, list):
            pango_lin = ','.join(pango_lin)

        if mutations:
            if isinstance(mutations, list):
                pass
            elif isinstance(mutations, str):
                 mutations = mutations.replace(" ", "")
                 mutations = list(mutations.split(","))
            mutations = '' + ' AND '.join(mutations) + ''   
          
        query = pango_lin + '&' + f'location_id={location}'
        
    if mutations:
        query = query + '&' + f'mutations={mutations}'
    if cumulative:
        query = query + '&' + 'cumulative=true'
    if datemin and datemax:
        query = query + f'&datemin={datemin}&datemax={datemax}'
   
    if lineage_crumbs:
        raw_data = get_outbreak_data('genomics/prevalence-by-location', query, collect_all=False)
    else:
        raw_data = get_outbreak_data('genomics/prevalence-by-location', f'pangolin_lineage={query}', collect_all=False)

    key_list = raw_data['results']
    key_list = list(key_list)
    
    if cumulative:
        for i in key_list:
            if i == key_list[0]:
                 data = {'Values' : raw_data['results'][i]}
                 df = pd.DataFrame(data) 
            else:
                newdf = {'Values' : raw_data['results'][i]}
                df = pd.concat([data, newdf], sort=False)  
    else:
        for i in key_list:
            if i == key_list[0]:
                df = pd.DataFrame(raw_data['results'][i])
            else:
                newdf = pd.DataFrame(raw_data['results'][i]) # append each df
                df = pd.concat([df, newdf], sort=False)  
                
    return df


def lineage_by_sub_admin(pango_lin, mutations=None, location=None, ndays=0, detected=None, server=server):
    """Cumulative prevalence of a PANGO lineage by the immediate admin level of a location

        Arguments:
        :param pangolin_lineage: (Required). A list or string. List of lineages separated by ,
        :param mutations: (Somewhat optional). A string or list of strings. Uses AND logic.
        :param location_id: (Somewhat optional). A string. If not specified, returns cumulative prevalence at the country level globally.
        :param ndays: (Somewhat optional). An integer. Specify number of days from current date to calculative cumuative counts. If not specified, there is no limit on the window.
        :param detected: (Somewhat optional). If true returns only if at least found in location
        :return: A pandas dataframe."""
        
    if isinstance(pango_lin, str):
        pango_lin = pango_lin.replace(" ", "")
    elif isinstance(pango_lin, list):
         pango_lin = ','.join(pango_lin)
    query = pango_lin
         
    if mutations:
        if isinstance(mutations, list):
            pass
        elif isinstance(mutations, str):
             mutations = mutations.replace(" ", "")
             mutations = list(mutations.split(","))
        mutations = '' + ' AND '.join(mutations) + ''   
    
    if mutations:
        query = '' + pango_lin + '&' + f'mutations={mutations}'
    if location:
        query = query + '&' + f'location_id={location}'
    if ndays > 0:
        query = query + '&' + f'ndays={ndays}'
        
    raw_data = get_outbreak_data('genomics/lineage-by-sub-admin-most-recent', f'pangolin_lineage={query}', collect_all=False)
    key_list = raw_data['results']
    key_list = list(key_list)
    
    for i in key_list:
        if i == key_list[0]:
            df = pd.DataFrame(raw_data['results'][i])
        else:
            newdf = pd.DataFrame(raw_data['results'][i]) # append each df
            df = pd.concat([df, newdf], sort=False)
    return df
    

def collection_date(pango_lin, mutations=None, location=None, server=server):
    """Most recent collection date by location

    Arguments:
     :param pango_lin: A string. (Required).
     :param mutations: (Somewhat optional). A string or list of strings. Comma separated list of mutations.
     :param location: (Somewhat optional). If not specified, return most recent date globally.
     :return: A pandas dataframe.
    """
    if mutations:
        if isinstance(mutations, list):
            mutations = ','.join(mutations)
        elif isinstance(mutations, str):
             mutations = mutations.replace(" ", "")
    
    query = pango_lin
    if mutations:
        query = query + '&' + f'mutations={mutations}'
    if location:
        query = query + '&' + f'location_id={location}'
        
    raw_data = get_outbreak_data('genomics/most-recent-collection-date-by-location', f'pangolin_lineage={query}', collect_all=False)
   
    data = {'Values' : raw_data['results']}
    df = pd.DataFrame(data) 
    return df


def submission_date(pango_lin, mutations=None, location=None, server=server):
    """Returns the most recent submission date by location

     Arguments:
     :param pango_lin: A string. (Required).
     :param mutations: (Somewhat optional). A string or list of strings. Comma separated list of mutations.
     :param location: (Somewhat optional). If not specified, return most recent date globally.
     :return: A pandas dataframe."""
    if mutations:
         if isinstance(mutations, list):
             mutations = ','.join(mutations)
         elif isinstance(mutations, str):
              mutations = mutations.replace(" ", "")
     
    query = pango_lin
    if mutations:
         query = query + '&' + f'mutations={mutations}'
    if location:
         query = query + '&' + f'location_id={location}'
         
    raw_data = get_outbreak_data('genomics/most-recent-submission-date-by-location', f'pangolin_lineage={query}', collect_all=False)
    
    data = {'Values' : raw_data['results']}
    df = pd.DataFrame(data) 
    return df
 
   
def mutation_details(mutations, server=server):
    """ Returns details of a mutation.
    
    Arguments:
     :param mutations: (Required). Comma separated list of mutations.
     :return: A pandas dataframe.
    """
    
    if isinstance(mutations, str):
         mutations = mutations.replace(" ", "")
    elif isinstance(mutations, list):
         mutations = ','.join(mutations)
   
    raw_data = get_outbreak_data('genomics/mutation-details', f'mutations={mutations}', collect_all=False)
    
    r = raw_data['results']
    keys = list(r[0])
   
    for i in r: # for each seperate result
        values = list(i.values())
        if i == r[0]:
            df=pd.DataFrame({"Key": keys,
                 "Values":values})
        else:
                newdf = pd.DataFrame({"Key": keys,
                     "Values":values}) # append each df
                df = pd.concat([df, newdf], axis=1, sort=False)
    return df


def daily_lag(location=None, server=server):
    """Return the daily lag between collection and submission dates by location

    Arguments:
     :param location_id: (Somewhat optional). If not specified, return lag globally.
     :return: A pandas dataframe.
    """
    query = ''
    if location:
        query =  '&' + f'location_id={location}'
        
    raw_data = get_outbreak_data('genomics/collection-submission', query, collect_all=False)
    
    r = raw_data['results']
    
    for i in r: # for each seperate result
        values = tuple(i.values())
        
        if i == r[0]:
            df=pd.DataFrame({"date_collected": values[0], "date_submitted": values[1], "total_count": values[2]})
        else:
                newdf = pd.DataFrame({"date_collected": values[0], "date_submitted": values[1], "total_count": values[2]}) # append each df
                df = pd.concat([df, newdf], sort=False)
    return df
    

def wildcard_lineage(name, server=server):
    """Match lineage name using wildcards. 

    Arguments:
    :param name: (Required). A string. Supports wildcards. Must use '*' at end of string. (Example: b.1*, ba.2*)
    :return: A pandas dataframe."""
    
    query = '' + '&' + f'name={name}'
    raw_data = get_outbreak_data('genomics/lineage', query, collect_all=False)
    r = raw_data['results']
    
    for i in r: # for each separate result
        values = tuple(i.values())
        if i == r[0]: # follow new procedure as found for daily_lag
            df=pd.DataFrame({"name": values[0],
                  "total_count":values[1]}, index=[0])
        else:
                newdf = pd.DataFrame({"name": values[0],
                      "total_count":values[1]}, index=[0]) # append each df
                df = pd.concat([df, newdf], sort=False)
    return df
     


def wildcard_location(name, server=server):
    """Match location name using wildcards. 

    Arguments:
    :param name: (Required). A string. Must use * at end of string. Supports wildcards. (Example: united*)
    :return: A pandas dataframe."""
    
    query = '' + '&' + f'name={name}'
    raw_data = get_outbreak_data('genomics/location', query, collect_all=False)
    r = raw_data['results']
   
    for i in r: # for each seperate result
        values = tuple(i.values())
        if i == r[0]:
            df=pd.DataFrame({"country": values[0], "country_id ": values[1],'id':values[2], "label":values[3],
                             "admin_level":values[4], "total_count":values[5]}, index = [0])
        else:
                newdf = pd.DataFrame({"country": values[0], "country_id ": values[1],'id':values[2], "label":values[3], 
                                      "admin_level":values[4], "total_count":values[5]}, index = [0]) # append each df
                df = pd.concat([df, newdf], sort=False)
    return df
     

def location_details(location, server=server):
    """Get location details using location ID.
     
    Arguments:
    :param location: A string. (Required).
    :return: Some pandas dataframes."""
   
    query = '' + '&' + f'id={location}'
    raw_data = get_outbreak_data('genomics/location-lookup', query, collect_all=False)
    data = {'Values' : raw_data['results']}
    df = pd.DataFrame(data) 
    return df

    
def wildcard_mutations(name, server=server):
    """Match mutations using wildcards.
    
     Arguments:
     :param name: (Required)  A string. Must use * at end of string. Supports wildcards. (Example: s:e484*)
     :return: A pandas dataframe."""

    query = '' + '&' + f'name={name}'
    raw_data = get_outbreak_data('genomics/mutations', query, collect_all=False)
    r = raw_data['results']
    
    for i in r: # for each seperate result
        values = tuple(i.values())
        if i == r[0]:
            df=pd.DataFrame({"name": values[0],
                  "total_count":values[1]}, index=[0])
        else:
                newdf = pd.DataFrame({"name": values[0],
                      "total_count":values[1]}, index=[0]) # append each df
                df = pd.concat([df, newdf], sort=False)
    return df

### Significance API enpoints: ###
    
def growth_rates(lineage, location='Global'):
    """Returns the growth rate score for a given lineage in a given location.
    
     Arguments:
     :param lineage: (Required)  A string. 
     :param location: (Required. Default: 'Global') A list or string. Separate multiple locations with ","
     :return: A pandas dataframe."""
    
    if isinstance(location, str):
        locations = location.replace(", " , "+OR+")
    elif isinstance(location, list):
             locations = '+OR+'.join(location)
    
    query = f'q=lineage:{lineage}+AND+location:{locations}'
    raw_data = get_outbreak_data('growth_rate/query', query, collect_all=False)
    df = pd.DataFrame(raw_data['hits'])
    
    return df

metadata_endpoint = "wastewater_metadata/query"
demix_endpoint = "wastewater_demix/query"
variants_endpoint = "wastewater_variants/query"

def get_ww_query(
    country=None,
    region=None,
    collection_site_id=None,
    date_range=None,
    sra_ids=None,
    viral_load_at_least=None,
    population_at_least=None,
    demix_success=True
):
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
      query_params.append(f"demix_success:{'true' if demix_success else 'false'}")

    return " AND ".join(query_params)

def get_wastewater_samples(**kwargs):
    """
    Retrieve wastewater samples data based on specified filters.

    Arguments:
    :param country: (Optional) Country name.
    :param region: (Optional) Region name.
    :param collection_site_id: (Optional) Site ID.
    :param date_range: (Optional) Date range in the format [start_date, end_date].
    :param sra_ids: (Optional) List of SRA IDs.
    :param viral_load_at_least: (Optional) Minimum viral load threshold.
    :param population_at_least: (Optional) Minimum population threshold.

    Returns:
    :return: A pandas DataFrame containing wastewater samples data.
    """
    query = get_ww_query(**kwargs)
    data = get_outbreak_data(metadata_endpoint, f"{dopage}&q=" + query, collect_all=True)
    try:
        df = pd.DataFrame(data['hits']).drop(columns=['_score', '_id'])
    except:
        raise KeyError("No data for query was found.")
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    df['ww_population'] = df['ww_population'].fillna(1000)
    return df

def get_wastewater_latest(**kwargs):
    query = get_ww_query(**kwargs)
    data = get_outbreak_data(metadata_endpoint, "size=1&sort=-collection_date&fields=collection_date&q=" + query) 
    try:
        return pd.DataFrame(data['hits'])['collection_date'][0]
    except:
        raise KeyError("No data for query was found.")

def fetch_ww_data(sample_metadata, endpoint):
    """
    Retrieve and join variants or demix info with sample metadata

    Arguments:
    :param sample_metadata: DataFrame containing metadata.
    :param endpoint: API endpoint to retrieve data.

    Returns:
    :return: A pandas DataFrame containing merged data with exploded nested data.
    """
    data = {"q": sample_metadata['sra_accession'].tolist(), "scopes": "sra_accession"}
    url = f'https://{server}/{endpoint}'
    response = requests.post(url, headers=get_auth(), json=data)
    df = pd.DataFrame(response.json()).drop(columns=['_score', '_id'])
    # Merge data with metadata
    merged_data = pd.merge(sample_metadata, df, on='sra_accession')
    # Explode nested data (old API structure)
    # exploded_data = pd.json_normalize(json.loads(merged_data.explode('variants' if 'variants' in data_df.columns else 'lineages').to_json(orient="records")))
    return merged_data.drop(columns='notfound', errors='ignore')

def get_wastewater_metadata(sample_metadata):
    """
    Add wastewater mutations data to a DF of samples

    Arguments:
    :param sample_metadata: DataFrame containing metadata.

    Returns:
    :return: A pandas DataFrame containing merged wastewater mutations data with metadata.
    """
    df = fetch_ww_data(sample_metadata, metadata_endpoint)
    df['viral_load'] = df['viral_load'].where(df['viral_load'] != -1, pd.NA)
    df['ww_population'] = df['ww_population'].fillna(1000)
    return df


def get_wastewater_mutations(sample_metadata):
    """
    Add wastewater mutations data to a DF of samples

    Arguments:
    :param sample_metadata: DataFrame containing metadata.

    Returns:
    :return: A pandas DataFrame containing merged wastewater mutations data with metadata.
    """
    return fetch_ww_data(sample_metadata, variants_endpoint)

def get_wastewater_lineages(sample_metadata):
    """
    Add wastewater demix results to a DF of samples

    Arguments:
    :param sample_metadata: DataFrame containing metadata.

    Returns:
    :return: A pandas DataFrame containing merged wastewater lineage abundance data with metadata.
    """
    return fetch_ww_data(sample_metadata, demix_endpoint)

def get_wastewater_samples_by_lineage(lineage):
    """
    Gets IDs of wastewater samples containing a certain lineage

    Arguments:
    :param lineage: Target lineage

    Returns:
    :return: A pandas series containing IDs of samples found to contain that lineage
    """
    data = get_outbreak_data(demix_endpoint, f"{dopage}&fields=sra_accession&q=name:{lineage}", collect_all=True)
    try:
        return pd.DataFrame(data['hits']).drop(columns=['_score', '_id'])
    except:
        raise KeyError("No data for query was found.")
    return data['sra_accession'].unique()

def get_wastewater_samples_by_mutation(site, alt_base=None):
    """
    Gets IDs of wastewater samples containing a mutation a certain site

    Arguments:
    :param site: Base pair index of mutations of interest
    :param alt_base: The new base at that site

    Returns:
    :return: A pandas series containing IDs of samples found to contain that lineage
    """
    alt_base = '' if alt_base is None else ' AND alt_base:'+alt_base
    data = get_outbreak_data(variants_endpoint, f"{dopage}&fields=sra_accession&q=site:{site}{alt_base}", collect_all=True)
    try:
        return pd.DataFrame(data['hits']).drop(columns=['_score', '_id'])
    except:
        raise KeyError("No data for query was found.")
    return data['sra_accession'].unique()

def normalize_viral_loads_by_site(df):
    site_vars = df.groupby('collection_site_id', observed=True)['viral_load'].std().rename('site_var')
    df = df.merge(site_vars, how='left', on='collection_site_id')
    df['normed_viral_load'] = df['viral_load'] / df['site_var']
    df['normed_viral_load'].where(~np.isfinite(df['normed_viral_load']), pd.NA)
    return df.drop(columns=['site_var'])

def datebin_and_agg(df, freq='7D', startdate=None, enddate=None, loaded=True):
    df = df.copy()
    if startdate is None: startdate = df['collection_date'].min()
    if enddate is None: enddate = df['collection_date'].max()
    bins = pd.date_range(startdate, enddate, freq=freq)
    df['date_bin'] = pd.cut(df['collection_date'], bins)
    df = df[~df['date_bin'].isna()]
    df['weight'] = df['normed_viral_load'] * df['ww_population'] if loaded else df['ww_population']
    agg_loads = lambda x: (x['normed_viral_load'] * x['ww_population']).sum() / x['ww_population'].sum()
    agg_abundance = lambda lin: lambda x: (x['abundance'] * (x['name'] == lin) * x['weight']).sum() / (x['abundance'] * x['weight']).sum()
    bins = df.groupby('date_bin', observed=True)
    agged_loads = bins.apply(agg_loads).rename('viral_load')
    agged_loads[agged_loads == float('inf')] = pd.NA
    agged_loads = agged_loads / agged_loads.max()
    agged_abundances = [ bins.apply(agg_abundance(lin)).rename(lin).fillna(0) for lin in df['name'].unique() ]
    agged_abundances = pd.DataFrame(agged_abundances).T
    agged_abundances = agged_abundances.rename(columns = {c:c.split('-like')[0] for c in agged_abundances.columns})
    agged_abundances = agged_abundances.T.groupby(agged_abundances.columns).sum().T
    agged_abundances = [agged_abundances[column] for column in agged_abundances.columns]
    return pd.concat( [agged_loads] + agged_abundances, axis=1)

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

def cluster_df(tree, clusters, df):
    (tree, lineage_key) = tree
    (U,V) = clusters
    if 'viral_load' in df.columns:
        abundances = df.drop(columns=['viral_load']).mul(df['viral_load'], axis=0)
    abundances = np.array(abundances.reindex(lineage_key.keys()).fillna(0))
    abundances = df.sum(axis=0)
    def get_agg_abundance(lin, abundances, W=set([])):
        cs = [get_agg_abundance(c, abundances, W) for c in lin['children'] if not c in W]
        return (abundances[lin['name']] if lin['name'] in abundances else 0) + np.sum(cs)
    viral_load = df['viral_load'] if 'viral_load' in df.columns else None
    df = df.drop(columns=['viral_load'])
    abundances_dated = [row for date,row in df.iterrows()]
    dates = [date for date,row in df.iterrows()]
    order = np.argsort([w['alias'] for w in list(U)+list(V)])
    lins = list(np.array(list(U)+list(V))[order])
    ulabels = [f'      {u["alias"]}.*' + (f' ({u["name"]})' if u["name"] != u["alias"] else '') for u in U]
    vlabels = [f'other {v["alias"]}.*' + (f' ({v["name"]})' if v["name"] != v["alias"] else '') for v in V]
    legend = list(np.array(ulabels+vlabels)[order])
    clustered_abundances = pd.DataFrame(
        { d: { label:get_agg_abundance(lin, a, U|V)
            for label, lin in zip(legend, lins) }
        for d,a in zip(dates, abundances_dated) } ).transpose()
    clustered_abundances[np.sum(clustered_abundances, axis=1) < 0.5] = pd.NA
    clustered_abundances['other *.*'] += 1 - clustered_abundances.sum(axis=1)
    if viral_load is not None: clustered_abundances = clustered_abundances.join(viral_load)
    return clustered_abundances, [lin['name'] for lin in lins], np.array([1]*len(U)+[0]*len(V))[order]
    
def get_colors(lins, isnatural, lineage_key):
    colors = np.searchsorted(
        sorted([lin['alias'] for lin in lineage_key.values()]),
        [lineage_key[lin]['alias'] for lin in lins] )
    colors = colors ** 2
    colors = (colors - np.min(colors)) / (np.max(colors)-np.min(colors)) * 0.75
    return hsv_to_rgb([(color, 1, 0.333 + 0.333*b) for color, b in zip(colors, isnatural)])

def get_riverplot_baseline(clustered_abundances, n=128):
    c = clustered_abundances.drop(columns=['viral_load']) \
         .mul(clustered_abundances['viral_load'].interpolate(), axis=0).dropna()
    d = c.div(clustered_abundances['viral_load'].dropna(), axis=0)
    def score(O):
        return (c.cumsum(axis=1).add(O, axis=0).rolling(window=2).apply(np.diff).mul(d)**2).sum().sum()
    Ot = -clustered_abundances['viral_load'].dropna()/2
    for n in range(128):
        O = np.random.normal(size=Ot.shape) / (n+48) * 2
        if score(O+Ot) < score(Ot):
            Ot += O
            Ot -= np.mean(Ot)
    return pd.Series(Ot, c.index).reindex(clustered_abundances.index).interpolate()