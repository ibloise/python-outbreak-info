import requests
import pandas as pd

server = 'api.outbreak.info'  # or 'dev.outbreak.info'
auth = ***REMOVED***  # keep this private!
nopage = 'fetch_all=true&page=0'  # worth verifying that this works with newer ES versions as well
covid19_endpoint = 'covid19/query'


def get_outbreak_data(endpoint, argstring, server=server, auth=auth):
    """
    Recieves raw data using outbreak API
       
    Arguments: 
        endpoint: directory in server the data is stored
        argstring: feature arguments to provide to API call
    
    Returns: 
        A request object containing the raw data
       
    """
    auth = {'Authorization': str(auth)}
    return requests.get(f'https://{server}/{endpoint}?q={argstring}', headers=auth)


def page_data(data_origin, collect, num_pages = None, server=server, auth=auth, covid19_endpoint=covid19_endpoint):
    """
    Used as a helper function to page through results and return more data.
    :param data_origin: Initial request containing the args
    :param num_pages: Numer of pages to scroll.
    :param server: Used to call which server in the API.
    :param auth: Used as a header in call to the API.
    :param covid19_endpoint: Used to specify the endpoint the covid19 data is stored.
    :param collect: if True, returns max pages & disables num_pages.
    :return: Pandas Dataframe
    """
    if collect:
        assert(num_pages is None)
    else:
        assert(num_pages is not None)

    auth = {'Authorization': str(auth)}
    scroll_id = data_origin.json()['_scroll_id']
    scroll_df = pd.DataFrame(columns=pd.Series(data_origin.json()['hits'][0]).index)

    fetching_page = '&fetch_all=True&page='
    curr_page = 0
    if collect:
        while 'success' not in data_origin.json().keys():
            # individual request df
            data = pd.DataFrame(data_origin.json()['hits'])
            scroll_df = scroll_df.append(data, ignore_index=True)

            page = fetching_page + str(curr_page)
            to_scroll = 'scroll_id=' + scroll_id + page
            data_origin = requests.get(f'https://{server}/{covid19_endpoint}?{to_scroll}', headers=auth)
            curr_page += 1
        # applying datetime to dates column and sorting in ascending
        scroll_df['date'] = scroll_df['date'].apply(lambda x: pd.to_datetime(x))
        scroll_df = scroll_df.sort_values(by='date', ascending=True)
        scroll_df.reset_index(drop=True, inplace=True)
        return scroll_df
    else:
        while curr_page <= num_pages:
            # individual request df
            data = pd.DataFrame(data_origin.json()['hits'])
            scroll_df = scroll_df.append(data, ignore_index=True)

            page = fetching_page + str(curr_page)
            to_scroll = 'scroll_id=' + scroll_id + page
            data_origin = requests.get(f'https://{server}/{covid19_endpoint}?{to_scroll}', headers=auth)
            curr_page += 1
        # applying datetime to dates column and sorting in ascending
        scroll_df['date'] = scroll_df['date'].apply(lambda x: pd.to_datetime(x))
        scroll_df = scroll_df.sort_values(by='date', ascending=True)
        scroll_df.reset_index(drop=True, inplace=True)
        return scroll_df
        
def get_prevalence_by_location(endpoint, argstring, server=server, auth=auth):
    
    """Used with prevalence_by_location. Works similarly to get_outbreak_data. 
       Prevalence_by_location() requires a different url string.
    
    Arguments: 
        endpoint: directory in server the data is stored
        argstring: feature arguments to provide to API call
    
    Returns: 
        A request object containing the raw data"""
    auth = {'Authorization': str(auth)}
    return requests.get(f'https://{server}/{endpoint}?{argstring}', headers=auth) 



     



                              


    

