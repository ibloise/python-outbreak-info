.. Py_Outbreak_API documentation master file, created by
   sphinx-quickstart on Sat Oct  8 22:16:44 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the Python Outbreak.info package docs!
========================================================
Here you can find information on the functions you will use to collect and analyze SARS-COV-2  data from the Outbreak.info API. 
Our package pulls data from the `Outbreak.info API  <https://api.outbreak.info/>`_ and is reflected on our `Outbreak.info web interface <https://outbreak.info/>`_

Installation
----------------
We recommend installing the package via pip using:

``pip install python-outbreak-info``

Alternatively, the package can be directly installed from source via pip:

``pip install git+https://github.com/outbreak-info/python-outbreak-info.git``

Getting Started
----------------
The Python Outbreak.info package contains key functions for accessing genomic and epidemiological data for SARS-CoV-2. Access to genomic data requires logging in using GISAID credentials to generate an API key, using the ``authenticate_new_user()`` function. To perform authentication, you'll need to first run

.. code-block:: python

    from outbreak_data.authenticate_user import authenticate_new_user
    authenticate_new_user()

and then you should be able access all of the functionality of the package. Most of the rest of the tools are available within the ``outbreak_data`` component of the package. For example: 

.. code-block:: python

   from outbreak_data import outbreak_data
   lin_list = ['B.1.1.7','B.1.351','B.1.617.2']
   # request lineages occurring with minimum frequency of 0.05 (5%)
   df = outbreak_data.known_mutations(lin_list,freq=0.05)
   # filter mutations and sort by codon number
   df = df[df['gene']=='S'].sort_values(by='codon_num')

For wastewater abundance analyses, users will need to supply the appropriate location code corresponding to their location of interest and a date range. To do this, users would first retrieve wastewater data from ``outbreak_data`` then aggregate the data by date and weight to get the abundances for each lineage using the ``outbreak_tools`` part of the package. An example lookup should look like: 

.. code-block:: python

   from outbreak_data import outbreak_data
   from outbreak_tools import outbreak_tools
   state = "Ohio"
   startdate, enddate = "2023-09-01", "2023-10-01"
   ww_samples = outbreak_data.get_wastewater_samples(region=state, date_range=[startdate, enddate])
   ww_samples = outbreak_data.get_wastewater_lineages(ww_samples)
   ww_abundances = outbreak_tools.datebin_and_agg(ww_samples, weights=outbreak_tools.get_ww_weights(ww_samples), startdate=startdate, enddate=enddate, freq='7D')
   # This data frame is large, so we'll focus on one lineage
   ww_abundances['B.1.1.191'] 

   (2023-08-31, 2023-09-07]    0.068126
   (2023-09-07, 2023-09-14]    0.017081
   (2023-09-14, 2023-09-21]    0.030141
   (2023-09-21, 2023-09-28]    0.031744
   Name: EG.2, dtype: float64
   

**About Clinical and Wastewater Tools**

Toward the beginning of the SARS-Cov-2 pandemic, viral genome sequencing data were collected through specimens that were obtained from clinical testing. Yet studies have shown that the method has introduced sampling bias due to systemic healthcare disparities, particularly in poor and underserved communities. 

In contrast, wastewater samples have been highly useful for tracking regional infection dynamics while providing less biased abundance estimates than clinical testing. Data collected by tracking viral genomic sequences in wastewater has also improved community prevalence estimates and detects emerging variants earlier on. 

The Andersen Lab has developed improved virus concentration protocols and deconvolution software that fully resolve multiple virus strains from wastewater. The resulting data is now deployed by Python-outbreak-info. In short, SARS-Cov-2 analysis can be done using both clinical and wastewater tools, yet data from the wastewater analysis tools may be more accurate in some situations.

`Click here: <https://www.nature.com/articles/s41586-022-05049-6>`_ for more information on wastewater analysis.


Table of Contents:
===================


Core Outbreak Data Tools
--------------------------
.. toctree::
   all_lineage_prevalences
   auth_setup
   cases_by_location
   daily_lag
   growth_rates
   gr_significance
   known_mutations
   lineage_by_sub_admin
   lineage_cl_prevalence
   most_recent_cl_data
   mutation_details
   mutation_prevalences
   seq_counts
   wildcard_location
   wildcard_lineage

Wastewater Analysis Tools
--------------------------
.. toctree::
   get_wastewater_latest
   get_wastewater_lineages
   get_wastewater_metadata
   get_wastewater_mutations
   get_wastewater_samples
   get_wastewater_samples_by_lineage
   get_wastewater_samples_by_mutation


Plotting and Organization Toolkit
----------------------------------
.. toctree::
   cluster_df
   const_idx
   datebin_and_agg
   first_date
   get_colors
   get_riverplot_baseline
   get_tree
   get_ww_weights


Example Applications and Analyses
----------------------------------
.. toctree::
   Epidemiological data analyses <epidem_analysis>
   Mutation Data Analyses <mutation_analysis>
   Dealing with Cryptic Variants <cryptic_vars> 
