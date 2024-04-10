import sys

from outbreak_data.authenticate_user import set_authentication

from outbreak_data import outbreak_data

set_authentication(sys.argv[1])

d = outbreak_data.get_wastewater_latest(region="Ohio")
s = outbreak_data.get_wastewater_samples(region="Ohio", date_range=["2023-06-01", d])
print(outbreak_data.get_wastewater_lineages(s))

s = outbreak_data.get_wastewater_samples(site_id="USA_OH_5f9e5487", viral_load_at_least=25000)
print(outbreak_data.get_wastewater_mutations(s))

s = outbreak_data.get_wastewater_samples_by_lineage('EG.5.1')
s = outbreak_data.get_wastewater_metadata(s)
print(outbreak_data.get_wastewater_mutations(s))

s = outbreak_data.get_wastewater_samples_by_mutation(site=1003, alt_base='G')
print(outbreak_data.get_wastewater_lineages(s))

