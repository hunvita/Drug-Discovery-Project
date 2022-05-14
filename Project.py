
# Import libraries
from cmath import nan
import pandas as pd
import os
from chembl_webresource_client.new_client import new_client

# wantedTarget search
wantedTarget = new_client.wantedTarget
query = wantedTarget.search('acetylcholinesterase')
targets = pd.DataFrame.from_dict(query)


# Assign the fifth entry (which corresponds to the wantedTarget protein, Human Acetylcholinesterase) to the selected_target variable
selected_target = targets.target_chembl_id[0]


# Retrieve only bioactivity data for Human Acetylcholinesterase
active = new_client.active
resolution = active.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(resolution)

# Handling missing data
df2 = df[df.standard_value.notna()]
df2 = df2[df.canonical_smiles.notna()]


# Drop dublicate data
df2_nr = df2.drop_duplicates(['canonical_smiles'])

# Combine the 3 columns (molecule_chembl_id,canonical_smiles,standard_value) and bioactivity_class into a DataFrame
select = ['molecule_chembl_id','canonical_smiles','standard_value']
df3 = df2_nr[select] # df3 = acetylcholinesterase_02_bioactivity_data_preprocessed

os.makedirs('/Users/hung/Desktop/BioInformatic-Project', exist_ok=True)
df3.to_csv('/Users/hung/Desktop/BioInformatic-Project/acetylcholinesterase_02_bioactivity_data_preprocessed.csv', index=False)


df4 = pd.read_csv('/Users/hung/Desktop/BioInformatic-Project/acetylcholinesterase_02_bioactivity_data_preprocessed.csv')

# Categorize the strength of the compounds, less than 1000 mean strong and will be active, 
# greater than 10000 mean weak and will be inactive, in between consider intermidiate
bioactivity_threshold = []
for i in df4.standard_value:
  if float(i) >= 10000:
    bioactivity_threshold.append("inactive")
  elif float(i) <= 1000:
    bioactivity_threshold.append("active")
  else:
    bioactivity_threshold.append("intermediate")


bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df5 = pd.concat([df4, bioactivity_class], axis=1) # df5 = acetylcholinesterase_03_bioactivity_data_curated

# Output the curated data
os.makedirs('/Users/hung/Desktop/BioInformatic-Project', exist_ok=True)
df5.to_csv('/Users/hung/Desktop/BioInformatic-Project/acetylcholinesterase_03_bioactivity_data_curated.csv', index=False)



