import pandas as pd
import numpy as nPy
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as seaB
import matplotlib.pyplot as plt

# Refine canonical smiles data by droping it
df = pd.read_csv('/Users/hung/Desktop/BioInformatic-Project/acetylcholinesterase_03_bioactivity_data_curated.csv')
df_no_smiles = df.drop(columns='canonical_smiles')

# Use Lipinski method and measure to evaluting the druglikeness compound
smiles = []

for i in df.canonical_smiles.tolist():
  cpd = str(i).split('.')
  cpd_longest = max(cpd, key = len)
  smiles.append(cpd_longest)

smiles = pd.Series(smiles, name = 'canonical_smiles')

df_clean_smiles = pd.concat([df_no_smiles,smiles], axis=1)

def lipinski(smiles, verbose=False):

    arrayData= []
    for element in smiles:
        mole=Chem.MolFromSmiles(element)
        arrayData.append(mole)
       
    base= nPy.arange(1,1)
    i=0  
    for mole in arrayData:
       
        MolWT = Descriptors.MolWt(mole)
        MolLogP = Descriptors.MolLogP(mole)
        Donors = Lipinski.NumHDonors(mole)
        NumAccept = Lipinski.NumHAcceptors(mole)
           
        rowT = nPy.array([MolWT,
                        MolLogP,
                        Donors,
                        NumAccept])
    
        if(i==0):
            base=rowT
        else:
            base=nPy.vstack([base, rowT])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=base,columns=columnNames)
    
    return descriptors

df_lipinski = lipinski(df_clean_smiles.canonical_smiles)


df_combined = pd.concat([df,df_lipinski], axis=1)


def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molarVar = i*(10**-9) # Converts nM to M
        pIC50.append(-nPy.log10(molarVar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x


def norm_value(input):
    normal = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        normal.append(i)

    input['standard_value_norm'] = normal
    x = input.drop('standard_value', 1)
        
    return x

dataF_norm = norm_value(df_combined)

dataF_Final = pIC50(dataF_norm)

dataF_2C = dataF_Final[dataF_Final['class'] != 'intermediate']

seaB.set(style='ticks')

plt.figure(figsize=(5.5, 5.5))

seaB.countplot(x='class', data=dataF_2C, edgecolor='black')

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_bioactivity_class.pdf')


#Plot plC50 value

plt.figure(figsize=(5.5, 5.5))

seaB.boxplot(x = 'class', y = 'pIC50', data = dataF_2C)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_ic50.pdf')

#MW graph

plt.figure(figsize=(5.5, 5.5))

seaB.boxplot(x = 'class', y = 'MW', data = dataF_2C)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_MW.pdf')

#LogP

plt.figure(figsize=(5.5, 5.5))

seaB.boxplot(x = 'class', y = 'LogP', data = dataF_2C)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_LogP.pdf')

#Num H Donors

plt.figure(figsize=(5.5, 5.5))

seaB.boxplot(x = 'class', y = 'NumHDonors', data = dataF_2C)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_NumHDonors.pdf')

# Num H Acceptors

plt.figure(figsize=(5.5, 5.5))

seaB.boxplot(x = 'class', y = 'NumHAcceptors', data = dataF_2C)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

plt.savefig('/Users/hung/Desktop/BioInformatic-Project/plot_NumHAcceptors.pdf')

