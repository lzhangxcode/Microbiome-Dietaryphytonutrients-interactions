# extract fragments from the revised smile file 
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments

# read file 
df_1 = pd.read_table('smiles_replace_invalid_ones.tsv', sep = "\t", header=None) #nrows = 20port pandas as pd
df_2 = pd.read_table('。/secondary_compounds_overlap.tsv')
df_1.columns =['CompoundID', 'smiles']

# filter the compounds in df_2 
targeted_compound_ids = df_2['compounds'].unique()

# Filter and merge the data with SMILES
filtered_data = df_1[df_1['CompoundID'].isin(targeted_compound_ids)] #186 rows 

# Get all fragment functions from rdkit.Chem.Fragments: 
fragment_functions = [func for func in dir(Fragments) if func.startswith('fr_')]

# Define a function to calculate all RDKit fragments
# some smiles are converted unsuccessfully, i need to find them and replace them with CID matched smiles 
def write_invalid_compounds_to_tsv(compound_id, smiles, filename):
    with open(filename, 'a') as f:
        f.write(f"{compound_id}\t{smiles}\n")



def calculate_rdkit_fragments(smiles, compound_id, invalid_compounds_file):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        fragment_counts = {}
        for func in fragment_functions:
            fragment_counts[func] = getattr(Fragments, func)(mol)
        return fragment_counts
    except Exception as e:
        print(f"Invalid SMILES string for compound {compound_id}: {smiles}")
        print(e)
        write_invalid_compounds_to_tsv(compound_id, smiles, invalid_compounds_file)
        return {}



# Calculate functional groups for each compound
invalid_compounds_file = 'invalid_compounds_second_times.tsv'
functional_groups_list = []
for _, row in filtered_data.iterrows(): # read every lines 
    smiles = row['smiles']
    compound_id = row['CompoundID']
    print(compound_id)
    functional_groups = calculate_rdkit_fragments(smiles,compound_id, invalid_compounds_file)
    functional_groups['CompoundID'] = compound_id
    functional_groups_list.append(functional_groups)

# Convert to DataFrame
functional_groups_df = pd.DataFrame(functional_groups_list)

# Merge with filtered_data
result = pd.merge(filtered_data, functional_groups_df, on='CompoundID')
result.to_csv('functional_group_result.tsv', sep='\t', index=False)



