# the script is used to highlight the enriched group 
# my-rdkit-env
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
from PIL import Image, ImageOps


# read files 
df_smarts = pd.read_table('FragmentDescriptors_manual.csv')
df_compounds = pd.read_table('functional_group_result.tsv')

# filtered only targeted columns 
# Extract columns from df_compounds that are present in df_smarts
relevant_columns = ['CompoundID', 'smiles'] + [col for col in df_compounds.columns if col in df_smarts['Code'].values]
df_compounds_filtered = df_compounds[relevant_columns]
df_compounds_filtered.to_csv('filtered_compounds_thick.tsv', sep='\t', index_label='Index')

# build dictionary 
smarts_dict = {row['Code']: row['SMARTS'] for _, row in df_smarts.iterrows()}

opts = Draw.MolDrawOptions()
opts.bondLineWidth = 5

# read through each compounds and highlight 
# # Process each compound and each functional group
for index, row in df_compounds_filtered.iterrows():
    compound_id = index  # Use the index as the compound ID for the file name
    smiles = row['smiles']
    mol = Chem.MolFromSmiles(smiles)
    
    # Check each SMARTS pattern
    for code, smarts in smarts_dict.items():
        if code in df_compounds_filtered.columns and row[code] > 0:
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            highlight_atoms = []
            for match in matches:
                highlight_atoms.extend(list(match))
            
            # Generate image with highlighted atoms for each functional group
            if highlight_atoms:
                highlight_atoms = list(set(highlight_atoms))  # Remove duplicates
                img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms, size=(1000, 1000), options=opts)
                border_thickness = 0  # Adjust this value for thicker border
                bordered_img = ImageOps.expand(img, border=border_thickness, fill='white')
                bordered_img.save(f"images_thick/{compound_id}_{code}_highlighted.png")



