# replace the manually corrected IDs 

library(tidyverse)
smiles_manual <- read_tsv("invalid_compounds_manual.tsv")

smiles_ori <- read_tsv("full_plant_compound_structures.smi", col_names = F)

for (i in smiles_manual$CompoundID){
  print(i)
  revised_smile <- smiles_manual$revised_smile[smiles_manual$CompoundID == i]
  smiles_ori$X2[smiles_ori$X1 == i] <- revised_smile
  
}

write_tsv(smiles_ori, file = "smiles_replace_invalid_ones.tsv", col_names = F)