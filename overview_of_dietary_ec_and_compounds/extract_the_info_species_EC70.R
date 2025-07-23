# extract the species and EC information 
library(tidyverse)
setwd("./minisp_70ECs")
load("overlap_and_unique_compounds_secondary_only_comp.RData")

secondary_compounds_overlap <- intersect(probiotic_compounds_ec_sep_comp$Compounds, modified_ec_profile_combine_simple_sel_compound$Compounds)
secondary_compounds_unique <- setdiff(modified_ec_profile_combine_simple_sel_compound$Compounds, probiotic_compounds_ec_sep_comp$Compounds)
# 116 overlapped compounds and 70 unique compounds 

# ec raw profile 
species_ec_link <- read_tsv("merge_Cohort.tsv")
colnames(species_ec_link)[1] <- "ec"
colnames(species_ec_link) <- str_remove_all(colnames(species_ec_link), "_Abundance-CPM") %>% str_remove_all(., "_T1") %>% str_remove_all("_Abundance-RPKs")

# combine metatable 
combine_metatable <- read_tsv("combine_metadata_sel.tsv", guess_max = 3000)

# metaphlan table - species table 
metaphlan_file <- read_tsv("merge_Cohort_metaphlan_species.tsv")
colnames(metaphlan_file) <- str_remove_all(colnames(metaphlan_file), "_metaphlan_bowtie2") %>% str_remove_all(., "_T1")
metaphlan_file_selv2 <- metaphlan_file %>%
  dplyr::select(., c("clade_name", combine_metatable$internal_sample_id))

#################################################

# species with compound with no pre-filter 

#################################################

secondary_tibble <- tibble(compounds = secondary_compounds_unique, group = "unique") %>% 
  add_row(compounds = secondary_compounds_overlap, group = "overlap")

ec_compounds_link <- modified_ec_profile_combine_simple_sel_compound %>% 
  select(ec_simplied, Compounds) %>% 
  unique() %>% 
  filter(Compounds %in% secondary_tibble$compounds) #186 compounds 

species_ec_link_sel <- species_ec_link %>%
  select(., c("ec", combine_metatable$internal_sample_id)) %>% #3068samples 
  add_column(ec_simplied = .$ec %>% str_split(., "\\|", simplify = T) %>% .[,1]) 

species_ec_link_sel_combine_sel <- species_ec_link_sel %>% 
  left_join(., ec_compounds_link) %>% 
  select(ec, ec_simplied, Compounds) %>% 
  unique() %>% 
  filter(!is.na(Compounds)) 

# remove all zero rows 
keep_items <- species_ec_link_sel %>% 
  select(-ec_simplied) %>% 
  column_to_rownames(., var = "ec") %>% 
  mutate(rowsum = rowSums(.)) %>% 
  filter(rowsum > 0)

# filtered out those without any linked compounds 
species_ec_link_sel_combine_sel_compounds <- species_ec_link_sel_combine_sel %>% 
  filter(ec %in% rownames(keep_items)) 

# linked the compounds with species name 
species_ec_link_sel_combine_sel_compounds_sp <- species_ec_link_sel_combine_sel_compounds %>% 
  add_column(species = str_split(.$ec, pattern = "\\|", simplify = TRUE)[,2]) %>%
  filter(species != "")

# only species and compounds and ec link 
ec_sp_compounds <- species_ec_link_sel_combine_sel_compounds_sp %>% 
  select(ec_simplied, species, Compounds) %>% # unclassified
  unique() 

# further mapped it with the phylum 
sp_prevalence <- metaphlan_file_selv2 %>% 
  pivot_longer(!clade_name) %>% 
  filter(value > 0) %>% 
  group_by(clade_name) %>% 
  summarise(sample_n = n()) %>% 
  ungroup() %>% 
  mutate(prevalence_sp = sample_n/3068) %>% 
  mutate(species = paste0(str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,6], ".", 
                          str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,7]))

# matched the name 
sp_compounds_prevalence_sel <- sp_prevalence %>% 
  left_join(ec_sp_compounds, .) %>% 
  filter(!is.na(clade_name)) %>% # remove the unmatched ones 
  filter(grepl("k__Bacteria", clade_name)) 


#################################################

# compound-species heatmap

#################################################

# prevalence filtering 
sp_compounds_prevalence_sel_0.05 <- sp_compounds_prevalence_sel %>% 
  filter(prevalence_sp > 0.05)

sp_compounds_prevalence_sel_0.05_wide <- sp_compounds_prevalence_sel_0.05 %>% 
  select(species, Compounds) %>% 
  unique() %>% 
  mutate(num = 1) %>% 
  pivot_wider(names_from = "Compounds", values_from = num) %>% 
  column_to_rownames(., var = "species") 

sp_compounds_prevalence_sel_0.05_wide[is.na(sp_compounds_prevalence_sel_0.05_wide)] <- 0

rowannotation <- tibble(species = rownames(sp_compounds_prevalence_sel_0.05_wide)) %>% 
  left_join(., sp_compounds_prevalence_sel %>% select(species, clade_name) %>% unique()) %>% 
  mutate(phylum = str_split(clade_name, pattern = "\\|", simplify = T) %>% .[,2]) %>% 
  select(species, phylum) %>% 
  column_to_rownames(., var = "species")
columnannotateion <- secondary_tibble %>% 
  filter(compounds %in% colnames(sp_abundance_prevalence_sel_0.05_wide)) %>% 
  column_to_rownames(., var = "compounds")


saveRDS(sp_abundance_prevalence_sel_0.05_wide, "sp_abundance_prevalence_sel_0.05_wide.rds")
saveRDS(rowannotation, "rowannotation.rds")
saveRDS(columnannotateion, "columnannotateion.rds")



