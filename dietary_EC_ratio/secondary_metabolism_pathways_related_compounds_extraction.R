# the aim of the script is to extract the secondary metabolism pathways related compounds in the KEGG database
# author: lu.zhang
# date: 2023.11.15 

# this script run in R

###############################################################################################################

#                                         Loading packages                                                    #

###############################################################################################################

library(tidyverse)

###############################################################################################################

#                                         prepare input file                                                  #

###############################################################################################################

# input dir 
input_dir <- "./extract_secondary_metablites/"

# nutrichem database 
nutrichem_database <- read_tsv(str_c(input_dir, "compound_associations.tsv"))

# combine metadata 
combine_metadata <- read_tsv(str_c(input_dir, "combine_metadata_sel.tsv"), guess_max = 3000)

# database file : before remove natural product. 
database_file <- readRDS(str_c(input_dir, "all_database_enzyme_compounds.rds"))

# modified ec profile : long format 
modified_ec_profile <- readRDS(str_c(input_dir, "modified_ec_profile_after_np.rds"))
modified_ec_profile_simple <- modified_ec_profile %>% 
	select(ec_simplied, Compounds) %>%
  	unique()

# input kegg pathway file 
input_pathway <- read_tsv(str_c(input_dir, "input_manual.txt"), col_names = FALSE)

# output_dir 
output_dir <- "./extract_secondary_metablites/result/"


###############################################################################################################

#                             extract all ecs in the pathways - secondary metabolites                         #

###############################################################################################################

library(KEGGREST)

#I use the link api function : reference :https://www.kegg.jp/kegg/docs/keggapi.html
#eg. https://rest.kegg.jp/link/ec/map00010 
#eg. https://rest.kegg.jp/link/cpd/map00010

input_pathway$map <- substr(input_pathway$X1, 1, 5)  

links_list <- list()
links_compound_list <- list()
for (map_id in input_pathway$map){
  
  
  links_list[[map_id]] <- keggLink("ec", paste0("map", map_id)) %>% data.frame(path = names(.), ec = .)
  links_compound_list[[map_id]] <- keggLink("compound", paste0("map", map_id)) %>% data.frame(path = names(.), compound = .) 
  
}

links_table <- do.call(rbind, links_list)
links_table$ec_simplied <- str_remove(links_table$ec, "ec:")

modified_ec_tibble <- tibble(ec = unique(modified_ec_profile$ec_simplied)) 
modified_ec_tibble_secondary <- modified_ec_tibble %>% filter(ec %in% links_table$ec_simplied) 


links_comp_table <- do.call(rbind, links_compound_list)

###############################################################################################################

#                             ec and compounds specific links - secondary metabolism related compounds        #

###############################################################################################################

# for each map, you extract the reactions that are in this map, 
map_rn <- keggLink("rn", "map00010") %>% data.frame(reaction = names(.), pathway = .) 
  
# then for each reaction, you extract all its related compounds and its related ecs 
rn_comp <- keggLink("compound", "rn:R00014") %>% data.frame(reaction = names(.), compound = .) 
rn_ec <- keggLink("ec", "rn:R00014") %>% data.frame(reaction = names(.), ec = .) 

# then you need to nest the compounds and ec links 
rn_ec_comp <- rn_comp %>% full_join(., rn_ec, by = "reaction")

# loop to extract all information 
## for some pathways, they don't have ec information recorded. 
nrow(input_pathway) #51 
links_table$path %>% unique() %>% length() #48
## not included pathways 
#"path:map00403" "path:map01052" "path:map01054"

links_comp_table$path %>% unique() %>% length() #51
# compounds part included all the pathways 

rn_ec_comp_table <- list()
for (each_path in links_table$path %>% unique()){
  print(each_path)
  input_path <- each_path %>% str_remove("path:")
  map_rn <- keggLink("rn", input_path) %>% data.frame(pathway = names(.), reaction = .) 
  rn_ec_comp <- list()
  for (each_reaction in unique(map_rn$reaction)){
    print(each_reaction)
    rn_comp <- keggLink("compound", each_reaction) %>% data.frame(reaction = names(.), compound = .) 
    rn_ec <- keggLink("ec", each_reaction) %>% data.frame(reaction = names(.), ec = .) 
    if (nrow(rn_ec) == 0){
      rn_ec <- data.frame(reaction = each_reaction, ec = "no") 
    }
    rn_ec_comp[[each_reaction]] <- rn_comp %>% full_join(., rn_ec, by = "reaction")
  }
  rn_ec_comp_table[[each_path]] <- do.call(rbind, rn_ec_comp)
}

# i loop for compounds unique pathway 

compound_only_path <- setdiff( links_comp_table$path, links_table$path)
for (each_path in compound_only_path){
  print(each_path)
  input_path <- each_path %>% str_remove("path:")
  map_rn <- keggLink("rn", input_path) %>% data.frame(pathway = names(.), reaction = .) 
  rn_ec_comp <- list()
  for (each_reaction in unique(map_rn$reaction)){
    print(each_reaction)
    rn_comp <- keggLink("compound", each_reaction) %>% data.frame(reaction = names(.), compound = .) 
    rn_ec <- keggLink("ec", each_reaction) %>% data.frame(reaction = names(.), ec = .) 
    if (nrow(rn_ec) == 0){
      rn_ec <- data.frame(reaction = each_reaction, ec = "no") 
    }
    rn_ec_comp[[each_reaction]] <- rn_comp %>% full_join(., rn_ec, by = "reaction")
  }
  
  
  rn_ec_comp_table[[each_path]] <- do.call(rbind, rn_ec_comp)
}


rn_ec_path_table <- do.call(rbind, rn_ec_comp_table)

rn_ec_path_table$pathway <- rownames(rn_ec_path_table)

rn_ec_path_table$modified_pathway <- substr(rn_ec_path_table$pathway, 1, 13)

# some pathways don't have a reaction matched: the compound ones have a tag : compound  
rn_ec_path_table_combine <- rn_ec_path_table %>% 
  full_join(., links_comp_table %>% mutate(tag = "compound"), by = c("modified_pathway" = "path", "compound" = "compound")) 

# combine with the ecs and pathways link , these pure ecs have a tag : ec 
rn_ec_path_table_combine_sel <- rn_ec_path_table_combine %>% filter(!is.na(tag)) %>%
  full_join(., links_table %>% mutate(tag_ec = "ec"), by = c("modified_pathway" = "path", "ec" = "ec")) #5888

# 
rn_ec_path_table_combine_sel$tag_ec[is.na(rn_ec_path_table_combine_sel$tag_ec)] <- "no"

rn_ec_path_table_combine_sel2 <- rn_ec_path_table_combine_sel %>% 
  filter(tag_ec!="no"|ec == "no"|(is.na(reaction) &tag == "compound")) # only keep the ones has ec recorded, or no ec in the pathway/reaction
# there are 4369 pairs of both ec-compound link 
# the logic here : since not all reactions related ecs are contained in the pathway, i first filter the ones which has ec recorded. or if there is no ec recorded in the reactions, then i extract all compounds. Or if there is no reaction, i also extract all compounds. 

###############################################################################################################

#                             convert the kegg compound ID into pubchem and chebis                            #

###############################################################################################################

sid <- keggConv("pubchem", links_comp_table$compound) %>% data.frame(path = names(.), compound = .)

saveRDS(sid, str_c(output_dir, "/secondary_sid.rds"))

# this part of analysis is run on the local linux computer 
cid_list <- list()
for (i in 1:length(sid$compound)){
  print(i)
  query <- sid$compound[i] %>% str_remove_all("pubchem:")
  cid_list[[i]] <- webchem::get_cid(query, from = "sid", domain = "substance")
  
}

chebi <- keggConv("chebi", links_comp_table$compound) %>% data.frame(path = names(.), chebi = .)

# the above small part of get cid information is obatained on local computer. 
cid_list <- readRDS(str_c(output_dir, "/secondary_cids_sid.rds"))
cid_list_table <- do.call(rbind, cid_list)

# match the cid to sid table then match the cid to links_comp_table 

sid_with_cid <- sid %>% mutate(compound_simple = str_remove(compound, "pubchem:")) %>% 
  left_join(., cid_list_table, by = c("compound_simple" = "query")) %>% 
  unique() %>% 
  mutate(kegg_comp = path) %>% 
  select(-path)
  
# match the chebi to links_comp_table 

chebi_modified <- chebi %>% mutate(kegg_comp = path) %>% select(-path) %>% unique()
  
links_comp_table_overall <- sid_with_cid %>% 
  select(cid, kegg_comp) %>% 
  left_join(links_comp_table, ., by = c("compound" = "kegg_comp")) %>% 
  left_join(., chebi_modified, by = c("compound" = "kegg_comp"))


######################################################################

#  total compounds in nutrichem 		   		     # 

######################################################################

# compounds contained in the pathways 
links_comp_table_overall$compound %>% unique() %>% length() #2423 

nutrichem_database[1,]$CompoundID %>% str_extract(., paste("CID:", 6453034, sep = ""))
nutrichem_database_withcid <- nutrichem_database %>% select(CompoundID) %>% unique() %>% mutate(cid = "na") 

for (i in 1:nrow(nutrichem_database_withcid)){
  print(i)
  nutrichem_database_withcid$cid[i] <- nutrichem_database_withcid$CompoundID[i] %>% str_extract_all(., "CID:\\d+") %>% unlist() %>% paste(., collapse = ",")
}

nutrichem_database_withcid_sep <- nutrichem_database_withcid %>% 
  separate_rows(., cid, sep = ",") %>% 
  left_join(., links_comp_table_overall %>% mutate(cid = paste("CID", cid, sep = ":")) %>% select(path, compound, cid))


# combine cid and chebi together 
nutrichem_database_withchebi <- nutrichem_database %>% select(CompoundID) %>% unique() %>% mutate(chebi = "na") 

for (i in 1:nrow(nutrichem_database_withchebi)){
  print(i)
  nutrichem_database_withchebi$chebi[i] <- nutrichem_database_withchebi$CompoundID[i] %>% str_extract_all(., "CHEBI:\\d+") %>% unlist() %>% paste(., collapse = ",")
}

nutrichem_database_withchebi_sep <- nutrichem_database_withchebi %>% 
  separate_rows(., chebi, sep = ",") %>% 
  left_join(., links_comp_table_overall %>% mutate(chebi = toupper(chebi)) %>% select(path, compound, chebi))

nutrichemid_chebi <- nutrichem_database_withchebi_sep %>% filter(chebi != "")  %>% .$CompoundID %>% unique()
nutrichemid_cid <- nutrichem_database_withcid_sep %>% filter(cid != "")  %>% .$CompoundID %>% unique()

union(nutrichemid_cid, nutrichemid_chebi) %>% length()

nutrichem_database_withcid_sep_sel <- nutrichem_database_withcid_sep %>% 
  filter(!is.na(path)) 


nutrichem_database_withchebi_sep_sel <- nutrichem_database_withchebi_sep %>% 
  filter(!is.na(path)) 

nutrichem_mapped_compounds <- bind_rows(nutrichem_database_withcid_sep_sel %>% select(CompoundID),
          nutrichem_database_withchebi_sep_sel %>% select(CompoundID)) %>% 
  unique() 


write_tsv(nutrichem_mapped_compounds, str_c(output_dir, "nutrichem_mapped_compounds.tsv"))
save.image(str_c(output_dir, "secondary_metabolites_extraction.RData"))

######################################################################

#  compounds modified by gut microbiota                              # 

######################################################################

# filter the modified ec by compounds 
modified_ec_profile_simple_sel <- modified_ec_profile_simple %>% filter(Compounds %in% nutrichem_mapped_compounds$CompoundID)
#186 compounds 


save.image(str_c(output_dir, "secondary_metabolites_extraction.RData"))


