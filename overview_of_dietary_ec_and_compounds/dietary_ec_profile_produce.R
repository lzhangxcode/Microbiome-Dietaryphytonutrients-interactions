# the aim of this script is to produce species accumulation curve of modified ECs 
# author: lu.zhang
# date: 16th,Oct,2023 
# the output will be the rarefied EC table as well as the modified EC profile


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Rscript script.R summary_database_rdata ec_profile_diff_continent manual_checked_list output_dir", call.=FALSE)
} 


############################################################

# loading package                                          #

############################################################

library(tidyverse)
library(gridExtra)



############################################################

# input data                                               #

############################################################

# RData of the summary database 
input_database_rdata <- args[1]

input_ec_profile <- args[2]

checked_np_db <- args[3]

output_dir <- args[4]


############################################################                             

# database file                                            #                             

############################################################ 

print("produce processed database file.")
load(input_database_rdata)

read_database_file <- function(database_rdata = NA, output_dir = NA){
    load(database_rdata)
    kegg_annotated$tag <- "KEGG"
    colnames(kegg_annotated)[2] <- "EC"
    pubchem_annotated$tag <- "pubchem"
    colnames(pubchem_annotated)[2] <- "EC"
    HMDB_CID_annotated$tag <- "HMDB"
    colnames(HMDB_CID_annotated)[2] <- "EC"
    Brenda_merge$tag <- "Brenda"
    colnames(Brenda_merge)[2] <- "EC"
    EZCAT_nutrichem_annotated$tag <- "EZCAT"
    colnames(EZCAT_nutrichem_annotated)[2] <- "EC"
    IntEnz_CHEBI_nutrichem_annotated$tag <- "IntEnz"
    colnames(IntEnz_CHEBI_nutrichem_annotated)[2] <- "EC"
    EAWAG_CID_nutrichem_annotated$tag <- "EAWAG"
    colnames(EAWAG_CID_nutrichem_annotated)[2] <- "EC"
    SFLD_CHEBI_nutrichem_annotated$tag <- "SFLD"
    colnames(SFLD_CHEBI_nutrichem_annotated)[2] <- "EC"
    transformer_cid_nutrichem_annotated$tag <- "transformer"
    colnames(transformer_cid_nutrichem_annotated)[2] <- "EC"
    M_CSA_nutrichem_annotated_update$tag <- "M_CSA"
    colnames(M_CSA_nutrichem_annotated_update)[2] <- "EC"


    all_database_file <- rbind(kegg_annotated,
                                           pubchem_annotated,
                                           HMDB_CID_annotated,
                                           Brenda_merge,
                                           EZCAT_nutrichem_annotated,
                                           IntEnz_CHEBI_nutrichem_annotated,
                                           EAWAG_CID_nutrichem_annotated,
                                           SFLD_CHEBI_nutrichem_annotated,
                                           transformer_cid_nutrichem_annotated,
             M_CSA_nutrichem_annotated_update
    )

    all_database_file_separate <- all_database_file %>% 
            separate_rows(., EC, sep = "/EC;") %>%  
            separate_rows(., EC, sep = ",") %>%  #HMDB
            separate_rows(., EC, sep = ";") %>% # excep HMDB, all other database have ; separation 
            separate_rows(., EC, sep = "AND") %>% # Pubchem used, to separate different reactions 
            filter(., EC != "") %>% 
            add_column(ec_simplied = str_remove(.$EC, "ec:"))

    write_delim(all_database_file_separate, paste(output_dir, "all_database_file_separate.txt", sep = "/"), delim = "\t")

    return(all_database_file_separate)

}

all_database_file_separate <- read_database_file(database_rdata = input_database_rdata,
                   output_dir = output_dir)


print(head(all_database_file_separate)) 


all_database_enzyme_compounds <- all_database_file_separate %>%
        select(ec_simplied, CompoundID, tag) %>%
        unique() %>%
        group_by(ec_simplied) %>%
        summarise(Compounds = paste(CompoundID, collapse = ","),
              Tag = paste(tag, collapse = ",")) #5130 enzymes in total 

saveRDS(all_database_enzyme_compounds, paste(output_dir, "all_database_enzyme_compounds.rds", sep = "/"))


############################################################

# combine the EC profiles with the database                #

############################################################

print("combine the ec profile with nutrichem database.")

# read files 
ec_profile <- list()

files <- list.files(input_ec_profile, pattern = "rds")

for (file in files){
    continent <- str_remove(file, "ec_unstratified_") %>% str_remove(".rds")
    print(continent)
    ec_profile[[continent]] <- readRDS(paste(input_ec_profile, file, sep = "/"))
}

ec_profile_sel <- list()
for (i in names(ec_profile)){
    print(i)
    ec_profile_sel[[i]] <- ec_profile[[i]] %>%
    left_join(., all_database_enzyme_compounds, by = c("Genefamilies" = "ec_simplied")) %>%
    filter(., !is.na(Compounds))

    saveRDS(ec_profile_sel[[i]], paste(output_dir, paste("modified_ec_profile_update_ramps_continent/ec_unstratified_modified_", i, ".rds", sep = ""), sep = "/"))
    write_tsv(ec_profile_sel[[i]], paste(output_dir, paste("modified_ec_profile_update_ramps_continent/ec_unstratified_modified_", i, ".tsv", sep = ""), sep = "/"))


}



############################################################

#	combine the ec profile together as one profile     # 

############################################################

combine_ec_profile <- function(input_ec_dir = NA, tag = "part"){
    modified_ec_profile <- list()

    files_modified <- list.files(input_ec_dir, pattern = "ec_unstratified.*rds")
    length <- 1 
    for (file_md in files_modified){
        continent <- str_remove(file_md, "ec_unstratified_") %>% str_remove(".rds")
        print(continent)
        if (tag == "all"){
            modified_ec_profile[[continent]] <- readRDS(paste(input_ec_dir, file_md, sep = "/"))
        }else if(tag == "part"){
            modified_ec_profile[[continent]] <- readRDS(paste(input_ec_dir, file_md, sep = "/")) %>% 
            select(-Compounds, -Tag)
        }
        
        if (length == 1){
            combine_ec_modified_profile <- modified_ec_profile[[continent]]
        }else if(length > 1){
            combine_ec_modified_profile <- full_join(combine_ec_modified_profile, modified_ec_profile[[continent]], by = "Genefamilies") 
        }
        length <- 2

        removed_na_dim <-  modified_ec_profile[[continent]] %>% drop_na() %>% 
            dim() %>% 
            as.vector()
        original_dim <- modified_ec_profile[[continent]] %>% dim() %>% 
            as.vector()
        if (identical(removed_na_dim, original_dim)){
            print(paste("no NA in original ec profile:",continent))
        }else{
            print(paste("wrong:", continent, sep = ""))
        }
    }
    
    combine_ec_modified_profile_remove_na <- combine_ec_modified_profile %>%
        replace(is.na(.), 0) %>%
        column_to_rownames(., var = "Genefamilies")

    combine_ec_modified_profile_remove_na <- combine_ec_modified_profile_remove_na[rowSums(abs(combine_ec_modified_profile_remove_na)) != 0,]

    return(combine_ec_modified_profile_remove_na)
}


modified_ec_profile <- combine_ec_profile(input_ec_dir = paste(output_dir, "modified_ec_profile_update_ramps_continent", sep = "/"))




############################################################

#  filter further manually natural product                 #

############################################################

modified_ec_profile_database <- modified_ec_profile %>%
    rownames_to_column(., var = "ec_simplied") %>% 
    left_join(., all_database_enzyme_compounds, by = "ec_simplied") %>%
    filter(., !is.na(Compounds))


compounds_num <- modified_ec_profile_database  %>% 
	separate_rows(Compounds, sep = ",") 

print(paste("this is the final compounds number after mapping to the gut microbiota:", length(unique(compounds_num$Compounds)), sep = " "))


checked_np_db_content <- readRDS(checked_np_db)

# filter the manually checked compounds 

modified_ec_profile_after_np <- compounds_num %>% 
	filter(Compounds %in% checked_np_db_content$CompoundID) 

# 775 compounds are left 
print(paste("this is the final compounds number after natural product database check:", length(unique(modified_ec_profile_after_np$Compounds)), sep = " "))

#####################################################################

#		save image 					    # 

##################################################################### 

print("save the rdata")
save.image(paste(output_dir, "step10_rarefaction_curve_modified_ec_profile_manual_check.RData",sep = "/"))



