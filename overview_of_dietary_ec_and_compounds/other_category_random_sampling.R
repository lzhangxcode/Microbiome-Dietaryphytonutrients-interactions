# this script is aim to do random sampling for ec & species
# author: l.zhang 
# date: 29,jan,2024 - 30,jan,2024 

########################################

#       packages needed                # 

########################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
	  stop("Rscript other_category_random_sampling.R combine_metatable metaphlan dietary_ec input_ec_profile_dir output_dir", call.=FALSE)
} 


library(tidyverse)
library(ggplot2)

########################################


#       input                          # 


########################################

# combined metatable 
combine_metatable <- args[1] %>% read_tsv(., guess_max = 3000)

# species table 
metaphlan <- args[2]
metaphlan_file <- read_tsv(metaphlan)
colnames(metaphlan_file) <- str_remove_all(colnames(metaphlan_file), "_metaphlan_bowtie2") %>% str_remove_all(., "_T1")
metaphlan_file_sel <- metaphlan_file %>%
    select(., c("clade_name", combine_metatable$internal_sample_id))  

# modified ec profile 
modified_ec <- readRDS(args[3])

# total ec profile 
input_ec_profile_dir <- args[4]

# output dir 
output_dir <- args[5]
#output_dir <- "othercategory/"


########################################

# combine ec profile                   # 

########################################

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


        # check all ec profile don't have NA in it 
        # all of the data doesn't have na values before join.
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


all_ec_profile <- combine_ec_profile(input_ec_dir = input_ec_profile_dir, tag = "all")



########################################


#  select only valid samples           # 
#  and obtain the unstratified table   #

########################################


all_ec_profile_clean <- all_ec_profile %>% 
        rownames_to_column(., var = "ec") %>%
        select(c("ec", combine_metatable$internal_sample_id))



separate_samples_into_continents <- function(input_profile = NA, continent_meta = NA){

        output_profile <- list()

        for (i in unique(continent_meta$continent)){
        
                print(i)
                target_samples <- continent_meta %>% 
                        filter(continent == i)
                id <- colnames(input_profile)[1]
                output_profile[[i]] <- input_profile %>%
                        select(colnames(input_profile)[1], target_samples$internal_sample_id)
        }
                 
        return(output_profile)

}



all_ec_sep <- list()
all_ec_sep <- separate_samples_into_continents(input_profile = all_ec_profile_clean, continent_meta = combine_metatable)
species_sep <- list()
species_sep <- separate_samples_into_continents(input_profile = metaphlan_file_sel, continent_meta = combine_metatable)

saveRDS(species_sep, paste(output_dir, "species_sep.rds", sep = "/"))
saveRDS(all_ec_sep, paste(output_dir, "all_ec_sep.rds", sep = "/"))


########################################


#  random sampling                     # 


########################################

print("start random sampling")
random_sample <- function(number_pair = NA,sample_list = NA,unique_tag = "unique_tag"){
    sample_set_df <- NULL
    sample_set_df <- sample(sample_list, number_pair, replace = FALSE) %>% sort()
    sample_set_df <- as.matrix(sample_set_df)
    i <- 2
    while (i <= 100) { # this is repeat numbers 
        print(i)
        if (unique_tag == "no_unique"){
                print(str_c("no unique for total sample size:", number_pair, " samples"))
                random_set <- sample(sample_list, number_pair, replace = FALSE)
        }else{
        
                random_set <- sample(sample_list, number_pair, replace = FALSE) %>% sort()
        }
        random_set <- as.matrix(random_set)
        sample_set_df <- cbind(sample_set_df, random_set)
        sample_set_df <- unique(sample_set_df, MARGIN = 2)
        i <- ncol(sample_set_df)+1
    }
    return(sample_set_df)

}



## count enzyme number 

count_enzymes <- function(enzyme_table = NA, sample_set_target = NA) {

    enzyme_count <- matrix(nrow = ncol(sample_set_target))
    

    for (i in 1:ncol(sample_set_target)) {
        nsample = nrow(sample_set_target)
        enzyme_table_sel <- enzyme_table %>% select(sample_set_target[ ,i]) %>%
            filter(rowSums(.)>0)

        enzyme_count[i, 1] <- nrow(enzyme_table_sel)
       
    }



return(enzyme_count)

}

print("all ec random sampling")
# all ec list 
sample_set_allec <- list()
enzyme_matrix_for_plot_allec <- list()

for (name in names(all_ec_sep)){
    print(name)
    sample_list_allec <- colnames(all_ec_sep[[name]])[-1] %>% .[1:(length(.))]
    enzyme_matrix_for_plot_allec[[name]] <- matrix(nrow = 100, ncol = length(sample_list_allec)) # the row is the repeat times : 100 

    for (sample_size in 1:length(sample_list_allec)){
        print(str_c(name, sample_size, sep = ":"))
        sample_size_name <- paste(name, sample_size, sep = "_")
        if (sample_size == length(sample_list_allec)){
                tag <- "no_unique"
        }else{
                tag <- "unique_tag"
        }
        sample_set_allec[[sample_size_name]] <-  random_sample(sample_list = sample_list_allec, number_pair = sample_size, unique_tag = tag)
        enzyme_count_list <- count_enzymes(enzyme_table = all_ec_sep[[name]], sample_set_target = sample_set_allec[[sample_size_name]])
        enzyme_matrix_for_plot_allec[[name]][,sample_size] <- enzyme_count_list[,1]
    }
}



all_ec_matrix <- enzyme_matrix_for_plot_allec
saveRDS(all_ec_matrix, paste(output_dir, "/enzyme_matrix_for_plot_allec.rds", sep = "/"))

# sp sample list 
print("species random sampling")
sample_set_species <- list()
enzyme_matrix_for_plot_species <- list()

for (name in names(species_sep)){
    print(name)
    sample_list_species <- colnames(species_sep[[name]])[-1] %>% .[1:(length(.))]
    enzyme_matrix_for_plot_species[[name]] <- matrix(nrow = 100, ncol = length(sample_list_species)) # the row is the repeat times : 100 

    for (sample_size in 1:length(sample_list_species)){
        print(str_c(name, sample_size, sep = ":"))
        sample_size_name <- paste(name, sample_size, sep = "_")
        if (sample_size == length(sample_list_species)){
                tag <- "no_unique"
        }else{
                tag <- "unique_tag"
        }
        sample_set_species[[sample_size_name]] <-  random_sample(sample_list = sample_list_species, number_pair = sample_size, unique_tag = tag)
        enzyme_count_list <- count_enzymes(enzyme_table = species_sep[[name]], sample_set_target = sample_set_species[[sample_size_name]])
        enzyme_matrix_for_plot_species[[name]][,sample_size] <- enzyme_count_list[,1]
    }
}

species_matrix <- enzyme_matrix_for_plot_species
saveRDS(species_matrix, paste(output_dir, "/enzyme_matrix_for_plot_species.rds", sep = "/"))

print("finish all random sampling")


save.image(str_c(output_dir, "/othercategory_sampleaccumulation.RData"))


