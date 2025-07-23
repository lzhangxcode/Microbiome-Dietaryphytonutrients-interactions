# this script is used for running random sampling for the modified ec profile 
# author: lu.zhang 
# date: 01.29,2024 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Rscript dietary_ec_profile_random_sampling.R dietary_ec_profile_full dietary_ec_profile_simple combine_metatable output_dir", call.=FALSE)
} 

############################################################

# loading package                                          #

############################################################

print("Loading packages ... ...")
library(tidyverse)


############################################################

# input data                                               #

############################################################

print("Loading input files ... ...")
# per ec-compound pair per line 
modified_ec_profile_after_np <- args[1] %>% readRDS()
# modified ec profile - per ec per line
modified_ec_profile_simple <- args[2] %>% readRDS()
# combine metatable 
combine_metatable <- args[3] %>% read_tsv(., guess_max = 3000)
# output_dir 
output_dir <- args[4]

#############################################################

# separate ec profile into different continents             # 

#############################################################

ec_profile_sel <- list() 
for (con in combine_metatable$continent %>% unique()){
	print(con)
	target_sample <- combine_metatable %>% filter(continent %in% con)
	print(length(target_sample$internal_sample_id %>% unique()))
	ec_profile_sel[[con]] <- modified_ec_profile_simple %>% 
		select(c("ec_simplied", target_sample$internal_sample_id))
}


#############################################################

# random sampling                                           # 

#############################################################
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


sample_set <- list()
enzyme_matrix_for_plot <- list()

for (name in names(ec_profile_sel)){
    print(name)
    sample_list <- colnames(ec_profile_sel[[name]])[-1] %>% .[1:(length(.))]
    enzyme_matrix_for_plot[[name]] <- matrix(nrow = 100, ncol = length(sample_list)) # the row is the repeat times : 100 

    for (sample_size in 1:length(sample_list)){
        print(str_c(name, sample_size, sep = ":"))
	sample_size_name <- paste(name, sample_size, sep = "_")
        if (sample_size == length(sample_list)){
		tag <- "no_unique"
	}else{
		tag <- "unique_tag"
	}
        sample_set[[sample_size_name]] <-  random_sample(sample_list = sample_list, number_pair = sample_size, unique_tag = tag) # this part code has been checked. 
        enzyme_count_list <- count_enzymes(enzyme_table = ec_profile_sel[[name]], sample_set_target = sample_set[[sample_size_name]]) # this part code has been checked 
        enzyme_matrix_for_plot[[name]][,sample_size] <- enzyme_count_list[,1]
    }
}


######################################################################

#               save the matrix                                      # 

######################################################################

print("save the matrix") 
enzyme_melt_Af <- reshape2::melt(enzyme_matrix_for_plot[["Africa"]]) # Var1: repeat time, Var2 : sampling_size 

enzyme_melt_Am <- reshape2::melt(enzyme_matrix_for_plot[["America"]])

enzyme_melt_As <- reshape2::melt(enzyme_matrix_for_plot[["Asia"]]) 

enzyme_melt_Eu <- reshape2::melt(enzyme_matrix_for_plot[["Europe"]])

enzyme_melt_Oc <- reshape2::melt(enzyme_matrix_for_plot[["Oceania"]])

enzyme_melt_Af$Var1 <- "Africa"

enzyme_melt_Am$Var1 <- "America"

enzyme_melt_As$Var1 <- "Asia"

enzyme_melt_Eu$Var1 <- "Europe"

enzyme_melt_Oc$Var1 <- "Oceania"


all_continents_melt <- rbind(enzyme_melt_Af, enzyme_melt_Am, enzyme_melt_As, enzyme_melt_Eu, enzyme_melt_Oc)
all_continents_melt$Var2 <- as.factor(all_continents_melt$Var2)

saveRDS(all_continents_melt, paste(output_dir,"/all_continents_modified_melt.rds", sep = ""))

all_continents_melt <- rbind(enzyme_melt_Af, enzyme_melt_Am, enzyme_melt_As, enzyme_melt_Eu, enzyme_melt_Oc)
all_continents_melt$Var2 <- as.factor(all_continents_melt$Var2)

saveRDS(all_continents_melt, paste(output_dir,"/all_continents_modified_melt_np.rds", sep = ""))

#####################################################################

#               save image                                          # 

##################################################################### 

print("save the rdata... ...")
save.image(paste(output_dir, "step10_rarefaction_curve_modified_ec_profile_manual_check_plot.RData",sep = "/"))




