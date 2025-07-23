# the aim of the script is to do 1% prevalence filtering per continent & draw upset plot for foods
# author: lu.zhang 
# date: 2024.mai.31

# this script run in terminal 

######################################################################

# Data input 

######################################################################

library(tidyverse)
# modified ec profile 
modified_ec_profile_after_np <- readRDS(args[1])
modified_ec_profile_after_np <- readRDS("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result1_rarefaction_curve/modified_ec_profile/modified_ec_profile_after_np.rds") # this file, column ec_simplied Compounds are the needed ones.

# full ec database :it is not used in the compounds figure 
full_ec_db <- read_tsv(args[2])
full_ec_db <- read_tsv("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result12_beneficial_food/nutrichem_md_manual_with_non_microbiotaEC_fullEC.tsv") 

# metatable
combine_metadata <- read_tsv(args[3], guess_max = 3000)
combine_metadata <- read_tsv("/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result11_the_smallest_species_the_mostEC/compounds/combine_metadata_sel.tsv", guess_max = 3000)

# output 
output_dir <- args[4]
output_dir <- "/sbidata/projects/lzhang/202012_small_molecule/Analysis/metaphlan_human3/Result/summary_result/Manual_check/Result8_related_smallmolecules/1_percent/food/"



combine_metadata_Af <- combine_metadata %>%
    dplyr::filter(., continent == "Africa")
combine_metadata_Am <- combine_metadata %>%
    dplyr::filter(., continent == "America")
combine_metadata_As <- combine_metadata %>%
    dplyr::filter(., continent == "Asia")
combine_metadata_Eu <- combine_metadata %>%
    dplyr::filter(., continent == "Europe")
combine_metadata_Oc <- combine_metadata %>%
    dplyr::filter(., continent == "Oceania")



######################################################################

# check again 5% prevalence + 1% prevalence per continent 

######################################################################

prevalence_filtering_percontinent <- function(input_modified_ec = modified_ec_profile_after_np, metatable = combine_metadata, prevalence = 0.01){

	# continent number 
	continent_number <- data.frame(table(metatable$continent))
	# filter each continent 
	input_modified_ec_long <- input_modified_ec %>%
    		dplyr::select(-Compounds,-Tag) %>%
    		unique() %>%
    		pivot_longer(!ec_simplied) %>%
    		add_column(num = ifelse(.$value>0,1,0)) %>%
		left_join(., metatable, by = c("name" = "internal_sample_id")) %>% 
    		group_by(continent, ec_simplied) %>%
    		summarize(sum = sum(num)) %>%
		left_join(., continent_number, by = c("continent" = "Var1")) %>%
    		mutate(freq = sum / Freq) %>%
    		ungroup()
	
	# filter ec > prevalence 
        keep_ec <- input_modified_ec_long %>%
    		filter(freq >= prevalence)

	
	# produce 
	keep_ec_list <- list()
	keep_ec_list[["af"]] <- keep_ec %>% filter(continent == "Africa")
	keep_ec_list[["am"]] <- keep_ec %>% filter(continent == "America")
	keep_ec_list[["as"]] <- keep_ec %>% filter(continent == "Asia")
	keep_ec_list[["oc"]] <- keep_ec %>% filter(continent == "Oceania")
	keep_ec_list[["eu"]] <- keep_ec %>% filter(continent == "Europe")

	# return
	return(keep_ec_list)

}



# 1% prevalence per continent 

keep_ec_list_for_continent <- prevalence_filtering_percontinent(input_modified_ec = modified_ec_profile_after_np, prevalence = 0.01)


######################################################################

# filter each continent ec  

######################################################################
# 1% prevalence per continent 

# africa
ec_profile_Af_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Af$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["af"]]$ec_simplied)
ec_profile_Af_con <- ec_profile_Af_con[rowSums(ec_profile_Af_con[-(1:2)])>0,]

# america
ec_profile_Am_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Am$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["am"]]$ec_simplied)
ec_profile_Am_con <- ec_profile_Am_con[rowSums(ec_profile_Am_con[-(1:2)])>0,]

# asia
ec_profile_As_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_As$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["as"]]$ec_simplied)
ec_profile_As_con <- ec_profile_As_con[rowSums(ec_profile_As_con[-(1:2)])>0,]

# europe
ec_profile_Eu_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Eu$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["eu"]]$ec_simplied)
ec_profile_Eu_con <- ec_profile_Eu_con[rowSums(ec_profile_Eu_con[-(1:2)])>0,]

# oceania
ec_profile_Oc_con <- modified_ec_profile_after_np %>%
    dplyr::select("ec_simplied", "Compounds", combine_metadata_Oc$internal_sample_id) %>%
    filter(ec_simplied %in% keep_ec_list_for_continent[["oc"]]$ec_simplied)
ec_profile_Oc_con <- ec_profile_Oc_con[rowSums(ec_profile_Oc_con[-(1:2)])>0,]


######################################################################

# acquire the food list - percontinent  

######################################################################

input_profile_ec <- ec_profile_Oc_con 
nutrichem_full_ec_db <- full_ec_db

combine_compounds_with_food <- function(input_profile_ec = NA, nutrichem_full_ec_db = full_ec_db){

	input_profile_ec_compounds_foods <- input_profile_ec %>% 
		left_join(., nutrichem_full_ec_db %>% select(CompoundID, PlantID), by = c("Compounds" = "CompoundID")) # here I need to separate the TAXID 
	input_profile_ec_compounds_foods_sep <- input_profile_ec_compounds_foods %>% 
		separate_rows(., PlantID, sep = ";")
       return(input_profile_ec_compounds_foods_sep)	
}


ec_profile_Af_con_food <- combine_compounds_with_food(input_profile_ec = ec_profile_Af_con)
ec_profile_Am_con_food <- combine_compounds_with_food(input_profile_ec = ec_profile_Am_con)
ec_profile_As_con_food <- combine_compounds_with_food(input_profile_ec = ec_profile_As_con)
ec_profile_Eu_con_food <- combine_compounds_with_food(input_profile_ec = ec_profile_Eu_con)
ec_profile_Oc_con_food <- combine_compounds_with_food(input_profile_ec = ec_profile_Oc_con)

# 1 percent 
listInput_con <- list(AF = ec_profile_Af_con_food$PlantID %>% unique(),
                  AM = ec_profile_Am_con_food$PlantID %>% unique(),
                  AS = ec_profile_As_con_food$PlantID %>% unique(),
                  EU = ec_profile_Eu_con_food$PlantID %>% unique(),
                  OC = ec_profile_Oc_con_food$PlantID %>% unique())

# calculate the compounds number (for the order in the upset plot)

africa_food_num_con <- length(ec_profile_Af_con_food$PlantID %>% unique())
america_food_num_con <- length(ec_profile_Am_con_food$PlantID %>% unique())
asia_food_num_con <- length(ec_profile_As_con_food$PlantID %>% unique())
europe_food_num_con <- length(ec_profile_Eu_con_food$PlantID %>% unique())
oceania_food_num_con <- length(ec_profile_Oc_con_food$PlantID %>% unique())


# draw plot 
library(UpSetR)

color_food_tibble_con <- tibble(continent = c("AF", "AM", "AS", "EU", "OC"),
                                        color = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF"),
                                        number = c(africa_food_num_con, america_food_num_con, asia_food_num_con, europe_food_num_con, oceania_food_num_con)) %>%
                arrange(desc(number))



pdf(paste(output_dir, "upset_food_1percent_con.pdf", sep = "/"))#, width = 6, height = 5)
text_scale_options3 <- c(2.5, 2, 2, 1.6, 2.5, 3)
foods_upset <- UpSetR::upset(fromList(listInput_con), sets = c("AF", "AM", "AS", "EU", "OC"),
                                   sets.x.label = "LMWMs",
                                   mainbar.y.label = "N of foods biotransf. \nby microbiome",
                                   mb.ratio = c(0.65, 0.35), point.size = 4,
                                   order.by = "freq", #empty.intersections = "on",
                                   text.scale = text_scale_options3, sets.bar.color = color_food_tibble_con$color)
print(foods_upset)
dev.off()


save.image(str_c(output_dir, "/food_number_1_percent.RData"))


# save the list into rds
saveRDS(keep_ec_list_for_continent, str_c(output_dir, "/keep_ec_list_for_each_continents_food.rds"))


